#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "gsa/gsacak.h"
#include "utils.h" 

// Compute the SA of a parsing and its associated BWT. 
// Remap the chars in the .last file according to the
// BWT/SA permutation producing the .bwlast file 

// From the BWT compute and output the .ilist file giving 
// for each parsing symbol (considered in alphabetical order)
// the list of BWT positions where that symbol occurs.
// Among the symbols the EOF symbol is included so the 
// the list has size |BWT| = |parse| + 1. 
// Ilist is used for the generation of the final BWT
// together with the .occ file giving for each symbol its number of 
// occurrences in the BWT. 
// Once the list is computed the BWT is no longer useful and is discarded

// The parsing symbols are assumed to be 32 bit uints
// Currently, the parsing is assumed to be of length at most 2^32-1
// so each ilist entry has size 4 bytes. 

// Note: to support in the future parsing longer than 2^32-2, the arrays
// occ[] F[] and Ilist[] should become of type sa_index_t. However, this
// would also require to store to file more than 4 bytes per entry for
// occ[] (file .occ) and Ilist[] (file .ilist) 


// -------------------------------------------------------
// type used to represent an entry in the SA
// this is currently 32 bit for gsacak and 64 bit for gsacak-64
// note that here we use sacak (SA computation for a single string of 32 bit symbols) 
typedef uint_t sa_index_t;


// read the parse file 
uint32_t *read_parse(char *basename,long *tsize) 
{  
  FILE *parse = open_aux_file(basename,"parse","rb");
  // get file size
  if(fseek(parse,0,SEEK_END)!=0) die("parse fseek");
  long nn = ftell(parse);
  // check input file is OK
  if(nn%4!=0) {
    printf("Invalid input file: size not multiple of 4\n");
    exit(1);
  }
  #if !M64
  // if in 32 bit mode, the number of words is at most 2^31-2
  if(nn/4 > 0x7FFFFFFE) {
    printf("Input containing more than 2^31-2 phrases!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
  #else
  // if in 64 bit mode, the number of words is at most 2^32-2 for now
  if(nn/4 > 0xFFFFFFFEu) {
    printf("Input containing more than 2^32-2 phrases!\n");
    printf("This is currently a hard limit\n");
    exit(1);
  }
  #endif
  printf("Parse file contains %ld words\n",nn/4);
  long n = nn/4;
  // ------ allocate and read text file
  uint32_t *Text = malloc((n+1)*sizeof(*Text));
  if(Text==NULL) die("malloc failed (Text)");
  rewind(parse);

  // read the array in one shot
  assert(sizeof(*Text)==4);
  size_t s = fread(Text,sizeof(*Text),n,parse);
  if(s!=n) {
    char *msg=NULL;
    int e= asprintf(&msg,"read parse error: %zu vs %ld\n", s,n); 
    (void) e; die(msg);
  }
  
  if(fclose(parse)!=0) die("parse file close");
  Text[n]=0; // sacak needs a 0 eos 
  *tsize= n;
  return Text;
}  



int main(int argc, char *argv[])
{
  uint32_t *Text; // array of parsing symbols 
  long n;         // length of Text[] (not including final 0 symbol)
  size_t s;

  time_t start_wc = time(NULL);
  // check input data
  if(argc<2){
    printf("Usage: %s basename\n\n", argv[0]);
    puts("Compute the BWT of basename.parse and store its inverted list occurrence");
    puts("Permutes the file basename.last according to the same permutation");
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  // read parse file
  Text = read_parse(argv[1],&n);
  
  // open .last file for reading and .bwlast for writing
  FILE *lastin = open_aux_file(argv[1],"last","rb");
  FILE *lastout = open_aux_file(argv[1],"bwlast","wb");

  // ------- compute alphabet size-1
  uint32_t k=0;
  for(long i=0;i<n;i++) {
    if(Text[i]>k) k = Text[i];
  }
  // if(k>=0x7FFFFFFE) die("Emergency exit! Alphabet larger than 2^31 (2)");
  // -------- alloc and compute SA
  sa_index_t *SA = malloc((n+1)*sizeof(*SA));
  if(SA==NULL) die("malloc failed  (SA)");
  printf("Computing SA of size %ld over an alphabet of size %u\n",n+1,k+1);
  int depth = sacak_int(Text, SA, n+1, k+1);
  if(depth>=0)
    printf("SA computed with depth: %d\n", depth);
  else
    die("Error computing the SA"); 
    
  // allocate and load the last array
  uint8_t *last = malloc(n);
  if(last==NULL) die("malloc failed (LAST)"); 
  s = fread(last,1,n,lastin);
  if(s!=n) die("last read");
  if(fclose(lastin)!=0) die("last file close");
  
  // transform SA->BWT inplace and write remapped last array
  sa_index_t *BWTsa = SA;
  assert(n>1);
  assert(SA[0]==n);
  BWTsa[0] = Text[n-1];
  if(fputc(last[n-2],lastout)==EOF) die("bwlast output 1");
  for(long i=1;i<=n;i++) {
    if(SA[i]==0) {
      assert(i==1);
      BWTsa[i] = 0;       // eof in BWT
      if(fputc(0,lastout)==EOF) die("bwlast output 2"); // dummy char 
    }
    else {
      if(SA[i]==1) {if(fputc(last[n-1],lastout)==EOF) die("bwlast output 3");}
      else    {if(fputc(last[SA[i]-2],lastout)==EOF) die("bwlast output 4");}
      BWTsa[i] = Text[SA[i]-1];
    }
  }
  if(fclose(lastout)!=0) die("bwlast close");
  printf("---- %ld bwlast chars written ----\n",n+1);
  free(last);
  
  // --- copy BWT to text array (symbol by symbol since sizeof could be different)
  uint32_t *BWT = Text;
  for(long i=0;i<=n;i++)
    BWT[i] = BWTsa[i];
 
  // read # of occ of each char from file .occ
  uint32_t *occ = malloc((k+1)*sizeof(*occ)); // extra space for the only occ of 0
  if(occ==NULL) die("malloc failed (OCC)");
  FILE *occin = open_aux_file(argv[1],"occ","rb");
  s = fread(occ+1,sizeof(*occ), k,occin);
  if(s!=k) die("not enough occ data!");
  occ[0] = 1; // we know there is somewhere a 0 BWT entry 
  fclose(occin);
  // create F vector
  uint32_t *F = malloc((k+1)*sizeof(*F));
  if(F==NULL) die("malloc failed (F)");
  // init F[] using occ[]
  F[0] = 0;
  for(int i=1;i<=k;i++)
    F[i]=F[i-1]+occ[i-1];
  assert(F[k]+occ[k]==n+1);
  puts("---- computing inverted list ----");
  // ----- compute inverse list overwriting SA
  uint32_t *IList = (uint32_t *) SA;
  for(long i=0;i<=n;i++) {
    IList[F[BWT[i]]++] = i;
    occ[BWT[i]]--;
  }
  // ---check
  assert(IList[0]==1); // EOF is in BWT[1] since P[0] = $xxx is the smallest word and appears once
  assert(BWT[IList[0]]==0);
  for(long i=0;i<=k;i++) 
    assert(occ[i]==0);
  // ---save Ilist   
  FILE *ilist = open_aux_file(argv[1],"ilist","wb");
  s = fwrite(IList,sizeof(*IList),n+1,ilist);
  if(s!=n+1) die("Ilist write");
  fclose(ilist);
  printf("---- %ld ilist positions written (%ld bytes) ----\n",n+1,(n+1)*4l);
  // deallocate
  free(F);
  free(occ);
  free(SA);
  free(Text);
  printf("==== Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));  
  return 0;
}
