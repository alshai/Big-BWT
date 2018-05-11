#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "gsacak.h"


// write message and exit
void die(const char *s)
{
  perror(s);
  exit(1);
}    


int main(int argc, char *argv[]){

  // assert(sizeof(int_t)==4);
  uint_t *Text;
  long n=0;
  int e;
  size_t s;

  // input data
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

  // read main file
  char *name;
  e = asprintf(&name,"%s.parse",argv[1]);
  if(e<1) die("asprint error");
  FILE *parse = fopen(name,"rb");
  free(name);
  if(parse==NULL) die("parse open");
  // get file size
  if(fseek(parse,0,SEEK_END)!=0) die("fseek");
  long nn = ftell(parse);
  // check input file is OK
  if(nn%4!=0) {
    printf("Invalid input file: size not multiple of %d\n",4);
    exit(1);
  }
  #if !M64
  // if in 32 bit mode, the number of words is at most 2^31
  if(nn/4 > 0x7FFFFFFE) {
    printf("Input containing more than 2^31-2 phrases!\n");
    printf("Please use %s64\n",argv[0]);
    exit(1);
  }
  #endif
  printf("Parse file contains %ld words\n",nn/4);
  n = nn/4;
  // open last file for reading and bwlast for writing
  e = asprintf(&name,"%s.last",argv[1]);
  if(e<1) die("asprint error");  
  FILE *lastin = fopen(name,"rb");
  free(name);
  if(lastin==NULL) die("unable to open .last file");
  e = asprintf(&name,"%s.bwlast",argv[1]);
  if(e<1) die("asprint error");  
  FILE *lastout = fopen(name,"wb");
  free(name);
  if(lastout==NULL) die("unable to open .bwlast file");  
  // ------ allocate and read text file
  Text = malloc((n+1)*sizeof(*Text));
  if(Text==NULL) die("malloc failed  (Text)");
  rewind(parse);
  #if M64
  // copy 32bit int's into a 64 bit array
  for(long i=0;i<n;i++) {
    uint32_t b;
    s = fread(&b,4,1,parse);
    if(s!=1) die("Error reading from parse file");
    Text[i]=b;
  }
  #else
  // read the array in one shot
  assert(sizeof(*Text)==4);
  s = fread(Text,sizeof(*Text),n,parse);
  if(s!=n) {printf("%zu vs %ld\n", s,n); die("read parse error");}
  #endif
  fclose(parse);
  Text[n]=0; // sacak needs a 0 eos 
  // ------- compute alphabet size-1
  uint32_t k=0;
  for(int i=0;i<n;i++)
    if(Text[i]>k) k = Text[i];
  if(k>0x7FFFFFFE) die("Alphabet larger than 2^31");
  // -------- alloc and compute SA
  uint_t *BWT = malloc((n+1)*sizeof(*BWT));
  if(BWT==NULL) die("malloc failed  (BWT)");
  uint_t *SA=BWT; 
  printf("Computing SA of size %ld over an alphabet of size %u\n",n+1,k+1);
  int depth = sacak_int((int_t *) Text, SA, n+1, k+1);
  if(depth>=0)
    printf("SA computed with depth: %d\n", depth);
  else
    die("Error computing the SA"); 
  // allocate and load the last array
  uint8_t *last = malloc(n);
  if(last==NULL) die("malloc failed (LAST)"); 
  s = fread(last,1,n,lastin);
  if(s!=n) die("last read");
  fclose(lastin);
  // transform SA->BWT and write remapped last array
  assert(n>1);
  assert(SA[0]==n);
  BWT[0] = Text[n-1];
  fputc(last[n-2],lastout);
  for(long i=1;i<=n;i++) {
    if(SA[i]==0) {
      assert(i==1);
      BWT[i] = 0;       // eof in BWT
      if(fputc(0,lastout)==EOF) die("bwlast output"); // dummy char 
    }
    else {
      if(SA[i]==1) putc(last[n-1],lastout);
      else    fputc(last[SA[i]-2],lastout);
      BWT[i] = Text[SA[i]-1];
    }
  }
  fclose(lastout);
  printf("---- %ld bwlast chars written ----\n",n+1);
  #if M64
  free(Text);
  #endif
  // read # of occ of each char from file .occ
  uint32_t *occ = malloc((k+1)*4); // extra space for the only occ of 0
  if(occ==NULL) die("malloc failed (OCC)");
  e = asprintf(&name,"%s.occ",argv[1]);
  if(e<1) die("asprint error");  
  FILE *occin = fopen(name,"rb");
  free(name);
  if(occin==NULL) die("Error opening occ file");
  s = fread(occ+1,4, k,occin);
  if(s!=k) die("not enough occ data!");
  occ[0] = 1; // we know there is somewhere a 0 BWT entry 
  fclose(occin);
  // create F vector
  uint32_t *F = malloc((k+1)*4);
  if(F==NULL) die("malloc failed (F)");
  F[0] = 0;
  for(int i=1;i<=k;i++)
    F[i]=F[i-1]+occ[i-1];
  assert(F[k]+occ[k]==n+1);
  puts("---- computing inverted list ----");
  // ----- compute inverse list
  #if M64
  uint32_t *IList = malloc((n+1)*sizeof(*IList));
  if(IList==NULL) die("malloc failed (IList)");
  #else
  uint32_t *IList = Text;
  #endif
  for(long i=0;i<=n;i++) {
    IList[F[BWT[i]]++] = i;
    occ[BWT[i]]--;
  }
  // ---check
  assert(BWT[IList[0]]==0);
  for(long i=0;i<=k;i++) assert(occ[i]==0);
  // ---save Ilist 
  e = asprintf(&name,"%s.ilist",argv[1]);
  if(e<1) die("asprint error");  
  FILE *ilist = fopen(name,"wb");
  free(name);
  if(ilist==NULL) die("Unable to open .ilist file");
  s = fwrite(IList,sizeof(*IList),n+1,ilist);
  if(s!=n+1) die("Ilist write");
  fclose(ilist);
  printf("---- %ld iflist positions written (%ld bytes) ----\n",n+1,(n+1)*4l);
  // deallocate
  free(F);
  free(occ);
  free(last);
  free(BWT);
  free(IList);
  return 0;
}

