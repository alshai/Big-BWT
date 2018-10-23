#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include "gsa/gsacak.h"
#include "utils.h" 

// Compute the SA of a parsing and its associated BWT. 
// Remap the chars in the .last file according to the
// BWT/SA permutation producing the .bwlast file 

// If the -s option is used, also remap the sa values in the .sai 
// file (that uses IBYTES bytes per entry) producing the .bwsai file

// The remapping is done as follows, let T[0] T[1] ... T[n-1] T[n] denote
// the input parse where T[n]=0 is an extra EOS symbol added here   
// For j=0,...,n, if BWT[j] = T[i] then:
//   bwisa[j]  is the ending position+1 of T[i] in the original text 
//   bwlast[j] is the char immediately before T[i] in the original text,
//             ie the char in position w from the end of T[i-1]
// Since T[n]=0, and T[0] = $abc... is the only phrase starting with $, 
// it is T[n]<T[0]<T[i] and BWT[0]=T[n-1], BWT[1]=T[n]=0. We have 
//   bwisa[0]  = ending position + 1 of T[n-1], 
//   bwlast[0] = char in position w from the end of T[n-2]
//   bwisa[1] and bwlast[1] are dummy values since T[n] is not a real phrase
// For some j it is BWT[j]=T[0], for that j we set
//   bwisa[j] = ending position + 1 of T[0] as expected
//   bwlast[j] = char in position w from the end in T[n-1] (it should formally be
//               T[n] but again that is a dummy symbol and we skip it)

// From the BWT we compute and output the .ilist file giving 
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

// -------------------------------------------------------------
// struct containing command line parameters and other globals
typedef struct {
   char *basename;
   bool SAinfo;
   int th;    // number of segments for the last and sa files
} Args;


// read the parse file, add a 0 EOS symbol and return a pointer 
// to a new allocate uint32_t array containing it. Store size in *tsize
// note *tsize is the number of element in the parsing, but we add a 0
// symbol at the end so the returned array has *tsize+1 elements 
static uint32_t *read_parse(char *basename, long *tsize) 
{  
  FILE *parse = open_aux_file(basename,EXTPARSE,"rb");
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
  // if in 64 bit mode, the number of words is at most 2^32-2 (for now)
  if(nn/4 > 0xFFFFFFFEu) {
    printf("Input containing more than 2^32-2 phrases!\n");
    printf("This is currently a hard limit\n");
    exit(1);
  }
  #endif
  printf("Parse file contains %ld words\n",nn/4);
  long n = nn/4;
  // ------ allocate and read text file, len is n+1 for the EOS 
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

static void print_help(char *name)
{
  printf("Usage: %s <basename> [options]\n\n", name);
  puts("Compute the BWT of basename.parse and store its inverted list occurrence");
  puts("Permute the file basename.last according to the same permutation");
  puts("  Options:");
  puts("\t-h  \tshow help and exit");
  puts("\t-s  \tpermute also sa info");
  exit(1);
}

static void parseArgs(int argc, char** argv, Args *arg ) {
  extern int optind, opterr, optopt;
  extern char *optarg;  
  int c;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("\n");

  arg->SAinfo = false;
  arg->th = 0;
  while ((c = getopt( argc, argv, "sht:") ) != -1) {
    switch(c) {
      case 's':
      arg->SAinfo = true; break;
      case 'h':
         print_help(argv[0]); exit(1);
      case 't':
         arg->th = atoi(optarg); break;
      case '?':
      puts("Unknown option. Use -h for help.");
      exit(1);
    }
  }
  // read base name as the only non-option parameter 
  if (argc!=optind+1)
    print_help(argv[0]);
  arg->basename = argv[optind];
}

static sa_index_t *compute_SA(uint32_t *Text, long n, long k) 
{
  sa_index_t *SA = malloc(n*sizeof(*SA));
  if(SA==NULL) die("malloc failed  (SA)");
  printf("Computing SA of size %ld over an alphabet of size %ld\n",n,k);
  int depth = sacak_int(Text, SA, n, k);
  if(depth>=0)
    printf("SA computed with depth: %d\n", depth);
  else
    die("Error computing the SA"); 
  return SA;
}

// load the last array produced by newscan into an array
static uint8_t *load_last(Args *arg, long n)
{  
  // open .last file for reading
  mFile *lastin = mopen_aux_file(arg->basename,EXTLST,arg->th);
  // allocate and load the last array
  uint8_t *last = malloc(n);
  if(last==NULL) die("malloc failed (LAST)"); 
  size_t s = mfread(last,1,n,lastin);
  if(s!=n) die("last read");
  if(mfclose(lastin)!=0) die("last file close");
  return last;
}

  
static uint8_t *load_sa_info(Args *arg, long n)
{  
  // maybe sa info was not required 
  if(arg->SAinfo==false) return NULL;
  // open .sa_info file for reading and .bwlast for writing
  mFile *fin = mopen_aux_file(arg->basename,EXTSAI,arg->th);
  // allocate and load the sa info array
  uint8_t *sai = malloc(n*IBYTES);
  if(sai==NULL) die("malloc failed (SA INFO)"); 
  size_t s = mfread(sai,IBYTES,n,fin);
  if(s!=n) die("sa info read");
  if(mfclose(fin)!=0) die("sa info file close");
  return sai;
}

static FILE *open_sa_out(Args *arg)
{
  if(arg->SAinfo==false) return NULL;
  return open_aux_file(arg->basename,EXTBWSAI,"wb"); 
}


int main(int argc, char *argv[])
{
  uint32_t *Text; // array of parsing symbols 
  long n;         // length of Text[] (not including final 0 symbol)
  size_t s;
  Args arg;

  // read arguments 
  parseArgs(argc,argv,&arg);
  // start measuring wall clock time 
  time_t start_wc = time(NULL);
  // read parse file
  Text = read_parse(arg.basename,&n);
  
  // ------- compute largest input symbol (ie alphabet size-1)
  long k=0;
  for(long i=0;i<n;i++) {
    if(Text[i]>k) k = Text[i];
  }
  // -------- alloc and compute SA of the parse
  sa_index_t *SA = compute_SA(Text,n+1,k+1);

  // load last file 
  uint8_t *last = load_last(&arg,n);
  // load sa info file, if requested
  uint8_t *sa_info = load_sa_info(&arg,n); 
  FILE *lastout = open_aux_file(arg.basename,EXTBWLST,"wb");   
  FILE *sa_out = open_sa_out(&arg);
  // note that lastout and sa_out files will have n+1 elements instead of n

  // transform SA->BWT inplace and write remapped last array, and possibly sainfo
  sa_index_t *BWTsa = SA; // BWT overlapping SA
  assert(n>1);
  // first BWT symbol
  assert(SA[0]==n);
  BWTsa[0] = Text[n-1];
  if(fputc(last[n-2],lastout)==EOF) die("bwlast output 1");
  if(arg.SAinfo) get_and_write_myint(sa_info,n,n-1,sa_out); // ending position+1 of Text[n-1] in original text T (is |T|+w) 
  // 2nd, 3rd etc BWT symbols 
  for(long i=1;i<=n;i++) {
    if(SA[i]==0) {  
      assert(i==1);  // Text[0]=$abc... is the second lex word 
      BWTsa[i] = 0;   // eos in BWT, there is no phrase in D corresponding to this symbol so we write dummy values
      if(fputc(0,lastout)==EOF) die("bwlast output 2"); // dummy char 
      if(arg.SAinfo) write_myint(0,sa_out); // dummy end of word position, it is never used an 0 does not appear elsewhere in sa_out
    }
    else {
      if(SA[i]==1) {
        // BWT[i] = Text[0] = $abcd... = first word in the parsing where $ now plays the role of the EOS in the original text  
        if(fputc(last[n-1],lastout)==EOF) die("bwlast output 3");
      }
      else    {if(fputc(last[SA[i]-2],lastout)==EOF) die("bwlast output 4");}
      if(arg.SAinfo) get_and_write_myint(sa_info,n,SA[i]-1,sa_out); // ending position of BWT symbol in original text
      BWTsa[i] = Text[SA[i]-1];
    }
  }
  if(fclose(lastout)!=0) die("bwlast close");
  printf("---- %ld bwlast chars written ----\n",n+1);
  free(last);
  if(arg.SAinfo) {
    if(fclose(sa_out)!=0) die("sa_out close");
    free(sa_info);
  } 
  
  // --- copy BWT to text array (symbol by symbol since sizeof could be different)
  uint32_t *BWT = Text;
  for(long i=0;i<=n;i++)
    BWT[i] = BWTsa[i];
 
  // read # of occ of each char from file .occ
  uint32_t *occ = malloc((k+1)*sizeof(*occ)); // extra space for the only occ of 0
  if(occ==NULL) die("malloc failed (OCC)");
  FILE *occin = open_aux_file(arg.basename,"occ","rb");
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
  FILE *ilist = open_aux_file(arg.basename,EXTILIST,"wb");
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
