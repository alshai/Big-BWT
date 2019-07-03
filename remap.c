/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   remap.c
   simple reporting/remapping tool for suffix array construction
    
   Report only usage:
   
      remap filename [-v]
   
   Report number of distinct symbols and entropy of input file
   With option [-v] report number of occurrences of each symbol
   
   Remap usage:
     
      remap filename k [-v]
      
   Try to remap the symbols in the input file so that the symbols 0..k
   do not appear in the remapped file. The remapping is done preserving the
   lexicographic order of the suffixes; the remapped content is written 
   to file filename.k
   
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

void pretty_putchar(int c);
double report_alpha(unsigned long c[], double len);
void remap_alpha(unsigned long count[], int k, int map[]);
int verbose = 0; // verbose output 

void usage(char *name) {
    fprintf(stderr, "Usage:\n\t %s filename [k] [-v] [-h]\n\n",name);
    fprintf(stderr, "If called with only a file name reports the number of occurrences\n");    
    fprintf(stderr, "of each symbol and the order-0 entropy of the file\n\n");
    fprintf(stderr, "If called with a second integer argument remap symbols so that\n");
    fprintf(stderr, "symbols 0,1,..k, do not appear in the text (if possible)\n"); 
    fprintf(stderr, "The remapping is done preserving the lexicographic order of substrings\n\n"); 
    fprintf(stderr, "Command line options:\n"); 
    fprintf(stderr, "  -h\tthis help message\n"); 
    fprintf(stderr, "  -v\treport number of occurrences of each symbol\n"); 
    exit(1);        
}


int main(int argc, char *argv[])
{
  //extern char *optarg;
  extern int optind, opterr, optopt;
  char *fnam=NULL;
  FILE *f;
  int c,k=-1;

  /* ------------- read options from command line ----------- */
  opterr = 0;
  while ((c=getopt(argc, argv, "hv")) != -1) {
    switch (c) {
      case 'h':
        usage(argv[0]); break;
      case 'v':
        verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c. Try -h for help\n", optopt);
        exit(1);
    }
  }
  
  if(optind==argc) 
    usage(argv[0]);
  else {
    fnam=argv[optind];
    if(optind+1<argc) {
      k = atoi(argv[optind+1]);
      if(k<0) {
        fprintf(stderr,"Remap parameter must be non negative");
        exit(2);
      }
    }
  }

  /* -------- open and read input file ------------- */
  if (! (f=fopen(fnam, "rb"))) {
    perror(fnam);
    return 1;
  }
  unsigned long n=0, count[256] = {0};
  while(1) {
    int c = getc(f);
    if(c==EOF) break;
    else count[c]++;
    n++;
  }
  if(n==0) {puts("Empty input file"); return 2;}
  // --- count # symbols 
  int tot=0, i=0;
  for(i=0;i<256;i++)
    if(count[i]!=0) tot++;
  // report alphabet
  double entropy = report_alpha(count,n);
  printf("Number of input symbols: %li\n",n);  
  printf("Number of distinct symbols: %d\n",tot);  
  printf("Entropy: %lf\n",entropy/n); 

  // --- remap symbols if requested ---
  if(k>=0) {
    if(tot+k>=256) {
      printf("Cannot remap symbols as requested\n");
      return(3);
    }
    {// check if remapping is really necessary
      unsigned int firstk = 0;
      for(i=0;i<=k;i++) firstk += count[i];
      if(firstk==0) {
        if(k>0)
          printf("Symbols 0..%d not present in input file. Remapping not necessary\n",k);
        else
          printf("Symbol 0 not present in input file. Remapping not necessary\n");
        return 0;
      }
    }
    int map[256];
    remap_alpha(count, k, map);
    char filename[FILENAME_MAX];
    snprintf(filename,FILENAME_MAX,"%s.%d",fnam,k);
    printf("Writing remapped file to %s\n",filename);
    FILE *g= fopen(filename,"wb");
    rewind(f);
    unsigned long gcount[256] = {0};
    while(1) {
      int c = getc(f);
      if(c==EOF) break;
      assert(map[c]>=0);
      assert(map[c]<256);
      gcount[map[c]]++;
      putc(map[c],g);
    }
    fclose(g);
    entropy = report_alpha(gcount,n);
    printf("Number of distinct symbols in remapped file: %d\n",tot);  
    printf("Entropy times length: %lf\n",entropy/n); 
  }
  fclose(f);
  return 0;
}


// remap alphabet in order preserving way making sure
// that no symbol is remapped to 0..k 
void remap_alpha(unsigned long count[], int k, int map[])
{
  int jump = k+1, i=0;
  for(i=0;i<256;i++) {
    if(count[i]>0) {
      map[i] = i + jump;
      if(jump>0) 
        printf("Remapping %3d -> %3d\n",i,map[i]); 
    }
    else {
      map[i] = -1;         // i not in input file 
      if(jump>0) jump--;   // if not already zero reduce jump
    }
  }
  assert(jump==0);
}

// report the number of occurrences of each symbol and return the 
// order zero entropy times length 
double report_alpha(unsigned long count[], double len)
{
  double entropia=0;
  // --- print
  int c=0;
  for(c=0; c<256; c++)
    if(count[c]) {
      if(verbose>0) {
        printf("%3d  %2x  ", c, c);
        if(c>31 && c<127)
          putchar(c);
        else
          putchar(' ');
        printf("  %20ld\n",count[c]);
      }
      entropia += count[c]*log2(len/count[c]);
    }
  return entropia;
}


// ----- this function prints any char in a readable form
void pretty_putchar(int c)
{
  
  if(c>=32 && c<127)      // printable char
    printf("  %c", c);
  else if(c=='\n')
    printf(" \\n");        // \n
  else if(c=='\t')
    printf(" \\t");        // \t
  else     
    printf(" %02x", c);      // print hex code
}
