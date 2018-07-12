/* ******************************************************************
 * simplebwt
 * Given a text T compute the its bwt via the suffix array using the 
 * SACA-K algorithm.
 * The input file cannot contain the 0x0 char that is used int the BWT as 
 * the EOF symbol. The output file has extension .Bwt and length |T|+1
 * 
 * If compiled with M64=1 uses 64 bit uints for the suffix array; 
 * the overall space is 9n bytes and the input can be as large as 2^63-1
 * Otherwise uses 32 bit uints for the SA: the overall space is 5n bytes 
 * but the input can be of length at most 2^31 -1
 * ****************************************************************** */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "gsa/gsacak.h"


void die(const char *s)
{
  perror(s);
  exit(1);
}    


int main(int argc, char *argv[])
{
  uint8_t *Text;
  long n=0;
  int e;
  time_t start_wc = time(NULL);

  // check input data
  if(argc<2){
    printf("\nUsage:\n\t %s name\n\n", argv[0]);
    puts("Compute the BWT of file name and output it to name.Bwt");
    puts("The input file cannot contain the zero character");
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");
    
  // read main file
  char *name;
  FILE *fin = fopen(argv[1],"rb");
  if(fin==NULL) die("file open");
  // get file size
  if(fseek(fin,0,SEEK_END)!=0) die("fseek");
  n = ftell(fin);
  #if !M64
  if(n > 0x7FFFFFFE) {
    printf("Text size greater than  2^31-2!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
  #endif
  // ------ allocate and read text file
  Text = malloc(n+1);
  if(Text==NULL) die("malloc 1");
  rewind(fin);
  size_t s = fread(Text,1,n,fin);
  if(s!=n) {
    char *msg=NULL;
    int e= asprintf(&msg,"read parse error: %zu vs %ld\n", s,n); 
    (void) e; die(msg);
  }

  Text[n]=0; // sacak needs a 0 eos
  if(fclose(fin)!=0) die("Error closing text file"); 
  // -------- alloc and compute SA
  uint_t *SA = malloc((n+1)*sizeof(uint_t));
  if(SA==NULL) die("malloc 2");
  int depth = sacak(Text,SA, n+1);
  printf("SA computed with depth: %d\n", depth);
  // ---- output BWT
  e = asprintf(&name,"%s.Bwt",argv[1]);
  if(e<1) die("asprint error");  
  FILE *fbwt = fopen(name,"wb");
  free(name);
  if(fbwt==NULL) die("BWT open error");
  assert(SA[0]==n);
  if(fputc(Text[n-1],fbwt)==EOF) die("Error writing Bwt (1)");
  for(long i=1;i<=n;i++) {
    if(SA[i]>0) 
      e = fputc(Text[SA[i]-1],fbwt); // text char 
    else
      e = fputc(0,fbwt); // eof char
    if(e==EOF) die("Error writing Bwt (2)");
  }
  if(fclose(fbwt)!=0) die("Error closing BWT");;
  // deallocate
  free(SA);
  free(Text);
  printf("** Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));  
  return 0;
}

