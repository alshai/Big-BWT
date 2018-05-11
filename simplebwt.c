#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "gsacak.h"


void die(char *s)
{
  perror(s);
  exit(1);
}    


int main(int argc, char *argv[])
{
  uint8_t *Text;
  long n=0;
  int e;

  // intput data
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
  // ------ allocate and read text file
  Text = malloc(n+1);
  if(Text==NULL) die("malloc 1");
  rewind(fin);
  size_t s = fread(Text,1,n,fin);
  if(s!=n) {printf("%zu vs %zu\n", s,n); die("read parse error");}
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
  return 0;
}

