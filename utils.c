#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include "utils.h" 


// write error message and exit
void die(const char *s)
{
  perror(s);
  exit(1);
}    

int fd_open_aux_file(const char *base, const char *ext, int flags)
{
  char *name;
  int e = asprintf(&name,"%s.%s",base,ext);
  if(e<1) die("asprint error");
  int fd = open(name,flags,00666);
  if(fd<0) die(__func__);  
  free(name);
  return fd;
}

// open and return an auxiliary file
// open file named base.ext with mode mode
FILE *open_aux_file(const char *base, const char *ext, const char *mode)
{
  char *name;
  int e = asprintf(&name,"%s.%s",base,ext);
  if(e<1) die("asprint error");
  FILE *f = fopen(name,mode);
  if(f==NULL) die(name);  
  free(name);
  return f;
}

// open file named base.num.ext with mode mode
FILE *open_aux_file_num(const char *base, const char *ext, const int num, const char *mode)
{
  char *name;
  int e = asprintf(&name,"%s.%d.%s",base,num,ext);
  if(e<1) die("asprint error");
  FILE *f = fopen(name,mode);
  if(f==NULL) die(name);  
  free(name);
  return f;
}

// ------ functions handling file consisting of multiple segments 

// open multiple file 
mFile *mopen_aux_file(const char *base, const char *ext, int nsegs)
{
  mFile *f = malloc(sizeof(mFile));
  if(f==NULL) die("Out of memory");
  f->base = strdup(base);
  f->ext = strdup(ext);
  f->nsegs = nsegs;
  f->cur = 0;
  f->f = open_aux_file_num(f->base, f->ext, f->cur,"rb");
  return f;
}

// close a multiple file 
int mfclose(mFile *f) {
  FILE *aux = f->f;
  free(f->ext);
  free(f->base);
  free(f);
  return fclose(aux);
}

size_t mfread(void *ptr, size_t size, size_t nmemb, mFile *f)
{
  // try reading from current file
  while(1) {
    size_t s = fread(ptr,size,nmemb,f->f);
    if(s==nmemb) return s;
    if(s!=0) die(__func__); // a partial reading should never occur
    if(++f->cur==f->nsegs) return 0; // we are at the last segment
    if(fclose(f->f)!=0) die("Error closing old parse segment");
    f->f = open_aux_file_num(f->base, f->ext, f->cur,"rb");
  }
  assert(0);
  return 0;
}



// -------------------------------------

// extract an integer from a length n array containing IBYTES bytes per element
uint64_t get_myint(uint8_t *a, long n, long i)
{
  assert(i<n);
  assert(a!=NULL);
  long offset = (i+1)*IBYTES-1;
  uint64_t ai = 0;
  for(long j=0;j<IBYTES;j++) 
    ai = (ai << 8) | a[offset-j];
  return ai;
}

// write write an integer to file always using IBYTES
void write_myint(uint64_t u, FILE *f)
{
  assert(f!=NULL);
  size_t s = fwrite(&u,IBYTES,1,f);
  if(s!=1) die(__func__);
}

// extract an integer as in get_myint and write to file using write_myint
void get_and_write_myint(uint8_t *a, long n, long i, FILE *f)
{
  uint64_t u = get_myint(a,n,i);
  write_myint(u,f);
}
