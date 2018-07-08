#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


// write error message and exit
void die(const char *s)
{
  perror(s);
  exit(1);
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

int fd_open_aux_file(const char *base, const char *ext, int mode)
{
  char *name;
  int e = asprintf(&name,"%s.%s",base,ext);
  if(e<1) die("asprint error");
  int fd = open(name,mode);
  if(fd<0) die(__func__);  
  free(name);
  return fd;
}
