/* ******************************************************************************
 * pfbwt.cpp
 *  * 
 * Usage:
 *   pfbwt.x wsize file
 *
 * See newscan.cpp for usage 
 * 
 */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
extern "C" {
#include "gsacak.h"
}

using namespace std;
using namespace __gnu_cxx;

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char 
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter 


struct SeqId {
  uint32_t id;       // id of the sequence
  int remaining;     // remaining copies to be considered  
  uint32_t *bwtpos;  // list of bwt positions for this id
  uint8_t char2write;// char to be writte 

  // constructor
  SeqId(uint32_t i, int r, uint32_t *b, int8_t c) : id(i), remaining(r), bwtpos(b) {
    char2write = c;
  }

  // got to the next bwt position, return false if there are no more positions 
  bool next() {
    remaining--;
    bwtpos += 1;
    return remaining >0;
  }
    
  bool operator<(const SeqId& a);
  
};

bool SeqId::operator<(const SeqId& a) {
    return *bwtpos > *(a.bwtpos);
}



void die(string s)
{
  perror(s.c_str());
  exit(1);
}


// binary search for x in an array a[0..n-1] that does not contain x
// return the lowest position that is larger than x
int binsearch(uint32_t x, uint32_t a[], int n)
{
  int lo=0; int hi = n-1;
  while(hi>lo) {
    assert( ((lo==0) || x>a[lo-1]) && x< a[hi]);
    int mid = (lo+hi)/2;
    assert(x!=a[mid]);  // x is not in a[]
    if(x<a[mid]) hi = mid;
    else lo = mid+1;
  }
  assert(((hi==0) || x>a[hi-1]) && x< a[hi]);
  return hi; 
}

// return the length of the suffix startin in position p.
// also write to seqid the id of the sequecne containing that suffix 
int getlen(uint32_t p, uint32_t eos[], int n, uint32_t *seqid)
{
  assert(p<eos[n-1]);
  *seqid = binsearch(p,eos,n);
  return eos[*seqid] - p;
}

void bwt(uint8_t d[], long dsize, uint32_t ilist[], uint8_t last[], long psize, uint32_t istart[], int dwords, int w, string name)
{
  
  // compute sa and bwt of d
  uint32_t *sa = new uint32_t[dsize];
  int32_t *lcp = new int32_t[dsize];
  (void) psize;

  cout << "Computing SA and LCP of dictionary" << endl; 
  gsacak(d,sa,lcp,NULL,dsize);
  // do some checking on the sa
  assert(d[dsize-1]==EndOfDict);
  assert(sa[0]==dsize-1);
  for(int i=0;i<dwords;i++) 
    assert(d[sa[i+1]]==EndOfWord); // there are dwords eos symbols  
  assert(sa[dwords]==dsize-2);  
  for(int i=0;i<=w;i++)
    assert(d[sa[i+dwords+1]]==Dollar); // there are wsize+1 $ symbols        
  // in sa[dwords+w+1] we have the first word in the parsing   
  assert(d[0]==Dollar);
  assert(sa[dwords+w+1]==0);
  assert(lcp[dwords+w+2]==0); // end of Dollar chars 
  // set d[0] ==0 as this it te EOF chars in the file BWT
  d[0]=0;
  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
  uint32_t *eos = sa+1;

  // open output file 
  FILE *fbwt = fopen(name.c_str(),"wb");
  if(fbwt==NULL) die("Open bwt file");
  
  // main loop
  int full_words = 0; 
  long next;
  uint32_t seqid;
  for(long i=dwords+w+1; i< dsize; i=next ) {
    // we are considering t[sa[i]....]
    next = i+1;  // prepare for next iteration  
    int suffixLen = getlen(sa[i],eos,dwords,&seqid);
    // ignore suffixes of lenght <= w
    if(suffixLen<=w) continue;
    // simple case: the suffix is a full word 
    if(sa[i]==0 || d[sa[i]-1]==EndOfWord) {
      full_words++;
      for(uint32_t j=istart[seqid];j<istart[seqid+1];j++)
        if(fputc(last[ilist[j]],fbwt)==EOF) die("BWT write error");
      continue;
    }
    // hard case: there can be a group of equal suffixes starting at i
    // push seqid
    vector<uint32_t> id2merge(1,seqid); 
    vector<uint8_t> char2write(1,d[sa[i]-1]);
    while(next<dsize && lcp[next]>=suffixLen) {
      assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
      assert(sa[next]>0 && d[sa[next]-1]!=EndOfWord); // sa[end] cannot be a full word
      int nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
      assert(nextsuffixLen>=suffixLen);
      if(nextsuffixLen==suffixLen) {
        id2merge.push_back(seqid);           // sequence to consider
        char2write.push_back(d[sa[next]-1]);  // corresponding char
        next++;
      }
      else break;
    }
    size_t numwords = id2merge.size(); 
    // case of a single word
    if(numwords==1) {
      uint32_t s = id2merge[0];
      for(uint32_t j=istart[s];j<istart[s+1];j++)
        if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      continue;   
    }
    // many words, same char
    bool samechar=true;
    for(size_t i=1;(i<numwords)&&samechar;i++)
      samechar = (char2write[i-1]==char2write[i]); 
    if(samechar) {
      for(size_t i=0; i<id2merge.size(); i++) {
        uint32_t s = id2merge[i];
        for(uint32_t j=istart[s];j<istart[s+1];j++)
          if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 2");
      }
      continue;
    }
    // many words, many chars...     
    // process stack
    {
      // create heap
      vector<SeqId> heap;
      for(size_t i=0; i<numwords; i++) {
        uint32_t s = id2merge[i];
        heap.push_back(SeqId(s,istart[s+1]-istart[s], ilist+istart[s], char2write[i]));
      }
      std::make_heap(heap.begin(),heap.end());
      while(heap.size()>0) {
        // output char for the top of the heap
        SeqId s = heap.front();
        if(fputc(s.char2write,fbwt)==EOF) die("BWT write error 3");
        // remove top 
        pop_heap(heap.begin(),heap.end());
        heap.pop_back();
        // if remaining positions, reinsert to heap
        if(s.next()) {
          heap.push_back(s);
          push_heap(heap.begin(),heap.end());
        }
      }
      /*
      for(size_t i=0; i<numwords; i++) {
        cout << "Id: " << heap.front().id << "  Bwtpos: " << *(heap.front().bwtpos) << endl;
        pop_heap(heap.begin(),heap.end()); 
        heap.pop_back();
      }
      cout << "----------\n"; */   
    }
  }  
  assert(full_words==dwords);
  fclose(fbwt);
  delete[] lcp;
  delete[] sa;
}  

// compute the numbber of words in a dictionary
long get_num_words(uint8_t *d, long n)
{
  long i,num=0;
  for(i=0;i<n;i++)
    if(d[i]==EndOfWord) num++;
  assert(d[n-1]==EndOfDict);
  return num;
}

int main(int argc, char** argv)
{
  // check command line
  if(argc!=3) {
    cerr << "Usage:\n\t";
    cerr << argv[0] << " wsize file\n" << endl;
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  // translate command line parameters
  int w = atoi(argv[1]);             // sliding window size 
  string fname(argv[2]);              // basename for input files   

  // read dictionary file 
  string name = fname + ".dict";
  FILE *g = fopen(name.c_str(),"rb");
  if(g==NULL) die("Open dict file");
  fseek(g,0,SEEK_END);
  long dsize = ftell(g);
  if(dsize<0) die("ftell");
  if(dsize<=1+w) die("invalid dictionary file");
  cout  << "Dictionary file size: " << dsize << endl;
  uint8_t *d = new uint8_t[dsize];  
  rewind(g);
  long e = fread(d,1,dsize,g);
  if(e!=dsize) die("fread");
  fclose(g);
  
  // read occ file
  name = fname + ".occ";
  g = fopen(name.c_str(),"rb");
  if(g==NULL) die("Open occ file");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell 2");
  if(e%4!=0) die("invalid occ file");
  int dwords = e/4;
  cout  << "Dictionary words: " << dwords << endl;
  uint32_t *occ = new uint32_t[dwords+1];  
  rewind(g);
  e = fread(occ,4,dwords,g);
  if(e!=dwords) die("fread 2");
  fclose(g);
  assert(dwords==get_num_words(d,dsize));

  // read ilist file 
  name = fname + ".ilist";
  g = fopen(name.c_str(),"rb");
  if(g==NULL) die("Open ilist file");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell 3");
  if(e%4!=0) die("invalid ilist file");
  long psize = e/4;
  cout  << "Parsing size: " << psize << endl;
  uint32_t *ilist = new uint32_t[psize];  
  rewind(g);
  e = fread(ilist,4,psize,g);
  if(e!=psize) die("fread 3");
  fclose(g);
  assert(ilist[0]==1); // EOF is in PBWT[1]

  // read bwlast file 
  name = fname + ".bwlast";
  g = fopen(name.c_str(),"rb");
  if(g==NULL) die("Open bwlast file");
  cout  << "bwlast file size: " << psize << endl;
  uint8_t *bwlast = new uint8_t[psize];  
  e = fread(bwlast,1,psize,g);
  if(e!=psize) die("fread 4");
  fclose(g);

  // convert occ entries into starting positions inside ilist
  uint32_t last=1;
  for(long i=0;i<dwords;i++) {
    uint32_t tmp = occ[i];
    occ[i] = last;
    last += tmp;
  }
  assert(last==psize);
  occ[dwords]=psize;
      
  bwt(d,dsize,ilist,bwlast,psize,occ,dwords,w,fname+".bwt");
  
  delete[] bwlast;
  delete[] ilist;
  delete[] occ;
  delete[] d;  
  return 0;
}


