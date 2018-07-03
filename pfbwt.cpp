/* **************************************************************************
 * pfbwt.cpp
 *  
 * Usage:
 *   pfbwt.x wsize file
 *
 * See newscan.cpp for usage 
 * 
 **************************************************************************** */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <semaphore.h>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
extern "C" {
#include "gsa/gsacak.h"
#include "utils.h"
}

using namespace std;
using namespace __gnu_cxx;

#define Buf_size 20
#define Min_bwt_range 100000

static size_t get_bwt_size(char *name);
long binsearch(uint_t x, uint_t a[], long n);
int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid);
void sa2da(uint_t sa[], int_t lcp[], uint8_t d[], long dsize, long dwords, int w);
void compute_dict_bwt_lcp(uint8_t *d, long dsize,long dwords, int w, uint_t **sap, int_t **lcpp);


// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
struct SeqId {
  uint32_t id;       // lex. id of the dictionary word to which the suffix belongs
  int remaining;     // remaining copies of the suffix to be considered  
  uint32_t *bwtpos;  // list of bwt positions of this dictionary word
  uint8_t char2write;// char to be written (is the one preceeding the suffix)

  // constructor
  SeqId(uint32_t i, int r, uint32_t *b, int8_t c) : id(i), remaining(r), bwtpos(b) {
    char2write = c;
  }

  // advance to the next bwt position, return false if there are no more positions 
  bool next() {
    remaining--;
    bwtpos += 1;
    return remaining>0;
  }
  bool operator<(const SeqId& a);
};

bool SeqId::operator<(const SeqId& a) {
    return *bwtpos > *(a.bwtpos);
}

// range in the suffix array of the dictionary 
typedef struct {
  long sa_start;   // starting position in the dictionary 
  long sa_end;     // end position in the dictionary
  uint8_t *bwt;    // starting position in the output bwt
  long count;      // chars to be written to the output bwt;
} sa_range;


// single struct representing the main data for bwt computation 
typedef struct {
  uint8_t *dict;     // dictionary
  uint_t *sa;        // suffix array for dict[]
  int_t *lcp;        // lcp array
  long dsize;        // size of dict[] sa[] lcp[] 
  uint8_t *last;     // array of last symbols 
  uint32_t *ilist;   // inverted list 
  uint32_t *istart;  // starting position inside inverted list 
  long dwords;       // number of words in the dictionary
  long full_words; 
  long easy_bwts; 
  long hard_bwts; 
  uint8_t *bwt;
  int w;
} main_data;

// working data to be passed to each consumer thread
typedef struct {
  uint8_t *dict;     // dictionary
  uint_t *sa;        // suffix array for d[]
  int_t *lcp;        // lcp array
  uint8_t *last;     // array of last symbols 
  uint32_t *ilist;   // inverted list 
  uint32_t *istart;  // starting position inside inverted list 
  uint32_t dwords;   // number of words in the dictionary 
  sa_range  buffer[Buf_size]; // shared producer/consumer buffer 
  int cindex;                 // consumer index in buffer
  pthread_mutex_t mutex_consumers; // mutex and semaphores 
  sem_t sem_free_slots;
  sem_t sem_data_items;   
} thread_data;


long add_sa_entries(main_data *d, long i)
{
  uint_t *eos = d->sa+1;
  long next = i+1;
  uint32_t seqid;
  // we are considering dict[sa[i]....]
  int_t suffixLen = getlen(d->sa[i],eos,d->dwords,&seqid);
  // ignore suffixes of lenght <= w
  if(suffixLen <= d->w) return next;
  // ----- simple case: the suffix is a full word 
  if(d->sa[i]==0 || d->dict[d->sa[i]-1]==EndOfWord) {
    d->full_words++;
    for(long j=d->istart[seqid];j<d->istart[seqid+1];j++)
      *(d->bwt++) = d->last[d->ilist[j]];
    return next; // done with this entry 
  }
  // ----- hard case: there can be a group of equal suffixes starting at i
  // save seqid and the corresponding char 
  vector<uint32_t> id2merge(1,seqid); 
  vector<uint8_t> char2write(1,d->dict[d->sa[i]-1]);
  while(next<d->dsize && d->lcp[next]>=suffixLen) {
    assert(d->lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
    assert(d->sa[next]>0 && d->dict[d->sa[next]-1]!=EndOfWord); // sa[next] cannot be a full word
    int_t nextsuffixLen = getlen(d->sa[next],eos,d->dwords,&seqid);
    assert(nextsuffixLen>=suffixLen);
    if(nextsuffixLen==suffixLen) {
      id2merge.push_back(seqid);           // sequence to consider
      char2write.push_back(d->dict[d->sa[next]-1]);  // corresponding char
      next++;
    }
    else break;
  }
  // numwords dictionary words contains the same suffix
  size_t numwords = id2merge.size(); 
  // possibly many words but same char?
  bool samechar=true;
  for(size_t i=1;(i<numwords)&&samechar;i++)
    samechar = (char2write[i-1]==char2write[i]); 
  if(samechar) {
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i];
      for(long j=d->istart[s];j<d->istart[s+1];j++)
        *(d->bwt++) = char2write[0];
      d->easy_bwts +=  d->istart[s+1]- d->istart[s]; 
    }
    return next;
  }
  // many words, many chars...     
  {
    // create heap
    vector<SeqId> heap;
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s,d->istart[s+1]-d->istart[s], d->ilist+d->istart[s], char2write[i]));
    }
    std::make_heap(heap.begin(),heap.end());
    while(heap.size()>0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      *(d->bwt++) = s.char2write;
      d->hard_bwts += 1;
      // remove top 
      pop_heap(heap.begin(),heap.end());
      heap.pop_back();
      // if remaining positions, reinsert to heap
      if(s.next()) {
        heap.push_back(s);
        push_heap(heap.begin(),heap.end());
      }
    }
  }
  return next;
}

#if 0
// code for the consumer threads 
void tbody(void *v)
{
  thread_data *d = v;

  // main loop
  while(true) {
    int e = sem_wait(d->sem_data_items);
    if(e!=0) die("tbody:sem wait");
    e = pthread_mutex_lock(d->mutex_consumers);
    if(e!=0) die("tbody:mutex lock");
    sa_range r  = d->buffer[d->cindex];
    d->cindex = (d->cindex + 1) % Buf_size;
    e = pthread_mutex_unlock(a->mutex_consumers);
    if(e!=0) die("tbody:mutex unlock");
    e = sem_post(a->sem_free_slots);
    if(e!=0) die("tbody:sem post");
    if(r.sa_start == r.sa_end) break
    // process range r
    uint32_t seqid; 
    long i, next;
    uint8_t *bwt_start = *bwt = r.bwt;
    uint_t *eos = d->sa + 1;
    for(long i=r.sa_start; i< r.sa_end; i=next ) {
      // we are considering d[sa[i]....]
      next = i+1;  // prepare for next iteration  
      // compute length of this suffix and sequence it belongs
      int_t suffixLen = getlen(d->sa[i],eos,d->dwords,&seqid);
      // ignore suffixes of lenght <= w
      if(suffixLen<=w) continue;
      // ----- simple case: the suffix is a full word 
      if(d->sa[i]==0 || d->dict[s->sa[i]-1]==EndOfWord) {
        for(long j=d->istart[seqid];j<d->istart[seqid+1];j++)
          *bwt++ = d->last[d->ilist[j]];
        continue; // proceed with next i 
      }
      // ----- hard case: there can be a group of equal suffixes starting at i
      // save seqid and the corresponding char 
      vector<uint32_t> id2merge(1,seqid); 
      vector<uint8_t> char2write(1,d[sa[i]-1]);
      while(next<dsize && lcp[next]>=suffixLen) {
        assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
        assert(sa[next]>0 && d[sa[next]-1]!=EndOfWord); // sa[next] cannot be a full word
        int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
        assert(nextsuffixLen>=suffixLen);
        if(nextsuffixLen==suffixLen) {
          id2merge.push_back(seqid);           // sequence to consider
          char2write.push_back(d[sa[next]-1]);  // corresponding char
          next++;
        }
        else break;
      }
      size_t numwords = id2merge.size(); 
      // numwords dictionary words contains the same suffix
      // case of a single word
      if(numwords==1) {
        uint32_t s = id2merge[0];
        for(long j=istart[s];j<istart[s+1];j++)
          if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
        easy_bwts +=  istart[s+1]- istart[s]; 
        continue;   
      }
      // many words, same char?
      bool samechar=true;
      for(size_t i=1;(i<numwords)&&samechar;i++)
        samechar = (char2write[i-1]==char2write[i]); 
      if(samechar) {
        for(size_t i=0; i<id2merge.size(); i++) {
          uint32_t s = id2merge[i];
          for(long j=istart[s];j<istart[s+1];j++)
            if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 2");
          easy_bwts +=  istart[s+1]- istart[s]; 
        }
        continue;
      }
      // many words, many chars...     
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
          hard_bwts += 1;
          // remove top 
          pop_heap(heap.begin(),heap.end());
          heap.pop_back();
          // if remaining positions, reinsert to heap
          if(s.next()) {
            heap.push_back(s);
            push_heap(heap.begin(),heap.end());
          }
        }
      }
    }
    assert(i==d.sa_end);                 // make sure we have the exact range 
    assert(bwt_start + r.count == bwt);  // make sure we wrote the expect number of chars 
  }
  pthread_exit(NULL);  
}  
#endif


void bwt_multi_thread(uint8_t *d, long dsize, // dictionary and its size  
         uint32_t *ilist, uint8_t *last, long psize, // ilist, last and their size  
         uint32_t *istart, long dwords, // starting point in ilist for each word and # words
         int w, char *name, int numt)
{
  (void) psize; // used only in assertions
  // compute sa and bwt of d[] and do some checking on them 
  uint_t *sa; int_t *lcp; 
  compute_dict_bwt_lcp(d,dsize,dwords,w,&sa,&lcp);
  // set d[0] ==0 as this is the EOF char in the final BWT
  assert(d[0]==Dollar);
  d[0]=0;
  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
  uint_t *eos = sa+1;
  for(int i=0;i<dwords-1;i++)
    assert(eos[i]<eos[i+1]);
  
  // save everything in the thread_data structure
  main_data data;
  data.dict = d; data.sa = sa; data.lcp = lcp; data.dsize = dsize; 
  data.last = last; data.ilist = ilist; data.istart = istart;
  data.dwords = dwords; 
  data.full_words = data.easy_bwts = data.hard_bwts = 0;
  data.w = w;
  //data.cindex = 0;
  //data.mutex_consumers = PTHREAD_MUTEX_INITIALIZER;
  //sem_init(&data.sem_data_items,0,0);
  //sem_init(&data.sem_free_slots,0,Buf_size);

  // get final bwt size from the size of the input file    
  size_t bwt_size= get_bwt_size(name);
  // open output file and map it to the bwt array 
  FILE *fbwt = open_aux_file(name,"bwt","wb+");
  // make the BWT file of the correct size (otherwise mmap fails)
  if(ftruncate(fileno(fbwt),bwt_size)<0) die("truncate failed");
  data.bwt = (uint8_t *) mmap(NULL,bwt_size,PROT_READ|PROT_WRITE,MAP_SHARED,fileno(fbwt), 0);
  if(data.bwt==MAP_FAILED) die("mmap failed");
  // main loop: consider each entry in the SA of dict
  time_t start = time(NULL);
  if(numt==0) { // no helper threads
    for(long i=dwords+w+1; i< dsize;)
      i = add_sa_entries(&data,i);
  }
  else {
    // create threads
    // partition BWT chars 
    uint32_t seqid;
    long full_words = 0;
    long next, written = dwords+w+1, entries = 0;
    for(long i=dwords+w+1; i< dsize; i=next ) {
      // ---- check if a batch is ready
      if(entries >= 100000) {
        printf("%ld -- %ld (%ld)\n",written,written+entries, i);
        //for(j=written; j< written+entries; )
        //  j =  add_sa_entries(&data,j);
        //assert(j==written+entries);
        written += entries; entries=0;
      }
      // we are considering d[sa[i]....]
      next = i+1;  // prepare for next iteration  
      int_t suffixLen = getlen(sa[i],eos,dwords,&seqid);
      // ignore suffixes of lenght <= w
      if(suffixLen<=w) continue;
      // simple case: the suffix is a full word 
      if(sa[i]==0 || d[sa[i]-1]==EndOfWord) {
        full_words++;
        entries += istart[seqid+1] - istart[seqid];
        continue; // proceed with next i 
      }
      // hard case: there can be a group of equal suffixes starting at i
      entries += istart[seqid+1] - istart[seqid];
      while(next<dsize && lcp[next]>=suffixLen) {
        int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
        if(nextsuffixLen==suffixLen) {
          next++;
          entries += istart[seqid+1] - istart[seqid];
        }
        else break;
      }
    }
    assert(full_words==dwords);
    printf("%ld == %ld\n",written,written+entries);
    data.full_words = full_words;
  }
  munmap(data.bwt, 0);
  assert(data.full_words==dwords);
  cout << "Full words: " << data.full_words << endl;
  cout << "Easy bwt chars: " << data.easy_bwts << endl;
  cout << "Hard bwt chars: " << data.hard_bwts << endl;
  cout << "Generating the final BWT took " << difftime(time(NULL),start) << " wall clock seconds\n";    
  fclose(fbwt);
  delete[] lcp;
  delete[] sa;
}

/* *******************************************************************
 * Computation of the final BWT
 * 
 * istart[] and islist[] are used together. For each dictionary word i 
 * (considered in lexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i 
 * ******************************************************************* */
void bwt(uint8_t *d, long dsize, // dictionary and its size  
         uint32_t *ilist, uint8_t *last, long psize, // ilist, last and their size 
         uint32_t *istart, long dwords, // starting point in ilist for each word and # words
         int w, char *name)             // window size and base name for output file
{  
  (void) psize; // used only in assertions
  
  // compute sa and bwt of d and do some checking on them 
  uint_t *sa; int_t *lcp; 
  compute_dict_bwt_lcp(d,dsize,dwords,w,&sa,&lcp);
  // set d[0] ==0 as this is the EOF char in the final BWT
  assert(d[0]==Dollar);
  d[0]=0;

  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
  uint_t *eos = sa+1;
  for(int i=0;i<dwords-1;i++)
    assert(eos[i]<eos[i+1]);

  // open output file 
  FILE *fbwt = open_aux_file(name,"bwt","wb");
    
  // main loop: consider each entry in the SA of dict
  time_t start = time(NULL);
  long full_words = 0; 
  long easy_bwts = 0;
  long hard_bwts = 0;
  long next;
  uint32_t seqid;
  for(long i=dwords+w+1; i< dsize; i=next ) {
    // we are considering d[sa[i]....]
    next = i+1;  // prepare for next iteration  
    // compute length of this suffix and sequence it belongs
    int_t suffixLen = getlen(sa[i],eos,dwords,&seqid);
    // ignore suffixes of lenght <= w
    if(suffixLen<=w) continue;
    // ----- simple case: the suffix is a full word 
    if(sa[i]==0 || d[sa[i]-1]==EndOfWord) {
      full_words++;
      for(long j=istart[seqid];j<istart[seqid+1];j++)
        if(fputc(last[ilist[j]],fbwt)==EOF) die("BWT write error");
      continue; // proceed with next i 
    }
    // ----- hard case: there can be a group of equal suffixes starting at i
    // save seqid and the corresponding char 
    vector<uint32_t> id2merge(1,seqid); 
    vector<uint8_t> char2write(1,d[sa[i]-1]);
    while(next<dsize && lcp[next]>=suffixLen) {
      assert(lcp[next]==suffixLen);  // the lcp cannot be greater than suffixLen
      assert(sa[next]>0 && d[sa[next]-1]!=EndOfWord); // sa[next] cannot be a full word
      int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
      assert(nextsuffixLen>=suffixLen);
      if(nextsuffixLen==suffixLen) {
        id2merge.push_back(seqid);           // sequence to consider
        char2write.push_back(d[sa[next]-1]);  // corresponding char
        next++;
      }
      else break;
    }
    size_t numwords = id2merge.size(); 
    // numwords dictionary words contains the same suffix
    // case of a single word
    if(numwords==1) {
      uint32_t s = id2merge[0];
      for(long j=istart[s];j<istart[s+1];j++)
        if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      easy_bwts +=  istart[s+1]- istart[s]; 
      continue;   
    }
    // many words, same char?
    bool samechar=true;
    for(size_t i=1;(i<numwords)&&samechar;i++)
      samechar = (char2write[i-1]==char2write[i]); 
    if(samechar) {
      for(size_t i=0; i<id2merge.size(); i++) {
        uint32_t s = id2merge[i];
        for(long j=istart[s];j<istart[s+1];j++)
          if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 2");
        easy_bwts +=  istart[s+1]- istart[s]; 
      }
      continue;
    }
    // many words, many chars...     
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
        hard_bwts += 1;
        // remove top 
        pop_heap(heap.begin(),heap.end());
        heap.pop_back();
        // if remaining positions, reinsert to heap
        if(s.next()) {
          heap.push_back(s);
          push_heap(heap.begin(),heap.end());
        }
      }
    }
  }  
  assert(full_words==dwords);
  cout << "Full words: " << full_words << endl;
  cout << "Easy bwt chars: " << easy_bwts << endl;
  cout << "Hard bwt chars: " << hard_bwts << endl;
  cout << "Generating the final BWT took " << difftime(time(NULL),start) << " wall clock seconds\n";    
  fclose(fbwt);
  delete[] lcp;
  delete[] sa;
}  


void bwt_new(uint8_t *d, long dsize, // dictionary and its size  
         uint32_t *ilist, uint8_t *last, long psize, // ilist, last and their size 
         uint32_t *istart, long dwords, // starting point in ilist for each word and # words
         int w, char *name)             // window size and base name for output file
{  
  (void) psize; // used only in assertions
  // open output file 
  FILE *fbwt = open_aux_file(name,"bwt","wb");
  time_t  start = time(NULL);  
  
  // compute sa and bwt of d and do some checking on them 
  uint_t *sa; int_t *lcp; 
  compute_dict_bwt_lcp(d,dsize,dwords,w,&sa,&lcp);
  // set d[0] ==0 as this is the EOF char in the final BWT
  assert(d[0]==Dollar);
  d[0]=0;

  // convert sa,lcp->da,suflen + bit
  sa2da(sa,lcp,d,dsize,dwords,w);
  uint_t *da = sa + (dwords+w+1);
  uint_t *eos = sa+1;
  long dasize= dsize - (dwords+w+1);
  int_t *suflen = lcp + (dwords+w+1);
  int_t *wlen = lcp+1;
  lcp = NULL; sa = NULL; // make sure these are not used

  // main loop: consider each entry in the DA[] of dict
  long full_words = 0; 
  long easy_bwts = 0;
  long hard_bwts = 0;
  long next;
  for(long i=0; i< dasize; i=next ) {
    // we are considering d[sa[i]....] belonging to da[i]
    next = i+1;  // prepare for next iteration  
    // discard if it is a small suffix 
    if(suflen[i]<=w) continue;
    uint32_t seqid = da[i]&0x7FFFFFFF;
    assert(seqid<dwords);

    // ----- simple case: the suffix is a full word 
    if(suflen[i]==wlen[seqid]) {
      full_words++;
      for(long j=istart[seqid];j<istart[seqid+1];j++)
        if(fputc(last[ilist[j]],fbwt)==EOF) die("BWT write error");
      continue; // proceed with next i 
    }
    // ----- hard case: there can be a group of equal suffixes starting at i
    // save seqid and the corresponding char 
    vector<uint32_t> id2merge(1,seqid); 
    vector<uint8_t> char2write(1,d[eos[seqid]-suflen[i]-1]);
    while(next<dasize && suflen[next]==suflen[i]) {
      seqid = da[next]&0x7FFFFFFF;
      if(da[next]&0x80000000u) {
        assert(suflen[next]!=wlen[seqid]);   // the lcp cannot be greater than suffixLen
        id2merge.push_back(seqid);           // sequence to consider
        char2write.push_back(d[eos[seqid]-suflen[next]-1]);  // corresponding char
        next++;
      }
      else break;
    }
    size_t numwords = id2merge.size(); 
    // numwords dictionary words contains the same suffix
    // case of a single word
    if(numwords==1) {
      uint32_t s = id2merge[0];
      for(long j=istart[s];j<istart[s+1];j++)
        if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      easy_bwts +=  istart[s+1]- istart[s]; 
      continue;   
    }
    // many words, same char?
    bool samechar=true;
    for(size_t i=1;(i<numwords)&&samechar;i++)
      samechar = (char2write[i-1]==char2write[i]); 
    if(samechar) {
      for(size_t i=0; i<id2merge.size(); i++) {
        uint32_t s = id2merge[i];
        for(long j=istart[s];j<istart[s+1];j++)
          if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 2");
        easy_bwts +=  istart[s+1]- istart[s]; 
      }
      continue;
    }
    // many words, many chars...     
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
        hard_bwts += 1;
        // remove top 
        pop_heap(heap.begin(),heap.end());
        heap.pop_back();
        // if remaining positions, reinsert to heap
        if(s.next()) {
          heap.push_back(s);
          push_heap(heap.begin(),heap.end());
        }
      }
    }
  }  
  if(full_words!=dwords) 
    cerr << "Dwords: " << dwords << endl;
  cout << "Full words: " << full_words << endl;
  cout << "Easy bwt chars: " << easy_bwts << endl;
  cout << "Hard bwt chars: " << hard_bwts << endl;
  cout << "Generating the final BWT took " << difftime(time(NULL),start) << " wall clock seconds\n";    
  fclose(fbwt);
  delete[] lcp;
  delete[] sa;
} 


// compute the number of words in a dictionary
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
  if(argc!=3 && argc!=4) {
    cerr << "Usage:\n\t";
    cerr << argv[0] << " wsize file [threads]\n" << endl;
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  // translate command line parameters
  int w = atoi(argv[1]);             // sliding window size   
  int num_threads=0;                 // number of helper threads
  if(argc==4) num_threads = atoi(argv[3]);
  if(num_threads<0) {
    cerr << "Number of helper threads cannot be negative!\n";
    exit(1);
  }

  // read dictionary file 
  FILE *g = open_aux_file(argv[2],"dict","rb");
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
  g = open_aux_file(argv[2],"occ","rb");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell 2");
  if(e%4!=0) die("invalid occ file");
  int dwords = e/4;
  cout  << "Dictionary words: " << dwords << endl;
  uint32_t *occ = new uint32_t[dwords+1];  // dwords+1 since istart overwrites occ
  rewind(g);
  e = fread(occ,4,dwords,g);
  if(e!=dwords) die("fread 2");
  fclose(g);
  assert(dwords==get_num_words(d,dsize));

  // read ilist file 
  g = open_aux_file(argv[2],"ilist","rb");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell 3");
  if(e%4!=0) die("invalid ilist file");
  long psize = e/4;
  cout  << "Parsing size: " << psize << endl;
  if(psize>0xFFFFFFFEL) die("More than 2^32 -2 words in the parsing");
  uint32_t *ilist = new uint32_t[psize];  
  rewind(g);
  e = fread(ilist,4,psize,g);
  if(e!=psize) die("fread 3");
  fclose(g);
  assert(ilist[0]==1); // EOF is in PBWT[1] 

  // read bwlast file 
  g = open_aux_file(argv[2],"bwlast","rb");  
  cout  << "bwlast file size: " << psize << endl;
  uint8_t *bwlast = new uint8_t[psize];  
  e = fread(bwlast,1,psize,g);
  if(e!=psize) die("fread 4");
  fclose(g);

  // convert occ entries into starting positions inside ilist
  // ilist also contains the position of EOF but we don't care about it since it is not in dict 
  uint32_t last=1; // starting position in ilist of the smallest dictionary word  
  for(long i=0;i<dwords;i++) {
    uint32_t tmp = occ[i];
    occ[i] = last;
    last += tmp;
  }
  assert(last==psize);
  occ[dwords]=psize;
  // extra check: the smallest dictionary word is d0 =$.... that occurs once
  assert(occ[1]==occ[0]+1);
  
  // compute and write the final bwt 
  bwt_new(d,dsize,ilist,bwlast,psize,occ,dwords,w,argv[2]);
  //bwt_multi_thread(d,dsize,ilist,bwlast,psize,occ,dwords,w,argv[2],num_threads);
  // old version not using threads and working correctly
  //bwt(d,dsize,ilist,bwlast,psize,occ,dwords,w,argv[2]);
  
  delete[] bwlast;
  delete[] ilist;
  delete[] occ;
  delete[] d;  
  return 0;
}

// --------------------- aux functions ----------------------------------


// compute size of the bwt adding 1 to the input size
static size_t get_bwt_size(char *name)
{
  FILE *f = fopen(name,"rb");
  if(f==NULL) die("Input file open");
  int e = fseek(f,0,SEEK_END);
  if(e<0) die("Input file seek");
  long s = ftell(f);
  if(s<0) die("Input file tell");
  if(fclose(f)!=0) die("Input file close");
  cerr << "input file size: " << s << endl;
  return 1 + s;
}


// binary search for x in an array a[0..n-1] that doesn't contain x
// return the lowest position that is larger than x
long binsearch(uint_t x, uint_t a[], long n)
{
  long lo=0; long hi = n-1;
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



// transform sa[],lcp[] -> da[], suflen[] + 
// extra bit telling whether suflen[i]==lcp[i]
void sa2da(uint_t sa[], int_t lcp[], uint8_t d[], long dsize, long dwords, int w)
{
  long words=0;
  if(dwords>0x7FFFFFFF) {
    cerr << "Too many words in the dictionary. Current limit: 2^31-1\n";
    exit(1);
  }
  cout << "Converting SA and LCP Array to DA and SufLen\n";
  
  time_t  start = time(NULL);  
  // create eos[] array with ending position in d[] of each word
  uint_t *eos = sa + 1;
  int_t *wlen = lcp + 1;
  // save length of word i in lcp[i+1]==suflen[i] (sa[i+1]=eos[i] is the position of its eos)
  wlen[0] = eos[0];
  for(long i=1;i<dwords;i++) {
    wlen[i] = eos[i]-eos[i-1] -1;
    assert(wlen[i]>0);
    assert(d[eos[i-1]]==EndOfWord);
    assert(d[eos[i-1]+wlen[i]+1]==EndOfWord);
  }
  // convert sa,lcp -> da,suflen  
  uint32_t seqid; 
  for(long i=dwords+w+1; i<dsize; i++) {     // we are considering d[sa[i]....]     
    int_t suffixLen = getlen(sa[i],eos,dwords,&seqid);
    assert(seqid<=0x7FFFFFFF);     // seqid uses at most 31 bits
    assert(suffixLen>=lcp[i]);     // suffix length cannot be shorter than lcp
    assert(suffixLen<=wlen[seqid]);// suffix length cannot be larger than word length
    if(suffixLen==wlen[seqid]) {   // test if full word
      words++;                     
      assert(lcp[i]<suffixLen);    // full words are not prefix of other suffixes
    }
    if(lcp[i]==suffixLen) {         // save seqid + possibily extra bit 
      sa[i] = seqid | (1u << 31);  // mark last bit if lcp==suffix_len;
    }
    else 
      sa[i] = seqid;               // save only seqid = da[i]
    lcp[i] = suffixLen;            // save suffix length overwriting lcp
  }
  cout << "Conversion took " << difftime(time(NULL),start) << " wall clock seconds\n";  
  cout << "DA has size: " << dsize-dwords-w-1;
  cout << ". Dictionary words found: " << words << endl; 
}


// return the length of the suffix starting in position p.
// also write to seqid the id of the sequence containing that suffix 
// n is the # of distinct words in the dictionary, hence the length of eos[]
int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid)
{
  assert(p<eos[n-1]);
  *seqid = binsearch(p,eos,n);
  assert(eos[*seqid]> p); // distance between position p and the next $
  return eos[*seqid] - p;
}

// compute the SA and LCP array for the set of (unique) dictionary words
// using gSACA-K. Also do some checking based on the number and order of the special symbols
void compute_dict_bwt_lcp(uint8_t *d, long dsize,long dwords, int w, 
                          uint_t **sap, int_t **lcpp) // output parameters
{
  uint_t *sa = new uint_t[dsize];
  int_t *lcp = new int_t[dsize];
  (void) dwords; (void) w;

  cout  << "Each SA entry: " << sizeof(*sa) << " bytes\n";
  cout  << "Each LCP entry: " << sizeof(*lcp) << " bytes\n";

  cout << "Computing SA and LCP of dictionary" << endl; 
  time_t  start = time(NULL);
  gsacak(d,sa,lcp,NULL,dsize);
  cout << "Computing SA/LCP took " << difftime(time(NULL),start) << " wall clock seconds\n";  
  // ------ do some checking on the sa
  assert(d[dsize-1]==EndOfDict);
  assert(sa[0]==(unsigned long)dsize-1);// sa[0] is the EndOfDict symbol 
  for(long i=0;i<dwords;i++) 
    assert(d[sa[i+1]]==EndOfWord); // there are dwords EndOfWord symbols 
  // EndOfWord symbols are in position order, so the last is d[dsize-2]    
  assert(sa[dwords]==(unsigned long)dsize-2);  
  // there are wsize+1 $ symbols: 
  // one at the beginning of the first word, wsize at the end of the last word
  for(long i=0;i<=w;i++)
    assert(d[sa[i+dwords+1]]==Dollar);         
  // in sa[dwords+w+1] we have the first word in the parsing since that $ is the lex. larger  
  assert(d[0]==Dollar);
  assert(sa[dwords+w+1]==0);
  assert(d[dwords+w+2]>Dollar);  // end of Dollar chars 
  assert(lcp[dwords+w+2]==0); 
  // copy sa and lcp address
  *sap = sa;  *lcpp = lcp;  
}
