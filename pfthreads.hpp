/* **************************************************************************
 * pfthreads.hpp
 * 
 **************************************************************************** */
extern "C" {
#include "xerrors.h"
}



#define Buf_size 20
#define Min_bwt_range 1000000
#define Sa_block      1000000


// ----- parallel conversion of sa/lcp ->da/sufLen
typedef struct {
  long start;
  long end;
} sarange;

typedef struct{
    sarange *buffer;                // buffer prod/consumer
    int *cindex;                    // consumer index in buffer
    pthread_mutex_t *cindex_m;      // mutex for c_index
    sem_t *free_slots, *data_items; // prod/consumer semaphores
    uint_t *sa, *eos;
    int_t *lcp, *wlen;
    long dwords;
    long full_words;                // full words found, shared: use mutex 
} sa2da_data;


// initialize/destroy semaphores and mutex for producer/consumer 
static void pc_init(sem_t *free_slots, sem_t *data_items, pthread_mutex_t *m)
{
  xpthread_mutex_init(m,NULL,__LINE__,__FILE__);
  xsem_init(free_slots,0,Buf_size,__LINE__,__FILE__);
  xsem_init(data_items,0,0,__LINE__,__FILE__);
}  
  
static void pc_destroy(sem_t *free_slots, sem_t *data_items, pthread_mutex_t *m)
{
  xpthread_mutex_destroy(m,__LINE__,__FILE__);
  xsem_destroy(free_slots,__LINE__,__FILE__);
  xsem_destroy(data_items,__LINE__,__FILE__);
}  
  
  
void *sa2da_body(void *v)
{
  sa2da_data *d = (sa2da_data *) v;
  uint32_t seqid; 
  long words=0;
  while(true) {
    // --- get starting position from buffer 
    int e = sem_wait(d->data_items);
    if(e) die("wait in sa2da_body");
    e = pthread_mutex_lock(d->cindex_m);
    if(e) die("lock in sa2da_body");
    sarange r = d->buffer[*(d->cindex)]; 
    *(d->cindex) = (*(d->cindex) + 1) % Buf_size;
    e = pthread_mutex_unlock(d->cindex_m);
    if(e) die("unlock in sa2da_body");
    e = sem_post(d->free_slots);
    if(e) die("post in sa2da_body");
    // exit if start is illegal
    if(r.start<0) break;
    // process range start-end
    for(long i=r.start;i<r.end;i++) {// see sa2da()
      int_t suffixLen = getlen(d->sa[i],d->eos,d->dwords,&seqid);
      assert(seqid<=0x7FFFFFFF);     // seqid uses at most 31 bits
      assert(suffixLen>=d->lcp[i]);     // suffix length cannot be shorter than lcp
      assert(suffixLen<=d->wlen[seqid]);// suffix length cannot be larger than word length
      if(suffixLen==d->wlen[seqid]) {   // test if full word
        words++;                     
        assert(d->lcp[i]<suffixLen);    // full words are not prefix of other suffixes
      }
      if(d->lcp[i]==suffixLen) {         // save seqid + possibily extra bit 
        d->sa[i] = seqid | (1u << 31);  // mark last bit if lcp==suffix_len;
      }
      else 
        d->sa[i] = seqid;               // save only seqid = da[i]
      d->lcp[i] = suffixLen;            // save suffix length overwriting lcp
    }    
  }
  d->full_words = words;
  return NULL;
}

// transform sa[],lcp[] -> da[], suflen[] + 
// extra bit telling whether suflen[i]==lcp[i]
void sa2da(uint_t sa[], int_t lcp[], uint8_t d[], long dsize, long dwords, int w, int numt)
{
  long words=0;
  if(dwords>0x7FFFFFFF) {
    cerr << "Too many words in the dictionary. Current limit: 2^31-1\n";
    exit(1);
  }
  cout << "Converting SA and LCP Array to DA and SufLen using " << numt << " threads\n";
  
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
  if(numt==0) {
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
  }
  else { // multithread
    pthread_t t[numt];
    sa2da_data d[numt];
    pthread_mutex_t m; sem_t free_slots, data_items;
    pc_init(&free_slots,&data_items, &m);
    int index = 0;
    sarange buf[Buf_size];
    // thread creation
    for(int i=0;i<numt;i++) {
      d[i].buffer=buf;
      d[i].cindex = &index;  d[i].cindex_m = &m;
      d[i].free_slots = &free_slots; d[i].data_items = &data_items;
      d[i].full_words = 0;       
      d[i].sa = sa; d[i].eos = eos;
      d[i].lcp = lcp; d[i].wlen = wlen;
      d[i].dwords=dwords;
      pthread_create(&t[i],NULL,sa2da_body,&d[i]);
    }
    // producer code
    sarange r; int pindex = 0;
    for(long i=dwords+w+1; i<dsize; ) {
      r.start = i;
      r.end = i+Sa_block;
      if(r.end>dsize) r.end = dsize;
      // write to the buffer
      int e = sem_wait(d->free_slots);
      if(e) die("wait in sa2da");
      d->buffer[(pindex++) % Buf_size] = r; 
      e = sem_post(d->data_items);
      if(e) die("post in sa2da");
      i = r.end;
    }
    // send terminate data
    r.start = -1;
    for(int i=0;i<numt;i++) {
      int e = sem_wait(d->free_slots);
      if(e) die("wait in sa2da");
      d->buffer[(pindex++) % Buf_size] = r; 
      e = sem_post(d->data_items);
      if(e) die("post in sa2da");
    }
    // wait for termination
    for(int i=0;i<numt;i++) {
      pthread_join(t[i],NULL);
      words += d[i].full_words;
    }
    // done
    pc_destroy(&free_slots,&data_items, &m);
  }
  cout << "Conversion took " << difftime(time(NULL),start) << " wall clock seconds\n";  
  cout << "DA has size: " << dsize-dwords-w-1;
  cout << ". Dictionary words found: " << words << endl; 
}




// --------------------------------------------------------------------

// range in the suffix array of the dictionary 
typedef struct {
  long sa_start;   // starting position in the dictionary 
  long sa_end;     // end position in the dictionary
  uint8_t *bwt;    // starting position in the output bwt
  long count;      // chars to be written to the output bwt;
} sa_range;

// working data to be passed to each consumer thread
typedef struct {
  uint8_t *dict;     // dictionary
  uint_t *sa;        // suffix array for d[]
  int_t *lcp;        // lcp array
  long dsize;        // size of dict[] sa[] lcp[]  
  uint8_t *last;     // array of last symbols 
  uint32_t *ilist;   // inverted list 
  uint32_t *istart;  // starting position inside inverted list 
  long dwords;       // number of words in the dictionary 
  int w;             // window size
  sa_range  buffer[Buf_size]; // shared producer/consumer buffer 
  int cindex;                 // consumer index in buffer
  pthread_mutex_t mutex_consumers; // mutex and semaphores 
  sem_t sem_free_slots;
  sem_t sem_data_items;
  uint8_t *bwt;
  long full_words;           // output parameters, access with a mutex_consumer
  long easy_bwts; 
  long hard_bwts; 
} thread_data;

long add_sa_entries(thread_data *d, long i)
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
static size_t get_bwt_size(char *name);


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
  thread_data data;
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
#endif

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

#endif

