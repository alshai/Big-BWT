/* ******************************************************************************
 * pscan.cpp
 * 
 * paralell parsing algorithm for bwt construction of repetitive sequences based 
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/)
 * 
 * Usage:
 *   pscan.x wsize modulus file
 * 
 * Accepts any kind of file that does not contain the chars 0x0, 0x1, 0x2 
 * which are used internally. If input file is gzipped use pscan.x which 
 * automatically extracts the content
 * 
 * The parameters wsize and modulus are used to define the prefix free parsing 
 * using KR-fingerprints (see paper)
 * 
 * The algorithm computes the prefix free parsing of 
 *     T = (0x2)file_content(0x2)^wsize
 * cresting a dictionary of words D and a parsing P of T in terms of the  
 * dictionary words. Consecutive words in the parsing overlap by wsize.
 *
 * Let d denote the number of words in D and p the number of phrases in 
 * the parsing P
 * 
 * pscan outputs the following files:
 * 
 * file.dict
 * containing the dictionary words in lexicographic order with a 0x1 at the end of
 * each word and a 0x0 at the end of the file. Size: |D| + d + 1 where
 * |D| is the sum of the word lengths
 * 
 * file.occ
 * the number of occurrences of each word in lexicographic order.
 * We assume the number of occurrences of each word is at most 2^32-1
 * so the size is 4d bytes
 * 
 * file.parse
 * containing the parse P with each word identified with its 1-based lexicographic 
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes
 * 
 * file.last 
 * contaning the charater in positon w+1 from the end for each dictionary word
 * Size: d
 * 
 * file.sai (if option -s is given on the command line) 
 * containing the ending position +1 of each dictionary word in the original
 * text written using IBYTES bytes for each entry (IBYTES defined in utils.h)
 * Size: d*IBYTES
 * 
 * The output of pscan must be processed by bwtparse, which invoked as
 * 
 *    bwtparse file
 * 
 * computes the BWT of file.parse and produces file.ilist of size 4p+4 bytes
 * contaning, for each dictionary word in lexicographic order, the list 
 * of BWT positions where that word appears (ie i\in ilist(w) <=> BWT[i]=w).
 * There is also an entry for the EOF word which is not in the dictionary 
 * but is assumed to be the smallest word.  
 * 
 * In addition, bwtparse permutes file.last according to
 * the BWT permutation and generates file.bwlast such that file.bwlast[i] 
 * is the char from P[SA[i]-2] (if SA[i]==0 , BWT[i]=0 and file.bwlast[i]=0, 
 * if SA[i]==1, BWT[i]=P[0] and file.bwlast[i] is taken from P[n-1], the last 
 * word in the parsing).  
 * 
 * If the option -s is given to bwtparse, it permutes file.sai according
 * to the BWT permutation and generate file.bwsai using again IBYTES
 * per entry.  file.bwsai[i] is the ending position+1 of BWT[i] in the 
 * original text 
 * 
 * The output of bwtparse (the files .ilist .bwlast) together with the
 * dictionary itself (file .dict) and the number of occurrences
 * of each word (file .occ) are used to compute the final BWT by the 
 * pfbwt algorithm.
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
#include "utils.h"
#include "xerrors.h"
}

using namespace std;
//using namespace __gnu_cxx;



// =============== algorithm limits =================== 
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   int w = 10;            // sliding window size and its default 
   int p = 100;           // modulus for establishing stopping w-tuples 
   bool SAinfo = false;   // compute SA information
   bool compress = false; // parsing called in compress mode 
   int th=4;              // number of helper threads
   int verbose=0;         // verbosity level
   FILE *tmp_parse_file, *last_file, *sa_file; 
};

// -----------------------------------------------------------
// struct containing the maps and the relative mutex
struct MTmaps {
   int mt_ratio = 3;                       // ratio between #maps and #threads 
   int n;                                  // number of maps 
   vector<map<uint64_t,word_stats>> maps;  // maps
   pthread_mutex_t *muts;                  // mutex for each map
   
   // constructor
   MTmaps(int numthreads) {
     // init number of maps
     n = mt_ratio*numthreads;
     // init maps
     maps.resize(n);
     // init mutexes
     muts = new pthread_mutex_t[n];
     for(int i=0;i<n;i++) 
       xpthread_mutex_init(&muts[i], NULL, __LINE__,__FILE__);
   }
   
   // destructor
   ~MTmaps() {
     delete[] muts;
   }
   
   // return the total size of the maps, ie total number of stored words
   uint64_t size() {
     uint64_t s=0;
     for(int i=0;i<n;i++)
       s += maps[i].size();
     return s;
   }
   
   // return the rank of the string associated to hash h
   word_int_t rank(uint64_t h) {
     return maps[h%n].at(h).rank;
   }
   
  // add the association hash->w to the map hash%n
  // using a mutex for exclusive write  
  void update(uint64_t hash, string &w);
      
};

void MTmaps::update(uint64_t hash, string &w) 
{  
  int i = hash % n;
  map<uint64_t,word_stats> *freq = &maps[i];
  pthread_mutex_t *m = &muts[i]; 
  xpthread_mutex_lock(m,__LINE__,__FILE__);
  // update frequency table for current hash
  if(freq->find(hash)==freq->end()) {
      (*freq)[hash].occ = 1; // new hash
      (*freq)[hash].str = w; 
  }
  else {
      word_stats *wfreq = &(*freq)[hash];  // pointer to the stats for w
      wfreq->occ += 1; // known hash
      if(wfreq->occ <=0) {
        cerr << "Emergency exit! Maximum # of occurences of dictionary word (";
        cerr<< MAX_WORD_OCC << ") exceeded\n";
        exit(1);
      }
      if(wfreq->str != w) {
        cerr << "Emergency exit! Hash collision for strings:\n";
        cerr << wfreq->str << "\n  vs\n" <<  w << endl;
        exit(1);
      }
  }
  xpthread_mutex_unlock(m,__LINE__,__FILE__);
}





// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int *window;
  int asize;
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime 
  
  KR_window(int w): wsize(w) {
    asize = 256;
    asize_pot = 1;
    for(int i=1;i<wsize;i++) 
      asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm  
    // alloc and clear window
    window = new int[wsize];
    reset();     
  }
  
  // init window, hash, and tot_char 
  void reset() {
    for(int i=0;i<wsize;i++) window[i]=0;
    // init hash value and related values
    hash=tot_char=0;    
  }
  
  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers 
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution  
    hash = (asize*hash + c) % prime;      //  add char i 
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash; 
  }
  // debug only 
  string get_window() {
    string w = "";
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }
  
  ~KR_window() {
    delete[] window;
  } 

};

// -----------------------------------------------------------


// compute 64-bit KR hash of a string 
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
// if collisions occur use a prime close to 2**63 and 128 bit variables 
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      assert(c>=0 && c< 256);
      hash = (256*hash + c) % prime;    //  add char k
    } 
    return hash; 
}


#include "pscan.hpp"



// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, MTmaps &mtmaps, vector<const string *> &sortedDict)
{
  FILE *fdict;
  // open dictionary and occ files
  if(arg.compress)
    fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICZ,"wb");
  else
    fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICT,"wb");
  FILE *focc = open_aux_file(arg.inputFileName.c_str(),EXTOCC,"wb");
  
  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    const char *word = (*x).data();       // current dictionary word
    int offset=0; size_t len = (*x).size();  // offset and length of word
    assert(len>(size_t)arg.w);
    if(arg.compress) {  // if we are compressing remove overlapping and extraneous chars
      len -= arg.w;     // remove the last w chars 
      if(word[0]==Dollar) {offset=1; len -= 1;} // remove the very first Dollar
    }
    size_t s = fwrite(word+offset,1,len, fdict);
    if(s!=len) die("Error writing to DICT file");
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to DICT file");
    uint64_t hash = kr_hash(*x);
    auto& wf = (mtmaps.maps[hash%mtmaps.n]).at(hash);
    assert(wf.occ>0);
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("Error writing to OCC file");
    assert(wf.rank==0);
    wf.rank = wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to DICT file");
  if(fclose(focc)!=0) die("Error closing OCC file");
  if(fclose(fdict)!=0) die("Error closing DICT file");
}

void remapParse(Args &arg, MTmaps &mtmaps)
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
  FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");

  // recompute occ as an extra check 
  vector<occ_int_t> occ(mtmaps.size()+1,0); // ranks are zero based 
  uint64_t hash;
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank = mtmaps.rank(hash);
    occ[rank]++;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to new parse file");
  }
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  // check old and recomputed occ coincide
  for(auto &m : mtmaps.maps) 
    for(auto& x : m)
      assert(x.second.occ == occ[x.second.rank]);
}
 



void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        << "\t-t M\tnumber of helper threads, def. 4 " << endl
        << "\t-h  \tshow help and exit" << endl
        << "\t-s  \tcompute suffix array info" << endl;
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
   int c;
   extern char *optarg;
   extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

   string sarg;
   while ((c = getopt( argc, argv, "p:w:sht:vc") ) != -1) {
      switch(c) {
        case 's':
        arg.SAinfo = true; break;
        case 'c':
        arg.compress = true; break;
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'p':
        sarg.assign( optarg );
        arg.p = stoi( sarg ); break;
        case 't':
        sarg.assign( optarg );
        arg.th = stoi( sarg ); break;
        case 'v':
           arg.verbose++; break;
        case 'h':
           print_help(argv, arg); exit(1);
        case '?':
        cout << "Unknown option. Use -h for help." << endl;
        exit(1);
      }
   }
   // the only input parameter is the file name 
   if (argc == optind+1) {
     arg.inputFileName.assign( argv[optind] );
   }
   else {
      cout << "Invalid number of arguments" << endl;
      print_help(argv,arg);
   }
   // check algorithm parameters 
   if(arg.w <4) {
     cout << "Windows size must be at least 4\n";
     exit(1);
   }
   if(arg.p<10) {
     cout << "Modulus must be at leas 10\n";
     exit(1);
   }
   if(arg.th<=0) {
     cout << "There must be at least one helper thread\n";
     exit(1);
   }
}



int main(int argc, char** argv)
{
  // translate command line parameters and store them to arg 
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;  
  
  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;
  // init multithread maps 
  MTmaps mtmaps(arg.th);
  uint64_t totChar;

  // ------------ parsing input file 
  try {
    totChar = mt_process_file(arg,mtmaps);
  }
  catch(const std::bad_alloc&) {
      cout << "Out of memory (parsing phase)... emergency exit\n";
      die("bad alloc exception");
  }
  // first report 
  uint64_t totDWord = mtmaps.size();
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" <<endl;
  cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  // check # distinct words
  if(totDWord>MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinct words (" << totDWord << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }

  // -------------- second pass  
  start_wc = time(NULL);
  // create array of dictionary words
  vector<const string *> dictArray;
  dictArray.reserve(totDWord);
  // fill array
  uint64_t sumLen = 0;
  uint64_t totWord = 0;
  
  // copy words from all maps to the dictionary
  for(auto& wordFreq: mtmaps.maps) {
    for(auto& x: wordFreq) {
      sumLen += x.second.str.size();
      totWord += x.second.occ;
      dictArray.push_back(&x.second.str);
    }
  }
  assert(dictArray.size()==totDWord);
  cout << "Sum of lenghts of dictionary words: " << sumLen << endl; 
  cout << "Total number of words: " << totWord << endl; 

  // sort dictionary
  sort(dictArray.begin(), dictArray.end(),pstringCompare);
  // write plain dictionary and occ file, also compute rank for each hash 
  cout << "Writing plain dictionary and occ file\n";
  writeDictOcc(arg, mtmaps, dictArray);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
    
  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file\n";
  remapParse(arg, mtmaps);
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";        
  return 0;
}
