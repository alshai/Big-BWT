/* ******************************************************************************
 * newscan.cpp
 * 
 * parsing algorithm for the bwt construction for repetitive sequences based 
 * on prefix free parsing
 * 
 * Usage:
 *   newscan.x wsize modulus file
 * 
 * Accept any kind of file that does not contain the chars 0x0, 0x1, 0x2 
 * which are used internally. If input file is gzipped use cnewscan.x which 
 * automatcally extracts the content
 * 
 * The parameters wsize and modulus are used to define the prefix free parsing 
 * using KR-fingerprints
 * 
 * The algorithm computes the prefix free parsing of 
 *     T = (0x2)file_content(0x2)^wsize
 * in a dictionary of words D and a parsing P of T in terms of the  
 * dictionary words. Note that the words in the parsing overlap by wsize.
 * Let d denote the number of words in D and p the number of phrases in 
 * the parsing P
 * 
 * newscan outputs the following files:
 * 
 * file.dict
 * containing the dictionary words in lexicogaphic order with a 0x1 at the end of
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
 * 
 * The output of newscan must be processed by bwtparse, that, invoked as
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
#ifdef GZSTREAM
#include <gzstream.h>
#endif


using namespace std;
using namespace __gnu_cxx;

#define Dollar    0x2  // special char for the parsing algorithm, must be the highest special char 
#define EndOfWord 0x1  // word delimiter for the plain dictionary file
#define EndOfDict 0x0  // end of dictionary char 

// =============== algorithm limits =================== 
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

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


// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};


void die(const string s)
{
  perror(s.c_str());
  exit(1);
}


// compute 64-bit KR hash of a string 
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



// save current word in the freq map and update it leaving only the 
// last minsize chars which is the overlap with next word  
void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, FILE *last)
{
  if(w.size() < minsize+1) return;
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");
  
  // update frequency table for current hash
  if(freq.find(hash)==freq.end()) {
      freq[hash].occ = 1; // new hash
      freq[hash].str = w; 
  }
  else {
      freq[hash].occ += 1; // known hash
      if(freq[hash].occ <=0) {
        cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
        cerr<< MAX_WORD_OCC << ") exceeded\n";
        exit(1);
      }
      if(freq[hash].str != w) {
        cerr << "Emergency exit! Hash collision for strings:\n";
        cerr << freq[hash].str << "\n  vs\n" <<  w << endl;
        exit(1);
      }
  }
  // output char w+1 from the end
  if(fputc(w[w.size()- minsize-1],last)==EOF) die(".last write");
  // keep only the overlapping part of the window
  w.erase(0,w.size() - minsize);
}


// prefix free parse of file fnam. w is the window size, p is the modulus 
// use a KR-hash as the word ID that is immediately written to the parse file
void process_file(int w, int p, string fnam, KR_window& krw, map<uint64_t,word_stats>& wordFreq, FILE *tmp_parse_file)
{
  //open a, possibly compressed input file
  #ifdef GZSTREAM 
  igzstream f(fnam.c_str());
  #else
  ifstream f(fnam.c_str());
  #endif  
  
  if(!f.rdbuf()->is_open()) {// is_open does not work on igzstreams 
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + string(fnam));
  }

  // open file that will contain the char at position -(w+1) of each word
  string fnamelast = fnam + ".last";
  FILE *last_file = fopen(fnamelast.c_str(),"wb");
  if(last_file==NULL) die(".last file open"); 
  
  // main loop on the chars of the input file
  int c;
  // init first word in the parsing with a NUL char 
  string word("");
  word.append(1,Dollar);
  while( (c = f.get()) != EOF ) {
    if(c<=Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
    word.append(1,c);
    uint64_t hash = krw.addchar(c);
    if(hash%p==0) {
      // end of word, save it and write its full hash to the output file
      // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" << krw.get_window() << ">" <<  endl;
      save_update_word(word,w,wordFreq,tmp_parse_file,last_file);
    }    
  }
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(w,Dollar);
  save_update_word(word,w,wordFreq,tmp_parse_file,last_file);
  if(fclose(last_file)!=0) die(".last close");  
  f.close();      
}

void old_style_report(map<uint64_t,word_stats>& wordFreq, uint64_t totChar,string fname,bool writeDict)
{  
  // init vectors for word statistics
  vector<int> lenFreq(10000,0); // lenFreq[i] = # distinct words of length i
  vector<int> lenTotw(10000,0); // lenTotw[i] =  # words of length i 
  
  // open dict file if requested
  FILE *f;
  if(writeDict) {
    string fdict = fname + string(".dict_old");
    cout << "Writing dictionary to file " << fdict << endl;
    f = fopen(fdict.c_str(),"wb");
    if(f==NULL) die("Open dictionary file");
  }
  else f=NULL;
  
  // process words in hash order: compute statistics
  uint64_t newLen = 0, totw = 0;
  int words_to_print = 0;   //  change this to send more or less words to stderr
  for (auto& x: wordFreq) {
    // select which words we want to be sent to stderr 
    if(words_to_print>0 && x.second.str.size()> 0 && x.second.occ>0) { 
      words_to_print--; 
      cerr << "~"<< x.second.str << "~ occ: " << x.second.occ << endl; 
    }
    uint32_t wlen = x.second.str.size();      // length of current word
    newLen += wlen;                          // update total len 
    if(wlen>= lenFreq.size()) {     
      lenFreq.resize(wlen+100,0);     // resize length-related arrays  
      lenTotw.resize(wlen+100,0);     
    }     
    lenFreq[wlen] += 1;              // one more distinct word of length wlen
    lenTotw[wlen] += x.second.occ;   // lenFreq[s] more occ's of a word of length wlen
    totw += x.second.occ;            // count total number of words
    if(f!=NULL) {
      if(fwrite(&x.first,sizeof(x.first),1,f)!=1) die("dict writing 0");    
      if(fwrite(&x.second.occ,sizeof(x.second.occ),1,f)!=1) die("dict writing 1");
      if(fwrite(&wlen,sizeof(wlen),1,f)!=1) die("dict writing 2");
      if(fwrite(x.second.str.c_str(),1,wlen,f)!=wlen) die("dict writing 3"); 
    }
  } 
  if(f!=NULL)
    if(fclose(f)!=0) die("dict close"); 
  cout << "Sum of lengths of distinct words: " << newLen <<" ("<<100.0*newLen/totChar<<"%)"<<endl;
  cout << "Total number of words: " << totw << ". Average word length: " << totChar/totw << endl;

  // output detailed statistics on word lengths
  printf("%8s  %10s %10s\t\t   %s\n","Length", "TotalWords", "DistWords", "Ave. # Occs");
  for(size_t i=0;i<lenFreq.size();i++) {
    if(lenFreq[i]!=0) {
      printf("%8zu: %10d %10d\t\t   %.2lf\n",i, lenTotw[i], lenFreq[i], 1.0*lenTotw[i]/lenFreq[i]);
      totw += i*lenFreq[i];
    }
  }
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict, string fname)
{
  assert(sortedDict.size() == wfreq.size());
  // open dictionary and occ files 
  string fdictname = fname + ".dict";
  FILE *fdict = fopen(fdictname.c_str(),"wb");
  if(fdict==NULL) die("dictionary open");
  string foccname = fname + ".occ";
  FILE *focc = fopen(foccname.c_str(),"wb");
  if(focc==NULL) die("occ file open");
  
  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    size_t s = fwrite((*x).data(),1,(*x).size(), fdict);
    if(s!=(*x).size()) die("dictionary word write");
    if(fputc(EndOfWord,fdict)==EOF) die("dictionary EndOfWord write");
    uint64_t hash = kr_hash(*x);
    auto& wf = wfreq.at(hash);
    assert(wf.occ>0);
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("occfile write");
    assert(wf.rank==0);
    wf.rank = wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("dictionary EndOfDict write");
  if(fclose(focc)!=0) die("occfile close");
  if(fclose(fdict)!=0) die("dictionary close");
}

void remapParse(map<uint64_t,word_stats> &wfreq, string old_parse, string new_parse)
{
  // double check occ
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based 
  // open parse files 
  FILE *oldp = fopen(old_parse.c_str(),"rb");
  if(oldp==NULL) die("old parse open");
  FILE *newp = fopen(new_parse.c_str(),"wb");
  if(newp==NULL) die("new parse open");
  uint64_t hash;
  while(!feof(oldp)) {
    size_t s = fread(&hash,sizeof(hash),1,oldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank = wfreq.at(hash).rank;
    occ[rank]++;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("New parse write");
  }
  if(fclose(newp)!=0) die("new parse close");
  if(fclose(oldp)!=0) die("old parse close");
  for(auto& x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
}


int main(int argc, char** argv)
{
  // check command line
  if(argc!=4) {
    cerr << "Usage:\n\t";
    cerr << argv[0] << " wsize modulus file\n" << endl;
    cerr << "Input file can be anything not containing chars 0x00 0x01 0x02\n";
    #ifdef GZSTREAM
    cerr << "If the input file is gzipped it is automatically extracted\n";
    #endif
    cerr <<endl;
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  time_t start_wc = time(NULL);
  
  // translate command line parameters
  int w = atoi(argv[1]);                 // sliding window size 
  int p = atoi(argv[2]);                 // modulus for stopping w-tuples 
  assert(w>1 && p>1);
  cout << "Windows size: " << w << endl;
  cout << "Stop word modulus: " << p << endl;  

  // init window-based karp-rabin fingerprint
  KR_window krw(w); // input is window size and alphabet 
  // init sorted map counting the number of occurrences of each word
  map<uint64_t,word_stats> wordFreq;  
  // open the parsing file 
  string fname(argv[3]);
  string fparse = fname + ".parse_old";
  cout << "Writing parsing to file " << fparse << endl;
  FILE *g = fopen(fparse.c_str(),"wb");
  if(g==NULL) die("Open parse file");

  // parsing input file 
  try {
      process_file(w,p,fname,krw,wordFreq, g);
  }
  catch(const std::bad_alloc&) {
      cout << "Out of memory... emergency exit\n";
      die("bad alloc exception");
  }
  if(fclose(g)!=0) die("parse_old close");  
  uint64_t totChar = krw.tot_char;
  uint64_t totDWord = wordFreq.size();
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" <<endl;
  cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  // check # distinct words
  if(totDWord>MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinc words (" << totDWord << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }

  // report statistics and, if last parameter is true, write old style dictionary 
  // old_style_report(wordFreq,totChar,fname,true);
  
  start_wc = time(NULL);
  // create array of dictionary words
  vector<const string *> dictArray;
  dictArray.reserve(totDWord);
  // fill array
  uint64_t sumLen = 0;
  uint64_t totWord = 0;
  for (auto& x: wordFreq) {
    sumLen += x.second.str.size();
    totWord += x.second.occ;
    dictArray.push_back(&x.second.str);
  }
  assert(dictArray.size()==totDWord);
  cout << "Sum of lenghts of dictionary words: " << sumLen << endl; 
  cout << "Total number of words: " << totWord << endl; 
  // sort dictionary
  sort(dictArray.begin(), dictArray.end(),pstringCompare);
  // write plain dictionary and occ file, also compute rank for each hash 
  cout << "Writing plain dictionary and occ file\n";
  writeDictOcc(wordFreq, dictArray,fname);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
    
  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file\n";
  remapParse(wordFreq, fparse, fname+".parse");
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  
  return 0;
}

