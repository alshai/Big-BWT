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
extern "C" {
#include "utils.h"
}


using namespace std;
using namespace __gnu_cxx;

// =============== algorithm limits =================== 
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;


// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   string parse0ext = EXTPARS0;    // extension tmp parse file 
   string parseExt =  EXTPARSE;    // extension final parse file  
   string occExt =    EXTOCC;      // extension occurrences file  
   string dictExt =   EXTDICT;     // extension dictionary file  
   string lastExt =   EXTLST;      // extension file containing last chars   
   string saExt =     EXTSAI;      // extension file containing sa info   
   int w = 10;            // sliding window size and its default 
   int p = 100;           // modulus for establishing stopping w-tuples 
   bool SAinfo = false;
};


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
void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos)
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
  if(fputc(w[w.size()- minsize-1],last)==EOF) die("Error writing to .last file");
  // compute ending position +1 of current word and write it to sa file 
  if(pos==0) pos = w.size()-1;
  else pos += w.size() -minsize; 
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
  // keep only the overlapping part of the window
  w.erase(0,w.size() - minsize);
}


// prefix free parse of file fnam. w is the window size, p is the modulus 
// use a KR-hash as the word ID that is immediately written to the parse file
void process_file(Args& arg, KR_window& krw, map<uint64_t,word_stats>& wordFreq)
{
  //open a, possibly compressed, input file
  string fnam = arg.inputFileName;
  #ifdef GZSTREAM 
  igzstream f(fnam.c_str());
  #else
  ifstream f(fnam.c_str());
  #endif    
  if(!f.rdbuf()->is_open()) {// is_open does not work on igzstreams 
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + string(fnam));
  }

  // open the 1st pass parsing file 
  string fparse = arg.inputFileName + "." + arg.parse0ext;
  // cout << "Writing tmp parsing to file " << fparse << endl;
  FILE *g = fopen(fparse.c_str(),"wb");
  if(g==NULL) die("Cannot open " + fparse);

  // open output file containing the char at position -(w+1) of each word
  string fnamelast = arg.inputFileName + "." + arg.lastExt;
  FILE *last_file = fopen(fnamelast.c_str(),"wb");
  if(last_file==NULL) die("Cannot open " + fnamelast); 
  
  // if requested open file containing the starting position of each word
  FILE *sa_file = NULL; string sa_name = "<not used>";
  if(arg.SAinfo) {
    sa_name = arg.inputFileName + "." + arg.saExt;
    sa_file = fopen(sa_name.c_str(),"wb");
    if(sa_file==NULL) die("Cannot open " + sa_name);
  } 
  
  // main loop on the chars of the input file
  int c;
  uint64_t pos = 0; // ending position +1 of current word in the original text, used for computing sa_info 
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are writte to the sa info file 
  // init first word in the parsing with a NUL char 
  string word("");
  word.append(1,Dollar);
  while( (c = f.get()) != EOF ) {
    if(c<=Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
    word.append(1,c);
    uint64_t hash = krw.addchar(c);
    if(hash%arg.p==0) {
      // end of word, save it and write its full hash to the output file
      // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" << krw.get_window() << ">" <<  endl;
      save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
    }    
  }
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg.w,Dollar);
  save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
  // close input and output files 
  if(sa_file) if(fclose(sa_file)!=0) die("Error closing "+ sa_name);
  if(fclose(last_file)!=0) die("Error closing "+ fnamelast );  
  if(fclose(g)!=0) die("Error closing "+ fparse);
  if(pos!=krw.tot_char+arg.w) cerr << "Pos: " << pos << " tot " << krw.tot_char << endl;
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
void writeDictOcc(Args &arg, map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict)
{
  assert(sortedDict.size() == wfreq.size());
  // open dictionary and occ files 
  string fdictname = arg.inputFileName + "." + arg.dictExt;
  FILE *fdict = fopen(fdictname.c_str(),"wb");
  if(fdict==NULL) die("Cannot open " + fdictname);
  string foccname = arg.inputFileName + "." + arg.occExt;
  FILE *focc = fopen(foccname.c_str(),"wb");
  if(focc==NULL) die("Cannot open " + foccname);
  
  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    size_t s = fwrite((*x).data(),1,(*x).size(), fdict);
    if(s!=(*x).size()) die("Error writing to " + fdictname);
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to " + fdictname);
    uint64_t hash = kr_hash(*x);
    auto& wf = wfreq.at(hash);
    assert(wf.occ>0);
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("Error writing to " + foccname);
    assert(wf.rank==0);
    wf.rank = wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to " + fdictname);
  if(fclose(focc)!=0) die("Error closing " + foccname);
  if(fclose(fdict)!=0) die("Error closing" + fdictname);
}

void remapParse(Args &arg, map<uint64_t,word_stats> &wfreq)
{
  // build file names
  string old_parse = arg.inputFileName + "." + arg.parse0ext;
  string new_parse = arg.inputFileName + "." + arg.parseExt;
  // recompute to double check occ
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based 
  // open parse files 
  FILE *oldp = fopen(old_parse.c_str(),"rb");
  if(oldp==NULL) die("Error opening " + old_parse);
  FILE *newp = fopen(new_parse.c_str(),"wb");
  if(newp==NULL) die("Error opening" + new_parse);
  uint64_t hash;
  while(!feof(oldp)) {
    size_t s = fread(&hash,sizeof(hash),1,oldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank = wfreq.at(hash).rank;
    occ[rank]++;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to " + new_parse);
  }
  if(fclose(newp)!=0) die("Error closing " + new_parse);
  if(fclose(oldp)!=0) die("Error closing " + old_parse);
  // check old and recomputed occ coincide 
  for(auto& x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
}
 



void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        << "\t-h  \tshow help and exit" << endl
        << "\t-s  \tcompute suffix array info" << endl;
  #ifdef GZSTREAM
  cout << "If the input file is gzipped it is automatically extracted\n";
  #endif
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
   while ((c = getopt( argc, argv, "p:w:sh") ) != -1) {
      switch(c) {
        case 's':
        arg.SAinfo = true; break;
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'p':
        sarg.assign( optarg );
        arg.p = stoi( sarg ); break;
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
}



int main(int argc, char** argv)
{
  
  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;  

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;  
  // init window-based karp-rabin fingerprint
  KR_window krw(arg.w); // input is window size and alphabet 
  // init sorted map counting the number of occurrences of each word
  map<uint64_t,word_stats> wordFreq;  

  // ------------ parsing input file 
  try {
      process_file(arg,krw,wordFreq);
  }
  catch(const std::bad_alloc&) {
      cout << "Out of memory... emergency exit\n";
      die("bad alloc exception");
  }
  // first report 
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

  // -------------- second pass  
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
  writeDictOcc(arg, wordFreq, dictArray);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
    
  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file\n";
  remapParse(arg, wordFreq);
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";        
  return 0;
}

