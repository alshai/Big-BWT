extern "C" {
#include "xerrors.h"
}
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

// struct returned by mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input 
  long start, end; // input
  long skipped, parsed, words;  // output
  FILE *parse, *last, *sa;
} mt_data;


void *mt_parse(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d->wordFreq;

  if(arg->verbose>1)
    printf("Scanning from %ld, size %ld\n",d->start,d->end-d->start);

  // open input file 
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }

  // prepare for parsing 
  f.seekg(d->start); // move to the begining of assigned region
  KR_window krw(arg->w);
  int c; string word = ""; 
  d->skipped = d->parsed = d->words = 0;
  if(d->start==0) {
    word.append(1,Dollar);// no need to reach the next kr-window 
  }
  else {   // reach the next breaking point 
    while( (c = f.get()) != EOF ) {
      if(c<=Dollar) die("Invalid char found in input file. Exiting...");
      d->skipped++;
      if(d->start + d->skipped == d->end + arg->w) {f.close(); return NULL;} 
      word.append(1,c);
      uint64_t hash = krw.addchar(c);
      if(hash%arg->p==0 && d->skipped >= arg->w) break;
    }
    if(c==EOF) {f.close(); return NULL;} // reached EOF without finding a breaking point nothing to do   
    d->parsed = arg->w;   // the kr-window is part of the next word
    d->skipped -= arg->w; // ... so w less chars have been skipped
    word.erase(0,word.size() - arg->w);// keep only the last w chars 
  }
  // cout << "Skipped: " << d->skipped << endl;
  
  // there is some parsing to do: open output files 
  d->parse = tmpfile();
  d->last = tmpfile();
  if(arg->SAinfo) d->sa = tmpfile();
  else d->sa = NULL;
  
  // do the parsing:
  uint64_t pos = d->start+d->skipped;  // starting position in text of current word
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file 
  while( (c = f.get()) != EOF ) {
    if(c<=Dollar) die("Invalid char found in input file. Exiting...");
    word.append(1,c);
    uint64_t hash = krw.addchar(c);
    d->parsed++;
    if(hash%arg->p==0 && d->parsed>arg->w) {
      // end of word, save it and write its full hash to the output file
      save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
      d->words++;
      if(d->start+d->skipped+d->parsed>=d->end+arg->w) {f.close(); return NULL;}
    }    
  }
  // end of file reached 
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg->w,Dollar);
  save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
  // close input file and return 
  f.close();
  return NULL;
}


// prefix free parse of file fnam. w is the window size, p is the modulus 
// use a KR-hash as the word ID that is written to the parse file
uint64_t mt_process_file(Args& arg, map<uint64_t,word_stats>& wf)
{
  // get input file size 
  ifstream f(arg.inputFileName, std::ifstream::ate);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " +arg.inputFileName);
  }
  long size = f.tellg();
  f.close();   

  // prepare and execute threads 
  assert(arg.th>0);
  pthread_t t[arg.th];
  mt_data td[arg.th];
  for(int i=0;i<arg.th;i++) {
    td[i].wordFreq = &wf;
    td[i].arg = &arg;
    td[i].start = i*(size/arg.th); // range start
    td[i].end = (i+1==arg.th) ? size : (i+1)*(size/arg.th); // range end
    assert(td[i].end<=size);
    xpthread_create(&t[i],NULL,&mt_parse,&td[i],__LINE__,__FILE__);
  }

  // open the 1st pass parsing file 
  FILE *parse = open_aux_file(arg.inputFileName.c_str(),arg.parse0ext.c_str(),"wb");
  // open output file containing the char at position -(w+1) of each word
  FILE *last = open_aux_file(arg.inputFileName.c_str(),arg.lastExt.c_str(),"wb");  
  // if requested open file containing the ending position+1 of each word
  FILE *sa = NULL;
  if(arg.SAinfo) 
    sa = open_aux_file(arg.inputFileName.c_str(),arg.saExt.c_str(),"wb");
  
  // wait for the threads to finish (in order) and copy data to output files
  long tot_char=0;
  for(int i=0;i<arg.th;i++) {
    char buf[BUFSIZ]; size_t size;
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
    if(arg.verbose) {
      cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:";
      cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" << td[i].words << endl;
    }
    if(td[i].words>0) {
      // copy parse
      rewind(td[i].parse);
      while ((size = fread(buf, 1, BUFSIZ, td[i].parse)))
        fwrite(buf, 1, size, parse);
      fclose(td[i].parse);
      // copy last
      rewind(td[i].last);
      while ((size = fread(buf, 1, BUFSIZ, td[i].last)))
        fwrite(buf, 1, size, last);
      fclose(td[i].last);
      if(arg.SAinfo) {
      // copy parse
      rewind(td[i].sa);
        while ((size = fread(buf, 1, BUFSIZ, td[i].sa)))
          fwrite(buf, 1, size, sa);
        fclose(td[i].sa);
      }
      // extra check
      assert(td[i].parsed>arg.w);
      tot_char += td[i].parsed - (i!=0? arg.w: 0); //parsed - overlapping 
    }
    else assert(i>0); // the first thread must produce some words
  }
  assert(tot_char==size);
  // close output files 
  if(sa) if(fclose(sa)!=0) die("Error closing SA file");
  if(fclose(last)!=0) die("Error closing last file");  
  if(fclose(parse)!=0) die("Error closing parse file");
  return size;   
}


