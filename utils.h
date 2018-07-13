
// special symbols used by the construction algorithm:
//   they cannot appear in the input file 
//   the 0 symbol is used in the final BWT file as the EOF char  

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char 
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter 
#define IBYTES 5     // bytes used to represent a large integer (at most 8)

// file name extensions
#define EXTPARSE "parse"
#define EXTPARS0 "parse_old"
#define EXTOCC   "occ"
#define EXTDICT  "dict"
#define EXTLST   "last"
#define EXTBWLST "bwlast"
#define EXTSAI   "sai"
#define EXTBWSAI "bwsai"
#define EXTILIST "ilist"


void die(const char *s);
FILE *open_aux_file(const char *base, const char *ext, const char *mode);
int fd_open_aux_file(const char *base, const char *ext, int mode);
uint64_t get_myint(uint8_t *a, long n, long i);
void write_myint(uint64_t u, FILE *f);
void get_and_write_myint(uint8_t *a, long n, long i, FILE *f);
