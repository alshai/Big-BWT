
// special symbols used by the construction algorithm:
//   they cannot appear in the input file 
//   the 0 symbol is used in the final BWT file as the EOF char  

#define Dollar 2     // special char for the parsing algorithm, must be the highest special char 
#define EndOfWord 1  // word delimiter for the plain dictionary file
#define EndOfDict 0  // end of dictionary delimiter 

void die(const char *s);
FILE *open_aux_file(const char *base, const char *ext, const char *mode);
int fd_open_aux_file(const char *base, const char *ext, int mode);
