#include <cstdio>
#include <chrono>
#include <getopt.h>
#include "hash.hpp"
#include "parse.hpp"

struct Timer {
    using clock = std::chrono::system_clock;
    using sec = std::chrono::duration<double>;

    Timer(std::string m) : msg(m) {
        start = clock::now();
    }

    ~Timer() {
        sec dur = (clock::now() - start);
        fprintf(stderr, "%s%.2fs\n", msg.data(), static_cast<double>(dur.count()));
    }

    std::chrono::time_point<clock> start;
    std::string msg;
};


uint8_t seq_nt4_ntoa_table[] = {
    /*   0 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /*  16 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /*  32 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5,
           /*                                        - */
    /*  48 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /*  64 */ 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 0, 5,
           /*    A  B  C  D        G  H        K     M  N */
    /*  80 */ 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
           /*       R  S  T     V  W  X  Y */
    /*  96 */ 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 0, 5,
           /*    a  b  c  d        g  h        k     m  n */
    /* 112 */ 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
           /*       r  s  t     v  w  x  y */
    /* 128 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 144 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 160 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 176 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 192 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 208 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 224 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    /* 240 */ 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
};

void print_help() { 
    fprintf(stderr, "Usage:\n\t./parse [options] <fasta file>\n");
    fprintf(stderr, "Options:\n-w\twindow size\n-p\tstop word modulus\n-s\tprint SA info\n");
}

Args parseArgs(int argc, char** argv) {
    Args args;
    int c;
    extern char *optarg;
    extern int optind;

    fputs("==== Command line:", stderr);
    for(int i=0;i<argc;i++)
        fprintf(stderr, " %s",argv[i]);
    fputs("\n", stderr);

    std::string sarg;
    while ((c = getopt( argc, argv, "p:w:hvsf") ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 'w':
                sarg.assign( optarg );
                args.w = stoi( sarg ); break;
            case 'p':
                sarg.assign( optarg );
                args.p = stoi( sarg ); break;
            case 's':
                args.sai = true; break;
            case 'h':
                print_help(); exit(1);
            case '?':
                fprintf(stderr, "Unknown option. Use -h for help.\n");
                exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
        args.in_fname.assign( argv[optind] );
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        print_help();
        exit(1);
    }
    // check algorithm parameters
    if ((args.w < 4) || (args.w > 32)) {
        fprintf(stderr, "Windows size must be between 4 and 31 (inclusive)\n");
        exit(1);
    }
    if (args.p<4) {
        fprintf(stderr, "Modulus must be at least 4\n");
        exit(1);
    }
    return args;
}


int main(int argc, char** argv) {
    Args args(parseArgs(argc, argv));
    // build the dictionary and populate .last, .sai and .parse_old
    Parser<WangHash> p(args.w, args.p);
    {
        Timer t("Task\tParsing\t");
        p.parse_fasta(args.in_fname.data(), args.sai);
    }
    {
        Timer t("Task\tSorting/Update\t");
        p.update_dict(args.in_fname.data());
    }
    {
        Timer t("Task\tParseDump\t");
        p.dump_parse(args.in_fname.data());
    }
    return 0;
}
