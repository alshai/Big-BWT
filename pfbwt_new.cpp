#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <getopt.h>
#include "pfbwt.hpp"

struct Args {
    std::string prefix;
    int w = 64;
    bool sa = false;
};

void print_help() {
    return;
}

Args parseArgs(int argc, char** argv) {
    Args args;
    int c;

    fputs("==== Command line:", stderr);
    for(int i=0;i<argc;i++)
        fprintf(stderr, " %s",argv[i]);
    fputs("\n", stderr);

    std::string sarg;
    while ((c = getopt( argc, argv, "w:hsf") ) != -1) {
        switch(c) {
            case 'f': // legacy
                break;
            case 's':
                args.sa = true; break;
            case 'w':
                args.w = atoi(optarg); break;
            case 'h':
                print_help(); exit(1);
            case '?':
                fprintf(stderr, "Unknown option. Use -h for help.\n");
                exit(1);
        }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
        args.prefix.assign( argv[optind] );
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        print_help();
        exit(1);
    }
    return args;
}

void run_pfbwt(Args args) {
    pfbwt p(args.prefix, args.w); // load dict, ilist, last, etc
    p.generate_bwt_lcp(args.prefix);
}

int main(int argc, char** argv) {
    // parse arguments here
    Args args(parseArgs(argc, argv));
    run_pfbwt(args);
}
