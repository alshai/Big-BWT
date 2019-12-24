#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <getopt.h>
#include "pfbwt.hpp"
// extern "C" {
// #include "utils.h"
// }

struct Args {
    std::string prefix;
    int w = 10;
    bool sa = false;
};

void print_help() {
    return;
}

Args parseArgs(int argc, char** argv) {
    Args args;
    int c;

    // fputs("==== Command line:", stderr);
    // for(int i=0;i<argc;i++)
    //     fprintf(stderr, " %s",argv[i]);
    // fputs("\n", stderr);

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
    FILE* bwt_fp = open_aux_file(args.prefix.data(),"bwt","wb");
    p.generate_bwt_lcp([&bwt_fp](uint8_t c) { fputc(c, bwt_fp); });
    fclose(bwt_fp);
}

int main(int argc, char** argv) {
    Args args(parseArgs(argc, argv));
    run_pfbwt(args);
}
