#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <string>
#include <getopt.h>
#include <functional>
#include <vector>
#include <algorithm>
#include <map>
#include <zlib.h>
#include <kseq.h>
KSEQ_INIT(gzFile, gzread);
#include "util.h"
#include "hash.h"
#include "rollinghashcpp/rabinkarphash.h"
#include "rollinghashcpp/cyclichash.h"

struct Args {
    std::string in_fname;
    size_t w = 10;
    size_t p = 100;
    bool sai = false;
};

void print_help(char** argc, Args args) {
    (void) argc;
    (void) args;
    fprintf(stderr, "help msg here\n");
}

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

/*
 * s: string to parse
 * l: length of string
 * k: window size
 * p: modulo to determine candidate windows
 * f: user-defined function to handle window
 * comments:
 *  Ns are treated as As. Anticipate weirdness if significant amounts of
 *  non-consecutive Ns appear in the input
 */
template<typename Func>
std::tuple<size_t, size_t> select_windows_2bit_wang(const char* s, const size_t l, const size_t k, const size_t p, Func f) {
    uint64_t kmer = 0;
    uint64_t mask = ((1ULL << 2 * k) - 1);
    // load the first kmer, then continue;
    for (size_t i = 0; i < k; ++i) {
        if (s[i] <= Dollar) die("invalid character found in input!\n");
        int c =  (int) seq_nt4_ntoa_table[(size_t) s[i]];
        // disallow non-ACGTN so we can pack in 2 bits
        if (c > 3) { fprintf(stderr, "unrecognized character %c in input! aborting...\n", s[i]); exit(1); }
        kmer = (kmer << 2 | c) & mask;
    } // we don't bother hashing the first kmer, since it's the start of a phrase anyway
    size_t start = 0;
    for (size_t i = 1; i < l - k + 1; ++i) {
        if (s[i+k-1] <= Dollar) die("invalid character found in input!\n");
        int c = (int) seq_nt4_ntoa_table[(size_t) s[i+k-1]];
        if (c > 3) { fprintf(stderr, "unrecognized character %c in input! aborting...\n", s[i+k-1]); exit(1); }
        kmer = (kmer << 2 | c) & mask;
        uint64_t kmer_hash = wang_hash(kmer);
        if (kmer_hash % p == 0) {
            f(start, i+k-start); // user can do what they please with a phrase boundary
            start = i;
        }
    }
    // deal with the last phrase
    return std::make_tuple(start, l-start);
}

template<typename Func, typename Hasher=CyclicHash<uint32_t>>
std::tuple<size_t, size_t> select_windows_rolling(const char* s, const size_t l, const size_t k, const size_t p, Func f) {
    Hasher hf(k);
    for (size_t i = 0; i < k; ++i) {
        hf.eat(s[i]);
    } // we don't bother hashing the first kmer, since it's the start of a phrase anyway
    size_t start = 0;
    for (size_t i = 1; i < l - k + 1; ++i) {
        hf.update(s[i+k-1], s[i]);
        if (hf.hashvalue % p == 0) {
            f(start, i+k-start); // user can do what they please with a phrase boundary
            start = i;
        }
    }
    // deal with the last phrase
    // TODO: add dollar signs to the last phrase
    return std::make_tuple(start, l-start);
}

template<typename Func>
std::tuple<size_t, size_t> select_windows_manzini_kr(const char* s, const size_t l, const size_t k, const size_t p, Func f) {
    KR_window krw(k);
    for (size_t i = 0; i < k; ++i) {
        krw.addchar(s[i]);
    } // we don't bother hashing the first kmer, since it's the start of a phrase anyway
    size_t start = 0;
    for (size_t i = 1; i < l - k + 1; ++i) {
        uint64_t hash = krw.addchar(s[i+k-1]);
        if (hash % p == 0) {
            f(start, i+k-start); // user can do what they please with a phrase boundary
            start = i;
        }
    }
    // deal with the last phrase
    // TODO: add dollar signs to the last phrase
    return std::make_tuple(start, l-start);
}

struct Freq {
    Freq(uint32_t x) : n(x) {}
    uint32_t n = 0;
    uint32_t r = 0;
};

using FreqMap = std::map<std::string, Freq>;
using OrderMap = std::vector<const char*>;

// TODO: multithreading
// we can expect most of the phrases that come out to have freq==1 for a typical genome
FreqMap parse_fasta(const Args& args, FreqMap& freqs, OrderMap& phrase_order) {
    gzFile fp = gzopen(args.in_fname.data(), "r");
    kseq_t* seq = kseq_init(fp);
    int l;
    // open relevant files
    FILE* last_fp = open_aux_file(args.in_fname.data(), EXTLST, "wb");
    FILE* sa_fp = NULL;
    if (args.sai)
        sa_fp = open_aux_file(args.in_fname.data(), EXTSAI, "wb");
    FILE* old_parse_fp = open_aux_file(args.in_fname.data(), "parse2", "w");
    std::string phrase;
    phrase.append(1, Dollar);
    size_t pos = 0;
    while (( l = kseq_read(seq) ) >= 0) {
        auto last_seq_idx = select_windows_2bit_wang(seq->seq.s, seq->seq.l, args.w, args.p,
                // process the phrases that select_windows spits out
                [&] (const size_t start, const size_t len) {
                    phrase.append(seq->seq.s + start, len);
                    auto ret = freqs.insert({phrase, Freq(1)});
                    if (!ret.second) ret.first->second.n += 1;
                    phrase_order.push_back(ret.first->first.data());
                    fprintf(old_parse_fp, "%s\n", ret.first->first.data());
                    pos += len-args.w;
                    if (fputc(phrase[phrase.size()-args.w-1], last_fp) == EOF) die("error writing to .last file\n");
                    if (args.sai) if(fwrite(&pos, IBYTES, 1, sa_fp) != 1) die("error writing to .sai file\n");
                    // next iter will return at start of phrase, so it's safe to clear
                    phrase.clear();
                });
        // prepend rest of sequence to next phrase
        phrase.append(seq->seq.s + std::get<0>(last_seq_idx), std::get<1>(last_seq_idx));
    }
    // deal with very last phrase
    phrase.append(args.w, Dollar);
    auto ret = freqs.insert({phrase, Freq(1)});
    if (!ret.second) ret.first->second.n += 1;
    phrase_order.push_back(ret.first->first.data());
    fprintf(old_parse_fp, "%s\n", ret.first->first.data());
    pos += phrase.size()-args.w;
    if (fputc(phrase[phrase.size()-args.w-1], last_fp) == EOF) die("error writing to .last file");
    if (args.sai) if (fwrite(&pos, IBYTES, 1, sa_fp) != 1) die("error writing to sai file");

    kseq_destroy(seq);
    gzclose(fp);
    if (fclose(last_fp)) die("error closing LAST file\n");
    if (args.sai) 
        if (fclose(sa_fp)) die("error closing SAI file\n");
    return freqs;
}

std::vector<const char*> write_dict_occ( const Args& args, FreqMap& freqs) {
    fprintf(stderr, "creating dict array\n");
    std::vector<const char*> dict_phrases;
    dict_phrases.reserve(freqs.size());
    for (auto it = freqs.begin(); it != freqs.end(); ++it) {
        dict_phrases.push_back(it->first.data()); // for some reason, dict_phrases.push_back does not work
    }
    fprintf(stderr, "sorting dict\n");
    std::sort(dict_phrases.begin(), dict_phrases.end(), 
            [](const char* l, const char* r) { return strcmp(l, r) <= 0; });

    fprintf(stderr, "writing dict and occ file...\n");
    FILE* dict_fp = open_aux_file(args.in_fname.data(), EXTDICT, "wb");
    FILE* occ_fp  = open_aux_file(args.in_fname.data(), EXTOCC, "wb");
    size_t rank = 1;
    for (auto x: dict_phrases) {
        if (fwrite(x, 1, strlen(x), dict_fp) != strlen(x)) 
            die("Error writing to DICT file\n");
        if (fputc(EndOfWord, dict_fp) == EOF) 
            die("Error writing EndOfWord to DICT file");
        // access dictionary and write occ to occ file
        auto& wf = freqs.at(x);
        if (fwrite(&wf.n, sizeof(wf.n), 1, occ_fp) != 1) 
            die("Error writing to OCC file\n");
        wf.r = rank++;
    }
    if (fputc(EndOfDict, dict_fp) == EOF) die("Error writing EndOfDict to DICT file");
    if (fclose(occ_fp)) die("Error closing OCC file");
    if (fclose(dict_fp)) die("Error closing DICT file");
    return dict_phrases;
}

void write_parse(const Args& args, FreqMap& freqs, const OrderMap& phrase_order) {
    FILE* parse_fp = open_aux_file(args.in_fname.data(), EXTPARSE, "wb");
    for (auto phrase: phrase_order) {
        auto wf = freqs.at(std::string(phrase));
        if (fwrite(&wf.r, sizeof(wf.r), 1, parse_fp) != 1) 
            die("Error writing to new parse file");
    }
    if (fclose(parse_fp)) die("Error closing new parse file");
    return;
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
                print_help(argv, args); exit(1);
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
    FreqMap freqs;
    OrderMap phrase_order;
    fprintf(stderr, "parsing...\n");
    auto freq = parse_fasta(args, freqs, phrase_order);
    fprintf(stderr, "writing dict...\n");
    {
    auto dict_phrases = write_dict_occ(args, freqs);
    }
    fprintf(stderr, "writing parse...\n");
    write_parse(args, freqs, phrase_order);
}
