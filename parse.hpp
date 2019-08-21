#ifndef PARSE_HPP
#define PARSE_HPP

#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <zlib.h>
#include <kseq.h>
KSEQ_INIT(gzFile, gzread);
extern "C" {
#include "utils.h"
}
#include "hash.h"

struct Args {
    std::string in_fname;
    size_t w = 10;
    size_t p = 100;
    bool sai = false;
};

extern uint8_t seq_nt4_ntoa_table[];

struct Freq {
    Freq(uint32_t x) : n(x) {}
    uint32_t n = 0;
    uint32_t r = 0;
};

struct WangHash {
    WangHash(size_t w) :
        k(w),
        mask((1ULL << 2 * k) - 1)
    {}

    uint64_t update(char c) {
        char x = seq_nt4_ntoa_table[(size_t) c];
        if (x > 3) { fprintf(stderr, "error, invalid character %c\n", x); exit(1);}
        kmer = ((kmer << 2) | x) & mask;
        hash = wang_hash(kmer);
        return hash;
    }

    uint64_t hashvalue() { return hash; }

    size_t k;
    uint64_t kmer = 0;
    uint64_t mask;
    uint64_t hash;
};

/* TODO finish this
struct KRHash {
    KRHash(size_t wsize) : 
        window(wsize,0), 
        k(wsize) { }

    uint64_t update(char c) {
        window[i] = c;
        i = (i+1) % k;
        return 1;
    }

    uint64_t hashvalue() { return hash; }

    std::vector<char> window;
    int k = 10;
    int i = 0;
    uint64_t hash = 0;
};
*/

using FreqMap=std::map<std::string, Freq>;

template <typename Hasher>
class Parser {

    public: 

    Parser(int wsize, int pmod) :
        w(wsize),
        p(pmod) { }

    void parse_fasta(const char* fname, const bool sai = false) {
        FILE* last_fp = open_aux_file(fname, EXTLST, "wb");
        FILE* sa_fp = NULL;
        if (sai) sa_fp = open_aux_file(fname, EXTSAI, "wb");
        // FILE* old_parse_fp = open_aux_file(fname, "parse2", "w");
        gzFile fp = gzopen(fname, "r");
        kseq_t* seq = kseq_init(fp);
        int l;
        size_t nseqs(0);
        uint64_t pos = 0;
        std::string phrase;
        phrase.append(1, Dollar);
        Hasher hf(w);
        while (( l = kseq_read(seq) ) >= 0) {
            for (size_t i = 0; i < seq->seq.l; ++i) {
                phrase.append(1, seq->seq.s[i]);
                hf.update(seq->seq.s[i]);
                if (hf.hashvalue() % p == 0) {
                    // fprintf(old_parse_fp, "%s\n", phrase.data());
                    pos = pos ? pos+phrase.size()-w : phrase.size()-1;
                    process_phrase(phrase, last_fp, sa_fp, &pos);
                    phrase.erase(0, phrase.size()-w);
                }
            }
            ++nseqs;
        }
        phrase.append(w, Dollar);
        // fprintf(old_parse_fp, "%s\n", phrase.data());
        pos = pos ? pos+phrase.size()-w : phrase.size()-1;
        process_phrase(phrase, last_fp, sa_fp, &pos);

        kseq_destroy(seq);
        gzclose(fp);
        // if (fclose(old_parse_fp)) die("error closing PARSE2 file\n");
        // else fprintf(stderr, "PARSE2 file written to %s.parse2\n", fname);
        if (fclose(last_fp)) die("error closing LAST file\n");
        else fprintf(stderr, "LAST file written to %s.%s\n", fname, EXTLST);
        if (sai) {
            if (fclose(sa_fp)) die("error closing SAI file\n");
            else (fprintf(stderr, "SAI file written to %s.%s\n", fname, EXTSAI));
        }
        return;
    }

    // assigns lexicographic rankings to items in dictionary
    // if fname is provided, dumps dictionary and occs to file
    void update_dict(const char* fname = NULL) {
        std::vector<const char*> dict_phrases;
        dict_phrases.reserve(freqs.size());
        for (auto it = freqs.begin(); it != freqs.end(); ++it) {
            dict_phrases.push_back(it->first.data());
        }
        fprintf(stderr, "sorting dict\n");
        std::sort(dict_phrases.begin(), dict_phrases.end(),
                [](const char* l, const char* r) { return strcmp(l, r) <= 0; });
        fprintf(stderr, "writing dict and occ file...\n");
        FILE* dict_fp = NULL;
        FILE* occ_fp = NULL;
        if (fname != NULL) {
            dict_fp = open_aux_file(fname, EXTDICT, "wb");
            occ_fp  = open_aux_file(fname, EXTOCC, "wb");
        }
        size_t rank = 1;
        for (auto x: dict_phrases) {
            // access dictionary and write occ to occ file
            auto& wf = freqs.at(x);
            wf.r = rank++;
            if (fname != NULL) {
                if (fwrite(x, 1, strlen(x), dict_fp) != strlen(x))
                    die("Error writing to DICT file\n");
                if (fputc(EndOfWord, dict_fp) == EOF)
                    die("Error writing EndOfWord to DICT file");
                if (fwrite(&wf.n, sizeof(wf.n), 1, occ_fp) != 1)
                    die("Error writing to OCC file\n");
            }
        }
        if (fname != NULL) {
            if (fputc(EndOfDict, dict_fp) == EOF) die("Error writing EndOfDict to DICT file");
            if (fclose(occ_fp)) die("Error closing OCC file");
            else fprintf(stderr, "OCC written to %s.%s\n", fname, EXTOCC);
            if (fclose(dict_fp)) die("Error closing DICT file");
            else fprintf(stderr, "DICT written to %s.%s\n", fname, EXTDICT);
        }
    }

    void dump_parse(const char* fname) {
        FILE* parse_fp = open_aux_file(fname, EXTPARSE, "wb");
        for (auto phrase: phrase_order) {
            auto wf = freqs.at(std::string(phrase));
            if (fwrite(&wf.r, sizeof(wf.r), 1, parse_fp) != 1)
                die("Error writing to new parse file");
        }
        if (fclose(parse_fp)) die("Error closing PARSE file");
        else fprintf(stderr, "PARSE written to %s.%s\n", fname, EXTPARSE);
        return;
    }

    private:

    void inline process_phrase(const std::string& phrase, FILE* last_fp=NULL, FILE* sa_fp=NULL, uint64_t* pos=NULL) {
        auto ret = freqs.insert({phrase, Freq(1)});
        if (!ret.second) ret.first->second.n += 1;
        phrase_order.push_back(ret.first->first.data());
        if (last_fp != NULL){
            if (fputc(phrase[phrase.size()-w-1], last_fp) == EOF) 
                die("error writing to .last file\n");
        }
        if (sa_fp != NULL) {
            if(fwrite(pos, IBYTES, 1, sa_fp) != 1) 
                die("error writing to .sai file\n");
        }
    }

    FreqMap freqs;
    std::vector<const char*> phrase_order;
    int w, p;
};


#endif
