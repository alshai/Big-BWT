#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <string>
#include <vector>
#include <fcntl.h>
#include "sdsl/bit_vectors.hpp"
extern "C" {
#include <sys/mman.h>
#include "utils.h"
#include "gsa/gsacak.h"
}

static constexpr int OCC_BYTES = sizeof(uint32_t);
static constexpr int ILIST_BYTES = sizeof(uint32_t);

struct SuffixT {
    SuffixT(uint8_t c, uint32_t i) : bwtc(c), bwtp(i) {}
    uint8_t bwtc;
    size_t bwtp;
    bool operator<(SuffixT& r);
};

bool SuffixT::operator<(SuffixT& r) {
    return this->bwtp < r.bwtp;
}

class pfbwt {

    public:

    pfbwt(std::string prefix, size_t win_size) : w(win_size) {
        load_files(prefix);
    }

    ~pfbwt() {
        free(dict);
        // if (munmap(dict, dsize) != 0) {
        //     die("error deallocating dict mmap");
        // }
    }


#define get_word_suflen(i, d, s) \
    d = dict_idx_rank(i); \
    s = dict_idx_select(d+1) - i;

     /* uses LCP of dict to build BWT (less memory, more time)
     */
    void generate_bwt_lcp(std::string prefix) {
        FILE* bwt_fp = open_aux_file(prefix.data(), "bwt", "wb");
        sort_dict_suffixes(true); // build SA and lcp
        // start from SA item that's not EndOfWord or EndOfDict
        size_t suff_len, wordi, ilist_pos;
        uint8_t bwtc, pbwtc = Dollar;
        size_t easy_cases = 0, hard_cases = 0;
        for (size_t i = dwords+w+1; i<dsize; ++i) { 
            get_word_suflen(sa[i], wordi, suff_len);
            if (suff_len <= w) continue; // ignore small suffixes
            // full word case
            if (sa[i] == 0 || dict_idx[sa[i]-1] == 1) {
                auto word_ilist = get_word_ilist(wordi);
                for (auto j: word_ilist) {
                    bwtc = bwlast[j];
                    fprintf(stderr, "%llu\n", j);
                    // TODO: do SA/DA/etc stuff here,
                    // TODO: uncomment this next line
                    // if (fputc(bwtc, bwt_fp) == EOF) die("error writing to bwt file");
                    pbwtc = bwtc; // for SA
                }
                ++easy_cases;
            } else { // hard case!
                // look at all the sufs that share LCP[suf]==this_suffixlen
                size_t nwordi, nsuff_len;
                std::vector<uint8_t> chars;
                std::vector<uint64_t> words;
                uint8_t c, pc = dict[sa[i]-1];
                chars.push_back(c);
                words.push_back(wordi);
                bool same_char = true;
                for (size_t j = i+1; j < dsize && lcp[j] == (int_t) suff_len; ++j) {
                    // get_word_suflen(sa[j], nwordi, nsuff_len);
                    if (nsuff_len != suff_len) die("something went wrong!\n");
                    c = dict[sa[j]-1];
                    // chars.push_back(c);
                    // words.push_back(nwordi);
                    same_char = c == pc;
                    pc = c;
                }
                if (same_char) {
                    // fprintf(stderr, "non-full-word, easy case\n");
                    ++easy_cases;
                } else {
                    /*
                    std::vector<SuffixT> suffs;
                    for (int j = 0; j < words.size(); ++j) {
                        // get ilist of each of these words, make a heap
                        for (auto k: get_word_ilist(words[j])) {
                            suffs.push_back(SuffixT(chars[j], k));
                        }
                    }
                    // std::sort(suffs.begin(), suffs.end());
                    // fprintf(stderr, "-----\n");
                    // for (auto s: suffs) {
                    //     fprintf(stderr, "%llu %c\n", s.bwtp, s.bwtc);
                    // }
                    */
                    ++hard_cases;
                }
                chars.clear();
                words.clear();
            }
        }
        return;
    }

    void generate_bwt_fm() {
        return;
    }

    private:

    /* run gSACAK on d 
     * populates sa, lcp, and dict_idx;
     */
    void sort_dict_suffixes(bool build_lcp = true) {
        if (dsize < 1) die("error: dictionary not loaded\n");
        fprintf(stderr, "running gsacak\n");
        sa.resize(dsize);
        lcp.resize(dsize);
        if (build_lcp) 
            gsacak(dict, &sa[0], &lcp[0], NULL, dsize);
        else { // for when memory is low
            gsacak(dict, &sa[0], NULL, NULL, dsize);
            // TODO: build FM index over dict
        }
        // make index of dict end positions
        dict_idx = sdsl::bit_vector(dsize, 0);
        for (size_t i = 1; i < dwords + 1; ++i) {
            dict_idx[sa[i]] = 1;
        }
        sdsl::util::init_support(dict_idx_rank, &dict_idx);
        sdsl::util::init_support(dict_idx_select, &dict_idx);
    }

    void load_files(std::string prefix) {
        load_dict(prefix);
        load_ilist_idx(prefix);
        load_ilist(prefix);
        load_bwlast(prefix);
    }

    void mmap_dict(std::string fname) {
        int dict_fd = fd_open_aux_file(fname.data(), EXTDICT, O_RDONLY);
        if (dict_fd < 0) die("error opening dict_fd");
        dict = (uint8_t*) mmap(NULL, dsize, PROT_READ, MAP_PRIVATE, dict_fd, 0);
    }

    void malloc_dict(FILE* dict_fp) {
        dict = (uint8_t*) malloc(sizeof(uint8_t) * dsize);
        if (fread(dict, sizeof(uint8_t), dsize, dict_fp) != (size_t) dsize)
            die("error reading dict_fp");
    }

    void load_dict(std::string fname) {
        FILE* dict_fp = open_aux_file(fname.data(), EXTDICT, "rb");
        // get dict size
        fseek(dict_fp, 0, SEEK_END);
        dsize = ftell(dict_fp);
        if (dsize <= 1 + w) die("dict is not long enough (< w)");
        rewind(dict_fp);
        malloc_dict(dict_fp); // TODO: check performance of mmap_dict for really large files
        fclose(dict_fp);
        int x = 0;
        for (size_t i = 0; i < dsize; ++i) {
            if (dict[i] == EndOfWord) {
                x += 1;
            }
        }
        fprintf(stderr, "dict[0]: %c\n", dict[0]);
        fprintf(stderr, "words found in dict: %d\n", x);
    }

    void load_ilist_idx(std::string fname) {
        FILE* occ_fp = open_aux_file(fname.data(), EXTOCC, "rb");
        // get occs size
        fseek(occ_fp, 0, SEEK_END);
        size_t osize = ftell(occ_fp);
        if (osize % 4 != 0) die("invalid occ file");
        dwords = osize / OCC_BYTES;
        rewind(occ_fp);
        std::vector<uint32_t> occs(dwords);
        if (fread(&occs[0], OCC_BYTES, occs.size(), occ_fp) != dwords) {
            die("error reading occs");
        }
        int total_occs = 0;
        for (auto o: occs) total_occs += o;
        ilist_idx = sdsl::bit_vector(total_occs+occs[dwords], 0);
        size_t o = 0;
        for (size_t i = 0; i < occs.size(); ++i) {
            o += occs[i];
            ilist_idx[o-1] = 1;
        }
        sdsl::util::init_support(ilist_idx_rank,   &ilist_idx);
        sdsl::util::init_support(ilist_idx_select, &ilist_idx);
        fclose(occ_fp);
    }

    void load_ilist(std::string fname) {
        ilist.clear();
        FILE* ilist_fp = open_aux_file(fname.data(), EXTILIST, "rb");
        // get ilist size
        fseek(ilist_fp, 0, SEEK_END);
        size_t ilsize = ftell(ilist_fp);
        if (ilsize % 4 != 0) die("invalid ilist file");
        size_t nelems = ilsize / ILIST_BYTES;
        rewind(ilist_fp);
        ilist.resize(nelems);
        if (fread(&ilist[0], ILIST_BYTES, ilist.size(), ilist_fp) != nelems) {
            die("error reading ilist file");
        }
        fclose(ilist_fp);
        fprintf(stderr, "ilist contains %lu elements\n", nelems);
    }

    void load_bwlast(std::string fname) {
        bwlast.clear();
        FILE* bwlast_fp = open_aux_file(fname.data(), EXTBWLST, "rb");
        // get bwlast size
        fseek(bwlast_fp, 0, SEEK_END);
        size_t bwlsize = ftell(bwlast_fp);
        // load occs as 1s in bit vector
        rewind(bwlast_fp);
        bwlast.resize(bwlsize); // TODO
        if (fread(&bwlast[0], sizeof(uint8_t), bwlast.size(), bwlast_fp) != bwlsize) {
            die("error reading bwlast file");
        }
        fclose(bwlast_fp);
    }

    std::vector<size_t> get_word_ilist(size_t wordi) {
        // get to the end of the previous word's list, then add one to get
        // to the start of the current word
        auto startpos = wordi ? ilist_idx_select(wordi) + 1 : 0;
        auto endpos = ilist_idx_select(wordi+1);
        std::vector<size_t> v;
        v.reserve(endpos - startpos + 1);
        for (size_t j = startpos; j < endpos+1; ++j)
            v.push_back(ilist[j]);
        return v;
    }

    size_t w; // word size of parser
    // note: can be malloc'd (malloc_dict) or mmap'd (mmap_dict)
    uint8_t* dict; // dict word array (word ends represented by EndOfWord)
    bool mmapped = false;
    uint64_t dsize; // number of characters in dict
    uint64_t dwords; // number of words in dict
    std::vector<uint8_t> bwlast; // parse-bwt char associated w/ ilist
    std::vector<uint32_t> ilist; // bwlast positions of dict words
    sdsl::bit_vector ilist_idx; // 1 on ends of dict word occs in ilist
    sdsl::bit_vector dict_idx; // 1 on word end positions in dict
    sdsl::bit_vector::rank_1_type   ilist_idx_rank;
    sdsl::bit_vector::rank_1_type   dict_idx_rank;
    sdsl::bit_vector::select_1_type ilist_idx_select;
    sdsl::bit_vector::select_1_type dict_idx_select;
    std::vector<uint_t> sa; // gSA of dict words
    std::vector<int_t> lcp; // gLCP of dict words
};
