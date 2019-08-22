#ifndef HASH_HPP
#define HASH_HPP
#include <cinttypes>
#include <cstdio>
#include <cstdlib>

extern uint8_t seq_nt4_ntoa_table[];

// copied from daniel baker's bonsai repo
// https://github.com/dnbaker/bonsai/blob/437c06218217e14704bc8903b25c9934f377642c/bonsai/include/hash.h#L13
inline uint64_t wang_hash(uint64_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

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

#endif
