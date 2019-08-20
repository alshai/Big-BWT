#include <cinttypes>

// copied from daniel baker's bonsai repo
// https://github.com/dnbaker/bonsai/blob/437c06218217e14704bc8903b25c9934f377642c/bonsai/include/hash.h#L13
constexpr inline uint64_t wang_hash(uint64_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

struct KR_window {
    int wsize;
    int *window;
    int asize;
    const uint64_t prime = 1999999973;
    uint64_t hash;
    uint64_t tot_char;
    uint64_t asize_pot;   // asize^(wsize-1) mod prime

    KR_window(int w): wsize(w) {
        asize = 256;
        asize_pot = 1;
        for(int i=1;i<wsize;i++)
            asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
        // alloc and clear window
        window = new int[wsize];
        reset();
    }

    // init window, hash, and tot_char
    void reset() {
        for(int i=0;i<wsize;i++) window[i]=0;
        // init hash value and related values
        hash=tot_char=0;
    }

    uint64_t addchar(int c) {
        int k = tot_char++ % wsize;
        // complex expression to avoid negative numbers
        hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
        hash = (asize*hash + c) % prime;      //  add char i
        window[k]=c;
        // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
        return hash;
    }
    // debug only
    std::string get_window() {
        std::string w = "";
        int k = (tot_char-1) % wsize;
        for(int i=k+1;i<k+1+wsize;i++)
            w.append(1,window[i%wsize]);
        return w;
    }

    ~KR_window() {
        delete[] window;
    }

};
