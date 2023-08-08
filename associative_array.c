#include "template.h"

const u64 hm_s = 20ull;
const u64 hm_t = 64ull - hm_s;
const u64 hm_m = 1514165138286539305ull;

typedef Pair_u64 Hashmap;

static Hashmap hm[1048577];
static bool hm_used[1048577];

#define hm_hash(x) (((u64)(x)) * hm_m >> hm_t)

__attribute__((constructor)) void _construct_hashmap_(void) {
    for (size_t i = 0; i < 1048577; i++) {
        hm[i] = (Hashmap){0, 0};
        hm_used[i] = false;
    }
}

void hm_insert(u64 k_, u64 v_) {
    int i = hm_hash(k_);
    if (!hm_used[i]) {
        hm[i].a = k_;
        hm[i].b = v_;
        hm_used[i] = true;
        return;
    }
    for (; hm_used[i]; i++) {
        if (hm[i].a == k_) {
            hm[i].b = v_;
            return;
        }
    }
    hm[i].a = k_;
    hm[i].b = v_;
    hm_used[i] = true;
}
u64 hm_find(u64 k_) {
    for (int i = hm_hash(k_); hm_used[i]; i++) {
        if (hm[i].a == k_) {
            return hm[i].b;
        }
    }
    return 0;
}

#undef hm_hash

void solve_associative_array(void) {
    unsigned Q;
    rd_uint(&Q);
    unsigned query;
    unsigned long long key, val;
    while (Q--) {
        rd_uint(&query);
        rd_ull(&key);
        if (query) {
            wt_ull(hm_find(key));
            wt_nl();
        } else {
            rd_ull(&val);
            hm_insert(key, val);
        }
    }
}
