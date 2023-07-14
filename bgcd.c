#include <assert.h>
#include <stdint.h>

uint32_t count_trailing_zeros(uint32_t x)
{
    assert(x != 0);
    const uint32_t debruijn_ctz[] = {0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};
    return debruijn_ctz[((x & -x) * 0x077cb531) >> 27];
}

uint32_t bgcd32(uint32_t a, uint32_t b)
{
    if (a == 0 || b == 0)
    {
        return a | b;
    }
    uint32_t s, t;
    s = count_trailing_zeros(a | b);
    a >>= count_trailing_zeros(a);
    do {
        b >>= count_trailing_zeros(b);
        if (a > b)
        {
            t = a;
            a = b;
            b = t;
        }
        b -= a;
    } while (b);
    return a << s;
}

uint64_t count_trailing_zeros(uint64_t x)
{
    assert(x != 0);
    const uint64_t debruijn_ctz[] = {0, 1, 2, 53, 3, 7, 54, 27, 4, 38, 41, 8, 34, 55, 48, 28, 62, 5, 39, 46, 44, 42, 22, 9, 24, 35, 59, 56, 49, 18, 29, 11, 63, 52, 6, 26, 37, 40, 33, 47, 61, 45, 43, 21, 23, 58, 17, 10, 51, 25, 36, 32, 60, 20, 57, 16, 50, 31, 19, 15, 30, 14, 13, 12};
    return debruijn_ctz[((x & -x) * 0x022fdd63cc95386d) >> 58];
}

uint64_t bgcd64(uint64_t a, uint64_t b)
{
    if (a == 0 || b == 0)
    {
        return a | b;
    }
    uint64_t s, t;
    s = count_trailing_zeros(a | b);
    a >>= count_trailing_zeros(a);
    do {
        b >>= count_trailing_zeros(b);
        if (a > b)
        {
            t = a;
            a = b;
            b = t;
        }
        b -= a;
    } while (b);
    return a << s;
}
