#include <stdint.h>

uint64_t gcd64(uint64_t a, uint64_t b)
{
    if (a == 0 || b == 0)
    {
        return a | b;
    }
    uint64_t s, t;
    s = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b)
        {
            t = a, a = b, b = t;
        }
        b -= a;
    } while (b);
    return a << s;
}
uint64_t powmod(uint64_t a, uint64_t k, uint64_t n)
{
    uint64_t ret = 1;
    while (k > 0)
    {
        if (k & 1)
        {
            ret = ((__uint128_t)ret * a) % n;
        }
        a = ((__uint128_t)a * a) % n;
        k >>= 1;
    }
    return ret;
}

int miller_rabin(uint64_t n)
{
    uint64_t s = __builtin_ctzll(n - 1);
    uint64_t d = (n - 1) >> s;
    uint64_t bases[7] = {2ull, 325ull, 9375ull, 28178ull, 450775ull, 9780504ull, 1795265022ull};
    for (int i = 0; i < 7; i++)
    {
        if (n <= bases[i])
        {
            return 1;
        }
        uint64_t a = pow_mod(bases[i], d, n);
        if (a == 1) { continue; }
        uint64_t r = 1;
        while (a != n - 1)
        {
            if (r == s) { return 0; }
            a = ((__uint128_t)a * a) % n;
            r++;
        }
    }
    return 1;
}

int is_prime(uint64_t n)
{
    if (n < 64ull) { return (1ull << n) & 2891462833508853932ull; }
    if (!(n & 1)) { return 0; }
    if (gcd64(n, 16294579238595022365ull) != 1ull) { return 0; }
    return miller_rabin(n);
}