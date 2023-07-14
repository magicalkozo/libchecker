#include <math.h>
#include <stdint.h>
#include <limits.h>

typedef uint64_t u64;

u64 pow_with_upper_bound(u64 x, u64 y)
{
    if (y == 0)
    {
        return 1ull;
    }
    u64 res = 1ull;
    while (y)
    {
        if (y & 1)
        {
            res = __builtin_mul_overflow_p(res, x, (u64)0) ? UINT64_MAX : res * x;
        }
        x = __builtin_mul_overflow_p(x, x, (u64)0) ? UINT64_MAX : x * x;
        y >>= 1;
    }
    return res;
}
u64 floor_kth_root_integer(u64 a, u64 k)
{
    if (a <= 1 || k == 1)
    {
        return a;
    }
    if (k >= 64)
    {
        return 1;
    }
    if (k == 2)
    {
        return sqrtl(a);
    }
    if (a == UINT64_MAX)
    {
        a--;
    }
    u64 res = (k == 3 ? cbrt(a) - 1 : pow(a, nextafter(1 / (double)k, 0)));
    while (pow_with_upper_bound(res + 1, k) <= a)
    {
        res++;
    }
    return res;
}

// binary search version
u64 kth_root_integer(u64 a, int k)
{
    if (k == 1)
    {
        return a;
    }
    u64 low = 0;
    u64 up = 1099511627776ull;
    while (up - low > 1)
    {
        u64 mid = low + ((up - low) >> 1);
        u128 mypow = 1ull;
        for (int i = 0; i < k; i++)
        {
            ret *= mid;
            if (ret >> 64)
            {
                break;
            }
        }
        if (mypow <= a)
        {
            low = mid;
        }
        else
        {
            up = mid;
        }
    }
    return low;
}
