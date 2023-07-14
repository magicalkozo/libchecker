#include <assert.h>
#include <stdint.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef __uint128_t u128;

static u64 m_b32, im_b32, div_b32, rem_b32;
void set_b32(u64 mod)
{
    m_b32 = mod;
    im_b32 = ((((u128)(1ull) << 64)) + mod - 1) / mod;
    div_b32 = 0;
    rem_b32 = 0;
}
void br32(u64 a)
{
    u64 x = (u64)(((u128)(a) * im_b32) >> 64);
    u64 y = x * m_b32;
    unsigned long long z;
    u32 w = __builtin_usubll_overflow(a, y, &z) ? m_b32 : 0;
    div_b32 = x;
    rem_b32 = z + w;
}
u32 add_b32(u32 a, u32 b)
{
    a += b;
    a -= (a >= (u32)m_b32 ? (u32)m_b32 : 0);
    return a;
}
u32 sub_b32(u32 a, u32 b)
{
    a += (a < b ? (u32)m_b32 : 0);
    a -= b;
    return a;
}
u32 min_b32(u32 a)
{
    return sub_b32(0, a);
}
u32 mul_b32(u32 a, u32 b)
{
    br32((u64)a * b);
    return (u32)rem_b32;
}
u32 squ_b32(u32 a)
{
    br32((u64)a * a);
    return (u32)rem_b32;
}
u32 pow_b32(u32 a, u32 k)
{
    u32 ret = a, deg = k;
    while (deg > 0)
    {
        if (deg & 1)
        {
            ret = mul_b32(ret, a);
        }
        a = squ_b32(a);
        deg >>= 1;
    }
    return ret;
}
u32 shl_b32(u32 a)
{
    return (a <<= 1) >= m_b32 ? a - m_b32 : a;
}
u32 shr_b32(u32 a)
{
    return (a & 1) ? ((a >> 1) + (m_b32 >> 1) + 1) : (a >> 1);
}