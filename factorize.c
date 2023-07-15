#pragma region template
// clang-format off

#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("unroll-loops")

#define _GNU_SOURCE
#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

typedef   int8_t      i8;
typedef   int16_t     i16;
typedef   int32_t     i32;
typedef   int64_t     i64;
typedef __int128_t    i128;
typedef   uint8_t     u8;
typedef   uint16_t    u16;
typedef   uint32_t    u32;
typedef   uint64_t    u64;
typedef __uint128_t   u128;
typedef   float       f32;
typedef   double      f64;

// https://gcc.gnu.org/onlinedocs/gcc/C-Extensions.html
// https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
#define clz32(a) ((a) ? __builtin_clz((a)) : (32))
#define clz64(a) ((a) ? __builtin_clzll((a)) : (64))
#define ctz32(a) ((a) ? __builtin_ctz((a)) : (32))
#define ctz64(a) ((a) ? __builtin_ctzll((a)) : (64))
#define pct32(a) __builtin_popcount((a))
#define pct64(a) __builtin_popcountll((a))
#define msb32(a) ((a) ? ((31) - __builtin_clz((a))) : (32))
#define msb64(a) ((a) ? ((63) - __builtin_clzll((a))) : (64))
#define bit_width32(a) ((a) ? ((32) - __builtin_clz((a))) : (0))
#define bit_width64(a) ((a) ? ((64) - __builtin_clzll((a))) : (0))
#define bit_ceil32(a) ((!(a)) ? (1) : ((pct32(a)) == (1) ? ((1u) << ((31) - clz32((a)))) : ((1u) << ((32) - clz32(a)))))
#define bit_ceil64(a) ((!(a)) ? (1) : ((pct64(a)) == (1) ? ((1ull) << ((63) - clz64((a)))) : ((1ull) << ((64) - clz64(a)))))
#define bit_floor32(a) ((!(a)) ? (0) : ((1u) << ((31) - clz32((a)))))
#define bit_floor64(a) ((!(a)) ? (0) : ((1ull) << ((63) - clz64((a)))))
#define flip_Nbit(a, n) ((a) ^ ((1) << (n)))
#define only_lsb(a) ((a) & (-(a)))
u32 bit_reverse32(u32 v) { v = (v & 0x55555555) << 1 | (v >> 1 & 0x55555555); v = (v & 0x33333333) << 2 | (v >> 2 & 0x33333333); v = (v & 0x0f0f0f0f) << 4 | (v >> 4 & 0x0f0f0f0f); v = (v & 0x00ff00ff) << 8 | (v >> 8 & 0x00ff00ff); v = (v & 0x0000ffff) << 16 | (v >> 16 & 0x0000ffff); return v; }
u64 bit_reverse64(u64 v) { v = (v & 0x5555555555555555) <<  1 | (v >>  1 & 0x5555555555555555); v = (v & 0x3333333333333333) <<  2 | (v >>  2 & 0x3333333333333333); v = (v & 0x0f0f0f0f0f0f0f0f) <<  4 | (v >>  4 & 0x0f0f0f0f0f0f0f0f); v = (v & 0x00ff00ff00ff00ff) <<  8 | (v >>  8 & 0x00ff00ff00ff00ff); v = (v & 0x0000ffff0000ffff) << 16 | (v >> 16 & 0x0000ffff0000ffff); v = (v & 0x00000000ffffffff) << 32 | (v >> 32 & 0x00000000ffffffff); return v; }
u32 rotr32(u32 x, i64 r) { if (r < 0) { u32 l = (u32)(((u64)(-r)) % 32); return x << l | x >> (-l & 31); } r %= 32; return x >> r | x << (-r & 31); }
u64 rotr64(u64 x, i64 r) { if (r < 0) { u64 l = ((u64)(-r)) % 64; return x << l | x >> (-l & 63); } r %= 64; return x >> r | x << (-r & 63); }
u32 rotl32(u32 x, i64 r) { return rotr32(x, -r); }
u64 rotl64(u64 x, i64 r) { return rotr64(x, -r); }

u32 in_u32(void) { u32 c, x = 0; while (c = getchar_unlocked(), c < '0' || c > '9') {} while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return x; }
u64 in_u64(void) { u64 c, x = 0; while (c = getchar_unlocked(), c < '0' || c > '9') {} while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return x; }
u128 in_u128(void) { u128 c, x = 0; while (c = getchar_unlocked(), c < '0' || c > '9') {} while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return x; }
i32 in_i32(void) { i32 c, x = 0, f = 1; while (c = getchar_unlocked(), c < '0' || c > '9') { if (c == '-') { f = -f; } } while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return f * x; }
i64 in_i64(void) { i64 c, x = 0, f = 1; while (c = getchar_unlocked(), c < '0' || c > '9') { if (c == '-') { f = -f; } } while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return f * x; }
i128 in_i128(void) { i128 c, x = 0, f = 1; while (c = getchar_unlocked(), c < '0' || c > '9') { if (c == '-') { f = -f; } } while ('/' < c && c < ':') { x = x * 10 + c - '0'; c = getchar_unlocked(); } return f * x; }
void out_u32(u32 x) { if (x >= 10) { out_u32(x / 10); } putchar_unlocked(x - x / 10 * 10 + '0'); }
void out_u64(u64 x) { if (x >= 10) { out_u64(x / 10); } putchar_unlocked(x - x / 10 * 10 + '0'); }
void out_u128(u128 x) { if (x >= 10) { out_u128(x / 10); } putchar_unlocked(x - x / 10 * 10 + '0'); }
void out_i32(i32 x) { if (x < 0) { putchar_unlocked('-'); x = -x; } out_u32((u32)x); }
void out_i64(i64 x) { if (x < 0) { putchar_unlocked('-'); x = -x; } out_u64((u64)x); }
void out_i128(i128 x) { if (x < 0) { putchar_unlocked('-'); x = -x; } out_u128((u128)x); }
void newline(void) { putchar_unlocked('\n'); }
void space(void) { putchar_unlocked(' '); }
void printb_32bit(u32 v) { u32 mask = (u32)1 << (sizeof(v) * CHAR_BIT - 1); do { putchar_unlocked(mask & v ? '1' : '0'); } while (mask >>= 1); }
void printb_64bit(u64 v) { u64 mask = (u64)1 << (sizeof(v) * CHAR_BIT - 1); do { putchar_unlocked(mask & v ? '1' : '0'); } while (mask >>= 1); }

// https://en.wikipedia.org/wiki/Linear_congruential_generator
u32 lcg_rand(void) { static u64 lcg_state = 14534622846793005ull; lcg_state = 6364136223846793005ull * lcg_state + 1442695040888963407ull; return (u32)lcg_state; }
u32 lcg_range(u32 l, u32 r) { return l + lcg_rand() % (r - l + 1); }
f32 lcg_randf(void) { u32 a = 0x3F800000u | (lcg_rand() >> 9); return (*((f32 *)(&a))) - 1; }
// https://en.wikipedia.org/wiki/Permuted_congruential_generator
// https://www.pcg-random.org/index.html
u32 pcg_rand(void) { static u64 pcg_state = 0x853c49e6748fea9bull; u64 t = pcg_state; pcg_state = t * 0x5851f42d4c957f2dull + 0xda3e39cb94b95bdbull; u32 sh = ((t >> 18u) ^ t) >> 27u; u32 ro = t >> 59u; return (sh >> ro) | (sh << ((-ro) & 31)); }
u32 pcg_range(u32 l, u32 r) { return l + pcg_rand() % (r - l + 1); }
f32 pcg_randf(void) { u32 a = 0x3F800000u | (pcg_rand() >> 9); return (*((f32 *)(&a))) - 1; }
// https://en.wikipedia.org/wiki/Middle-square_method
// https://doi.org/10.48550/arXiv.1704.00358
u64 msws_rand(void) { static u64 msws_state1 = 0; static u64 msws_state2 = 0; static u64 msws_state3 = 0xb5ad4eceda1ce2a9ul; static u64 msws_state4 = 0; static u64 msws_state5 = 0; static u64 msws_state6 = 0x278c5a4d8419fe6bul; u64 ret; msws_state1 *= msws_state1; ret = msws_state1 += (msws_state2 += msws_state3); msws_state1 = (msws_state1 >> 32) | (msws_state1 << 32); msws_state4 *= msws_state4; msws_state4 += (msws_state5 += msws_state6); msws_state4 = (msws_state4 >> 32) | (msws_state4 << 32); return ret ^ msws_state4; }
u64 msws_range(u64 l, u64 r) { return l + msws_rand() % (r - l + 1); }
f64 msws_randf(void) { u64 a = 0x3FF0000000000000ull | (msws_rand() >> 12); return (*((f64 *)(&a))) - 1; }
// https://ja.wikipedia.org/wiki/Xorshift
// https://prng.di.unimi.it/xoroshiro128plus.c
u64 xrsr128p_rand(void) { static u64 xrsr128p_state1 = 0x1ull; static u64 xrsr128p_state2 = 0x2ull; const u64 s0 = xrsr128p_state1; u64 s1 = xrsr128p_state2; const u64 ret = s0 + s1; s1 ^= s0; xrsr128p_state1 = rotl64(s0, 24) ^ s1 ^ (s1 << 16); xrsr128p_state2 = rotl64(s1, 37); return ret; }
u64 xrsr128p_range(u64 l, u64 r) { return l + xrsr128p_rand() % (r - l + 1); }
f64 xrsr128p_randf(void) { u64 a = 0x3FF0000000000000ull | (xrsr128p_rand() >> 12); return (*((f64 *)(&a))) - 1; }
// https://prng.di.unimi.it/xoroshiro128starstar.c
u64 xrsr128ss_rand(void) { static u64 xrsr128ss_state1 = 0x1ull; static u64 xrsr128ss_state2 = 0x2ull; const u64 s0 = xrsr128ss_state1; u64 s1 = xrsr128ss_state2; const u64 ret = rotl64(s0 * 5, 7) * 9; s1 ^= s0; xrsr128ss_state1 = rotl64(s0, 24) ^ s1 ^ (s1 << 16); xrsr128ss_state2 = rotl64(s1, 37); return ret; }
u64 xrsr128ss_range(u64 l, u64 r) { return l + xrsr128ss_rand() % (r - l + 1); }
f64 xrsr128ss_randf(void) { u64 a = 0x3FF0000000000000ull | (xrsr128ss_rand() >> 12); return (*((f64 *)(&a))) - 1; }

u64 saturate_pow_u64(u64 x, u64 y) { if (y == 0) { return 1ull; } u64 res = 1ull; while (y) { if (y & 1) { res = __builtin_mul_overflow_p(res, x, (u64)0) ? UINT64_MAX : res * x; } x = __builtin_mul_overflow_p(x, x, (u64)0) ? UINT64_MAX : x * x; y >>= 1; } return res; }
u64 floor_kth_root_integer(u64 a, u64 k) { if (a <= 1 || k == 1) { return a; } if (k >= 64) { return 1; } if (k == 2) { return sqrtl(a); } if (a == UINT64_MAX) { a--; } u64 res = (k == 3 ? cbrt(a) - 1 : pow(a, nextafter(1 / (double)k, 0))); while (saturate_pow_u64(res + 1, k) <= a) { res++; } return res; }

// https://en.wikipedia.org/wiki/Jacobi_symbol
int jacobi_symbol(long long a, long long n) { int j = 1; long long t; while (a) { if (a < 0) { a = -a; if ((n & 3) == 3) { j = -j; } } int s = ctz64(a); a >>= s; if ((((n & 7) == 3) || ((n & 7) == 5)) && (s & 1)) { j = -j; } if ((a & n & 3) == 3) { j = -j; } t = a, a = n, n = t; a %= n; if ((a << 1) > n) { a -= n; } } return n == 1 ? j : 0; }

// https://en.wikipedia.org/wiki/Binary_GCD_algorithm#Algorithm
u32 gcd32(u32 a, u32 b) { if (!a || !b) { return a | b; } u32 t, s = ctz32(a | b); a >>= ctz32(a); do { b >>= ctz32(b); if (a > b) { t = a, a = b, b = t; } b -= a; } while (b); return a << s; }
u64 gcd64(u64 a, u64 b) { if (!a || !b) { return a | b; } u64 t, s = ctz64(a | b); a >>= ctz64(a); do { b >>= ctz64(b); if (a > b) { t = a, a = b, b = t; } b -= a; } while (b); return a << s; }

// https://ja.wikipedia.org/wiki/%E3%82%B3%E3%83%A0%E3%82%BD%E3%83%BC%E3%83%88#%E6%94%B9%E8%89%AF%E7%89%88%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0
void combsort11_32(int a_len, u32 *a) { int g = a_len; u32 t; while (true) { int flag = 1; g = (((g * 10) / 13) > 1) ? ((g * 10) / 13) : 1; if (g == 9 || g == 10) { g = 11; } for (int i = 0; i + g < a_len; i++) { if (a[i] > a[i + g]) { t = a[i], a[i] = a[i + g], a[i + g] = t; flag = 0; } } if (g == 1 && flag) { break; } } }
void combsort11_64(int a_len, u64 *a) { int g = a_len; u64 t; while (true) { int flag = 1; g = (((g * 10) / 13) > 1) ? ((g * 10) / 13) : 1; if (g == 9 || g == 10) { g = 11; } for (int i = 0; i + g < a_len; i++) { if (a[i] > a[i + g]) { t = a[i], a[i] = a[i + g], a[i + g] = t; flag = 0; } } if (g == 1 && flag) { break; } } }

// a * x + b * y = gcd(a, b)
typedef struct { i32 a, b; u32 d; } Bezout32;
Bezout32 bezout32(u32 x, u32 y) { u32 t; bool swap = x < y; if (swap) { t = x; x = y; y = t; } if (y == 0) { if (x == 0) { return (Bezout32){0, 0, 0}; } else if (swap) { return (Bezout32){0, 1, x}; } else { return (Bezout32){1, 0, x}; } } i32 s0 = 1, s1 = 0, t0 = 0, t1 = 1; while (true) { u32 q = x / y, r = x % y; if (r == 0) { if (swap) { return (Bezout32){t1, s1, y}; } else { return (Bezout32){s1, t1, y}; } } i32 s2 = s0 - (i32)(q) * s1, t2 = t0 - (i32)(q) * t1; x = y, y = r; s0 = s1, s1 = s2, t0 = t1, t1 = t2; } }
u32 modinv32(u32 x, u32 mod) { Bezout32 abd = bezout32(x, mod); return abd.a < 0 ? mod + abd.a : (u32)abd.a; }
typedef struct { i64 a, b; u64 d; } Bezout64;
Bezout64 bezout64(u64 x, u64 y) { u64 t; bool swap = x < y; if (swap) { t = x; x = y; y = t; } if (y == 0) { if (x == 0) { return (Bezout64){0, 0, 0}; } else if (swap) { return (Bezout64){0, 1, x}; } else { return (Bezout64){1, 0, x}; } } i64 s0 = 1, s1 = 0, t0 = 0, t1 = 1; while (true) { u64 q = x / y, r = x % y; if (r == 0) { if (swap) { return (Bezout64){t1, s1, y}; } else { return (Bezout64){s1, t1, y}; } } i64 s2 = s0 - (i64)(q) * s1, t2 = t0 - (i64)(q) * t1; x = y, y = r; s0 = s1, s1 = s2, t0 = t1, t1 = t2; } }
u64 modinv64(u64 x, u64 mod) { Bezout64 abd = bezout64(x, mod); return abd.a < 0 ? mod + abd.a : (u64)abd.a; }

// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
/* mod < 1073741824 */
static u32 n_m32, n2_m32, ni_m32, r1_m32, r2_m32, r3_m32;
void set_m32(u32 mod) { if (mod == n_m32) { return; } n_m32 = mod; n2_m32 = mod << 1; ni_m32 = mod; ni_m32 *= 2 - ni_m32 * mod; ni_m32 *= 2 - ni_m32 * mod; ni_m32 *= 2 - ni_m32 * mod; ni_m32 *= 2 - ni_m32 * mod; r1_m32 = (u32)(i32)-1 % mod + 1; r2_m32 = (u64)(i64)-1 % mod + 1; r3_m32 = (u32)(((u64)r1_m32 * (u64)r2_m32) % mod); }
u32 mr32(u64 a) { u32 y = (u32)(a >> 32) - (u32)(((u64)((u32)a * ni_m32) * n_m32) >> 32); return (i32)y < 0 ? y + n_m32 : y; }
u32 to_m32(u32 a) { return mr32((u64)a * r2_m32); }
u32 from_m32(u32 a) { return mr32((u64)a); }
u32 add_m32(u32 a, u32 b) { a += b; a -= (a >= n_m32 ? n_m32 : 0); return a; }
u32 sub_m32(u32 a, u32 b) { a += (a < b ? n_m32 : 0); a -= b; return a; }
u32 min_m32(u32 a) { return sub_m32(0, a); }
u32 mul_m32(u32 a, u32 b) { return mr32((u64)a * b); }
u32 squ_m32(u32 a) { return mr32((u64)a * a); }
u32 shl_m32(u32 a) { return (a <<= 1) >= n_m32 ? a - n_m32 : a; }
u32 shr_m32(u32 a) { return (a & 1) ? ((a >> 1) + (n_m32 >> 1) + 1) : (a >> 1); }
u32 pow_m32(u32 a, u32 k) { u32 ret = r1_m32, deg = k; while (deg > 0) { if (deg & 1) { ret = mul_m32(ret, a); } a = squ_m32(a); deg >>= 1; } return ret; }
u32 inv_m32(u32 a) { return mr32((u64)r3_m32 * modinv32(a, n_m32)); }
u32 div_m32(u32 a, u32 b) { return mul_m32(a, inv_m32(b)); }

/* mod < 4611686018427387904 */
static u64 n_m64, n2_m64, ni_m64, r1_m64, r2_m64, r3_m64;
void set_m64(u64 mod) { if (mod == n_m64) { return; } n_m64 = mod; n2_m64 = mod << 1; ni_m64 = mod; ni_m64 *= 2 - ni_m64 * mod; ni_m64 *= 2 - ni_m64 * mod; ni_m64 *= 2 - ni_m64 * mod; ni_m64 *= 2 - ni_m64 * mod; ni_m64 *= 2 - ni_m64 * mod; r1_m64 = (u64)(i64)-1 % mod + 1; r2_m64 = (u128)(i128)-1 % mod + 1; r3_m64 = (u64)(((u128)r1_m64 * (u128)r2_m64) % mod); }
u64 mr64(u128 a) { u64 y = (u64)(a >> 64) - (u64)(((u128)((u64)a * ni_m64) * n_m64) >> 64); return (i64)y < 0 ? y + n_m64 : y; }
u64 to_m64(u64 a) { return mr64((u128)a * r2_m64); }
u64 from_m64(u64 a) { return mr64((u128)a); }
u64 add_m64(u64 a, u64 b) { a += b; a -= (a >= n_m64 ? n_m64 : 0); return a; }
u64 sub_m64(u64 a, u64 b) { a += (a < b ? n_m64 : 0); a -= b; return a; }
u64 min_m64(u64 a) { return sub_m64(0, a); }
u64 mul_m64(u64 a, u64 b) { return mr64((u128)a * b); }
u64 squ_m64(u64 a) { return mr64((u128)a * a); }
u64 shl_m64(u64 a) { return (a <<= 1) >= n_m64 ? a - n_m64 : a; }
u64 shr_m64(u64 a) { return (a & 1) ? ((a >> 1) + (n_m64 >> 1) + 1) : (a >> 1); }
u64 pow_m64(u64 a, u64 k) { u64 ret = r1_m64, deg = k; while (deg > 0) { if (deg & 1) { ret = mul_m64(ret, a); } a = squ_m64(a); deg >>= 1; } return ret; }
u64 inv_m64(u64 a) { return mr64((u128)r3_m64 * modinv64(a, n_m64)); }
u64 div_m64(u64 a, u64 b) { return mul_m64(a, inv_m64(b)); }

// https://en.wikipedia.org/wiki/Barrett_reduction
static u64 m_b32, im_b32, div_b32, rem_b32;
void set_b32(u64 mod) { m_b32 = mod; im_b32 = ((((u128)(1ull) << 64)) + mod - 1) / mod; div_b32 = 0; rem_b32 = 0; }
void br32(u64 a) { u64 x = (u64)(((u128)(a) * im_b32) >> 64); u64 y = x * m_b32; unsigned long long z; u32 w = __builtin_usubll_overflow(a, y, &z) ? m_b32 : 0; div_b32 = x; rem_b32 = z + w; }
u32 add_b32(u32 a, u32 b) { a += b; a -= (a >= (u32)m_b32 ? (u32)m_b32 : 0); return a; }
u32 sub_b32(u32 a, u32 b) { a += (a < b ? (u32)m_b32 : 0); a -= b; return a; }
u32 min_b32(u32 a) { return sub_b32(0, a); }
u32 mul_b32(u32 a, u32 b) { br32((u64)a * b); return (u32)rem_b32; }
u32 squ_b32(u32 a) { br32((u64)a * a); return (u32)rem_b32; }
u32 shl_b32(u32 a) { return (a <<= 1) >= m_b32 ? a - m_b32 : a; }
u32 shr_b32(u32 a) { return (a & 1) ? ((a >> 1) + (m_b32 >> 1) + 1) : (a >> 1); }
u32 pow_b32(u32 a, u32 k) { u32 ret = 1u, deg = k; while (deg > 0) { if (deg & 1) { ret = mul_b32(ret, a); } a = squ_b32(a); deg >>= 1; } return ret; }

static u128 m_b64, im_b64, div_b64, rem_b64;
void set_b64(u128 mod) { if (mod == m_b64) { return; } m_b64 = mod; im_b64 = (~((u128)0ull)) / mod; div_b64 = 0; rem_b64 = 0; }
void br64(u128 x) { if (m_b64 == 1) { div_b64 = x; rem_b64 = 0; return; } u8 f; u128 a = x >> 64; u128 b = x & 0xffffffffffffffffull; u128 c = im_b64 >> 64; u128 d = im_b64 & 0xffffffffffffffffull; u128 ac = a * c; u128 bd = (b * d) >> 64; u128 ad = a * d; u128 bc = b * c; f = (ad > ((u128)((i128)(-1L)) - bd)); bd += ad; ac += f; f = (bc > ((u128)((i128)(-1L)) - bd)); bd += bc; ac += f; u128 q = ac + (bd >> 64); u128 r = x - q * m_b64; if (m_b64 <= r) { r -= m_b64; q += 1; } div_b64 = q; rem_b64 = r; }
u64 add_b64(u64 a, u64 b) { a += b; a -= (a >= (u64)m_b64 ? (u64)m_b64 : 0); return a; }
u64 sub_b64(u64 a, u64 b) { a += (a < b ? (u64)m_b64 : 0); a -= b; return a; }
u64 min_b64(u64 a) { return sub_b64(0, a); }
u64 mul_b64(u64 a, u64 b) { br64((u128)a * b); return (u64)rem_b64; }
u64 squ_b64(u64 a) { br64((u128)a * a); return (u64)rem_b64; }
u64 shl_b64(u64 a) { return (a <<= 1) >= m_b64 ? a - m_b64 : a; }
u64 shr_b64(u64 a) { return (a & 1) ? ((a >> 1) + (m_b64 >> 1) + 1) : (a >> 1); }
u64 pow_b64(u64 a, u64 k) { u64 ret = 1ull, deg = k; while (deg > 0) { if (deg & 1) { ret = mul_b64(ret, a); } a = squ_b64(a); deg >>= 1; } return ret; }

// clang-format on
#pragma endregion template

#pragma region miller rabin
// clang-format off

bool miller_rabin_small(u32 n)
{
    set_m32(n);
    u32 s = ctz32(n - 1);
    u32 d = (n - 1) >> s;
    u32 bases[3] = {2u, 7u, 61u};
    for (int i = 0; i < 3; i++)
    {
        if (n <= bases[i]) { return true; }
        u32 a = pow_m32(to_m32(bases[i]), d);
        if (a == r1_m32) { continue; }
        u32 r = 1;
        while (a != n - r1_m32)
        {
            if (r == s) { return false; }
            a = squ_m32(a);
            r++;
        }
    }
    return true;
}
bool miller_rabin(u64 n)
{
    set_m64(n);
    u64 s = ctz64(n - 1);
    u64 d = (n - 1) >> s;
    u64 bases[7] = {2ull, 325ull, 9375ull, 28178ull, 450775ull, 9780504ull, 1795265022ull};
    for (int i = 0; i < 7; i++)
    {
        if (n <= bases[i]) { return true; }
        u64 a = pow_m64(to_m64(bases[i]), d);
        if (a == r1_m64) { continue; }
        u64 r = 1;
        while (a != n - r1_m64)
        {
            if (r == s) { return false; }
            a = squ_m64(a);
            r++;
        }
    }
    return true;
}

// n < 2^62
bool is_prime(u64 n)
{
    if (n < 64ull) { return (1ull << n) & 2891462833508853932ull; }
    if (!(n & 1)) { return false; }
    if (gcd64(n, 15) != 1ull) { return false; }
    if (n < 1073741824ull) { return miller_rabin_small((uint32_t)n); }
    return miller_rabin(n);
}

// clang-format on
#pragma endregion miller rabin

#pragma region factorize
// clang-format off

u32 pollard_brent32(u32 n, u32 c) { set_m32(n); const u32 one = r1_m32; const u32 two = to_m32(2u); const u32 cc = to_m32(c); const u32 m = 1ull << ((31 - clz32(n)) / 5); u32 x = one, y = two, z = one, q = one; u32 g = 1u; for (u32 r = 1; g == 1; r <<= 1) { x = y; for (u32 i = 0; i < r; ++i) { y = add_m32(squ_m32(y), cc); } for (u32 k = 0; k < r && g == 1u; k += m) { z = y; for (u32 _ = 0; _ < ((m < (r - k)) ? m : (r - k)); _++) { y = add_m32(squ_m32(y), cc); q = mul_m32(q, sub_m32(x, y)); } g = gcd32(from_m32(q), n); } } if (g == n) { do { z = add_m32(squ_m32(z), cc); g = gcd32(from_m32(sub_m32(x, z)), n); } while (g == 1); } return g; }
u64 pollard_brent64(u64 n, u64 c) { set_m64(n); const u64 one = r1_m64; const u64 two = to_m64(2ull); const u64 cc = to_m64(c); const u64 m = 1ull << ((63 - clz64(n)) / 5); u64 x = one, y = two, z = one, q = one; u64 g = 1ull; for (u64 r = 1; g == 1; r <<= 1) { x = y; for (u64 i = 0; i < r; ++i) { y = add_m64(squ_m64(y), cc); } for (u64 k = 0; k < r && g == 1ull; k += m) { z = y; for (u64 _ = 0; _ < ((m < (r - k)) ? m : (r - k)); _++) { y = add_m64(squ_m64(y), cc); q = mul_m64(q, sub_m64(x, y)); } g = gcd64(from_m64(q), n); } } if (g == n) { do { z = add_m64(squ_m64(z), cc); g = gcd64(from_m64(sub_m64(x, z)), n); } while (g == 1); } return g; }
typedef struct { u64 x, z; } MC;
u64 check(MC p) { return gcd64(from_m64(p.z), n_m64); }
static u64 a24;
MC create_curve_and_point(void) { while (true) { u64 a = msws_range(0, n_m64 - 1ull); u64 x = msws_range(0, n_m64 - 1ull); u64 m1 = r1_m64; u64 y2 = mul_m64(x, add_m64(mul_m64(x, add_m64(x, a)), m1)); if (jacobi_symbol(from_m64(y2), n_m64) == 1) { a24 = div_m64(add_m64(a, to_m64(2ull)), to_m64(4ull)); return (MC){x, m1}; } } }
MC add_MC(MC p, MC q, MC diff) { u64 u = mul_m64(sub_m64(p.x, p.z), add_m64(q.x, q.z)); u64 v = mul_m64(add_m64(p.x, p.z), sub_m64(q.x, q.z)); u64 upv = add_m64(u, v); u64 umv = sub_m64(u, v); u64 new_x = mul_m64(mul_m64(diff.z, upv), upv); u64 new_z = mul_m64(mul_m64(diff.x, umv), umv); return (MC){new_x, new_z}; }
MC dbl_MC(MC p) { u64 s = add_m64(p.x, p.z); u64 d = sub_m64(p.x, p.z); u64 s2 = squ_m64(s); u64 d2 = squ_m64(d); u64 t = sub_m64(s2, d2); u64 new_x = mul_m64(s2, d2); u64 new_z = mul_m64(t, add_m64(d2, mul_m64(a24, t))); return (MC){new_x, new_z}; }
MC mul_MC(MC p, size_t k) { MC p0 = p; MC p1 = dbl_MC(p); for (int b = bit_width64(k) - 2; b >= 0; --b) { MC t = add_MC(p1, p0, p); if ((k >> b) & 1) { p1 = dbl_MC(p1); p0 = t; } else { p0 = dbl_MC(p0); p1 = t; } } return p0; }
u64 ecm(u64 n) { for (int k = 2; k < bit_width64(n); ++k) { u64 r = floor_kth_root_integer(n, k); u64 pw = r; for (int i = 1; i < k; ++i) { pw *= r; } if (pw == n) { return r; } } u64 t; set_m64(n); u64 ecm_blocks[10] = { 5690199479092128000ull, 810162134158954261ull, 326580695497527083ull, 13784092967194631821ull, 1107997261359193637ull, 6532397423431938467ull, 96265407405451883ull, 260006624961107813ull, 707992818804600227ull, 22417030981ull }; while (true) { MC point = create_curve_and_point(); u64 f = 1ull; for (size_t block = 0; block < 10; ++block) { MC new_point = mul_MC(point, ecm_blocks[block]); f = check(new_point); if (f != 1) { if (f != n_m64) { return f; } else { break; } } point = new_point; } if (f == n_m64) { continue; } MC six = dbl_MC(add_MC(dbl_MC(point), point, point)); MC q0 = six; MC q1 = dbl_MC(six); for (int _ = 6; _ < 400; _ += 6) { q0 = add_MC(q1, six, q0); t = q0.x; q0.x = q1.x; q1.x = t; t = q0.z; q0.z = q1.z; q1.z = t; } u64 xprod = r1_m64; u64 x_norm = div_m64(point.x, point.z); for (int i = 396; i < 3000; i += 6) { xprod = mul_m64(xprod, sub_m64(q0.x, mul_m64(q0.z, x_norm))); if (i % 300 == 0) { f = gcd64(from_m64(xprod), n_m64); if (f != 1) { if (f != n_m64) { return f; } else { break; } } } q0 = add_MC(q1, six, q0); t = q0.x; q0.x = q1.x; q1.x = t; t = q0.z; q0.z = q1.z; q1.z = t; } if (f == 1) { f = gcd64(from_m64(xprod), n_m64); if (f != 1 && f != n_m64) { return f; } } } }
u64 find_prime_factor(u64 n) { if (is_prime(n)) { return n; } for (int _ = 0; _ < 200; ++_) { u64 m; if (n < 1073741824ull) { m = pollard_brent32((u32)n, pcg_range(1u, (u32)n - 1)); if (is_prime(m)) { return m; } } else if (n < 10000000000000000ull) { m = pollard_brent64(n, xrsr128p_range(1ull, n - 1)); if (is_prime(m)) { return m; } } else { m = ecm(n); if (is_prime(m)) { return m; } } n = m; } return -1; }
u64 *factorize(u64 n) { u64 *ret = (u64 *)calloc(65, sizeof(u64)); if (ret == NULL) { exit(EXIT_FAILURE); } ret[64] = 0; for (int i = 0; i < 64; ++i) { ret[i] = (u64)(i64)-1; } u64 s = ctz64(n); n >>= s; ret[64] += s; for (u64 i = 0; i < s; ++i) { ret[i] = 2; } for (u64 i = 3ull; i <= 100ull && i * i <= n; i += 2ull) { if (n % i == 0) { do { n /= i; ret[ret[64]++] = i; } while (n % i == 0); } } while (n > 1) { u64 p = find_prime_factor(n); do { n /= p; ret[ret[64]++] = p; } while (n % p == 0); } combsort11_64(64, ret); return ret; }
typedef struct { u64 p; u64 k; } Factors;
void combsort11_factors(int a_len, Factors *a) { int g = a_len; Factors t; while (true) { bool flag = 1; g = (((g * 10) / 13) > 1) ? ((g * 10) / 13) : 1; if (g == 9 || g == 10) { g = 11; } for (int i = 0; i + g < a_len; i++) { if (a[i].p > a[i + g].p) { t = a[i], a[i] = a[i + g], a[i + g] = t; flag = 0; } } if (g == 1 && flag) { break; } } }
Factors *factors(u64 n) { Factors *ret = (Factors *)calloc(17, sizeof(Factors)); if (ret == NULL) { exit(EXIT_FAILURE); } for (int i = 0; i < 17; i++) { ret[i] = (Factors){UINT64_MAX, UINT64_MAX}; } ret[16].k = 0; u64 s = ctz64(n); n >>= s; if (s != 0) { ret[0] = (Factors){2, s}; ret[16].k += 1; } for (u64 i = 3ull; i <= 100ull && i * i <= n; i += 2ull) { if (n % i == 0) { u64 cnt = 0; do { n /= i; cnt++; } while (n % i == 0); ret[ret[16].k++] = (Factors){i, cnt}; } } while (n > 1) { u64 p = find_prime_factor(n); u64 cnt = 0; do { n /= p; cnt++; } while (n % p == 0); ret[ret[16].k++] = (Factors){p, cnt}; } combsort11_factors(16, ret); return ret; }

// clang-format on
#pragma endregion factorize
