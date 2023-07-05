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

#define minimum(a, b) (((a) < (b)) ? (a) : (b))
#define maximum(a, b) (((a) > (b)) ? (a) : (b))
#define swap(a, b) ((a) ^= (b) ^= (a) ^= (b))

#define clz32(a) ((a) ? __builtin_clz((a)) : 32)
#define clz64(a) ((a) ? __builtin_clzll((a)) : 64)
#define ctz32(a) ((a) ? __builtin_ctz((a)) : 32)
#define ctz64(a) ((a) ? __builtin_ctzll((a)) : 64)
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


int main(void)
{
    return 0;
}