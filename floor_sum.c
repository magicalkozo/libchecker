#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("unroll-loops")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#pragma region fastio
// clang-format off

#define BUFFER_SIZE 1048576
char num[10000][4];
__attribute__((constructor)) void pre(void) { for (int i = 0; i < 10000; i++) { int n = i; for (int j = 3; j >= 0; j--) { num[i][j] = n % 10 | '0'; n /= 10; } } }
char *ibuf, obuf[BUFFER_SIZE], out[12];
int ibufi, outi, obufi;
__attribute__((constructor)) void _c() { struct stat sb; fstat(0, &sb); ibufi = sb.st_size; ibuf = (char *)mmap(0, ibufi, 1, 32769, 0, 0); }
inline void flush(void) { write(1, obuf, obufi); obufi = 0; }
void rd_char(char *x) { *x = *ibuf++; }
void rd_int(int *x) { char c; for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';) { *x = *x * 10 + (c & 15); } }
void rd_uint(unsigned *x) { char c; for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';) { *x = *x * 10 + (c & 15); } }
void rd_long(long long *x) { char c; for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';) { *x = *x * 10 + (c & 15); } }
void rd_ulong(unsigned long long *x) { char c; for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';) { *x = *x * 10 + (c & 15); } }
void wt_char(char c) { obuf[obufi++] = c; }
void wt_ulong(unsigned long long x) { if (obufi > BUFFER_SIZE - 32) { flush(); } if (x >= 1e16) { long long q0 = x / 100000000; int r0 = x % 100000000; int q1 = q0 / 100000000, r1 = q0 % 100000000; if (x >= 1e18) { memcpy(obuf + obufi, num[q1] + 1, 3); memcpy(obuf + obufi + 3, num[r1 / 10000], 4); memcpy(obuf + obufi + 7, num[r1 % 10000], 4); memcpy(obuf + obufi + 11, num[r0 / 10000], 4); memcpy(obuf + obufi + 15, num[r0 % 10000], 4); obufi += 19; } else if (x >= 1e17) { int q2 = (q1 * 103) >> 10; obuf[obufi] = q2 | '0'; obuf[obufi + 1] = (q1 - q2 * 10) | '0'; memcpy(obuf + obufi + 2, num[r1 / 10000], 4); memcpy(obuf + obufi + 6, num[r1 % 10000], 4); memcpy(obuf + obufi + 10, num[r0 / 10000], 4); memcpy(obuf + obufi + 14, num[r0 % 10000], 4); obufi += 18; } else { obuf[obufi] = q1 | '0'; memcpy(obuf + obufi + 1, num[r1 / 10000], 4); memcpy(obuf + obufi + 5, num[r1 % 10000], 4); memcpy(obuf + obufi + 9, num[r0 / 10000], 4); memcpy(obuf + obufi + 13, num[r0 % 10000], 4); obufi += 17; } } else { for (outi = 8; x >= 10000; outi -= 4) { memcpy(out + outi, num[x % 10000], 4); x /= 10000; } if (x >= 1000) { memcpy(obuf + obufi, num[x], 4); obufi += 4; } else if (x >= 100) { memcpy(obuf + obufi, num[x] + 1, 3); obufi += 3; } else if (x >= 10) { int q = (x * 103) >> 10; obuf[obufi] = q | '0'; obuf[obufi + 1] = (x - q * 10) | '0'; obufi += 2; } else { obuf[obufi++] = x | '0'; } memcpy(obuf + obufi, out + outi + 4, 8 - outi); obufi += 8 - outi; } }
void wt_long(long long x) { if (obufi > BUFFER_SIZE - 32) { flush(); } if (x < 0) { wt_char('-'); x = -x; } wt_ulong(x); }
void wt_uint(unsigned x) { if (obufi > BUFFER_SIZE - 32) { flush(); } for (outi = 8; x >= 10000; outi -= 4) { memcpy(out + outi, num[x % 10000], 4); x /= 10000; } if (x >= 1000) { memcpy(obuf + obufi, num[x], 4); obufi += 4; } else if (x >= 100) { memcpy(obuf + obufi, num[x] + 1, 3); obufi += 3; } else if (x >= 10) { int q = (x * 103) >> 10; obuf[obufi] = q | '0'; obuf[obufi + 1] = (x - q * 10) | '0'; obufi += 2; } else { obuf[obufi++] = x | '0'; } memcpy(obuf + obufi, out + outi + 4, 8 - outi); obufi += 8 - outi; }
void wt_int(int x) { if (obufi > BUFFER_SIZE - 32) { flush(); } if (x < 0) { wt_char('-'); x = -x; } wt_uint(x); }
__attribute__((destructor)) void _d() { flush(); }

// clang-format on
#pragma endregion fastio

long long floor_sum(int n, int m, int a, int b)
{
    long long ret = 0;
    if (a >= m)
    {
        ret += ((long long)n - 1) * n / 2 * (a / m);
        a %= m;
    }
    if (b >= m)
    {
        ret += (long long)(n) * (b / m);
        b %= m;
    }
    long long c = (long long)(a) * n + b;
    if (c < m)
    {
        return ret;
    }
    return ret + floor_sum(c / m, a, m, c % m);
}


int main(void)
{
    unsigned T;
    rd_uint(&T);
    int n, m, a, b;
    while (T--)
    {
        rd_int(&n);
        rd_int(&m);
        rd_int(&a);
        rd_int(&b);
        wt_long(floor_sum(n, m, a, b));
        wt_char('\n');
    }
    return 0;
}