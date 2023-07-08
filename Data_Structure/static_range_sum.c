#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("unroll-loops")

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#define BUFFER_SIZE 1048576

char num[10000][4];

__attribute__((constructor)) void pre(void)
{
    for (int i = 0; i < 10000; i++)
    {
        int n = i;
        for (int j = 3; j >= 0; j--)
        {
            num[i][j] = n % 10 | '0';
            n /= 10;
        }
    }
}
char *ibuf, obuf[BUFFER_SIZE], out[12];
int ibufi, outi, obufi;

__attribute__((constructor)) void _c()
{
    struct stat sb;
    fstat(0, &sb);
    ibufi = sb.st_size;
    ibuf = (char *)mmap(0, ibufi, 1, 32769, 0, 0);
}

inline void flush(void)
{
    write(1, obuf, obufi);
    obufi = 0;
}
void rd_char(char *x)
{
    *x = *ibuf++;
}
void rd_int(int *x)
{
    char c;
    for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';)
    {
        *x = *x * 10 + (c & 15);
    }
}
void rd_uint(unsigned *x)
{
    char c;
    for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';)
    {
        *x = *x * 10 + (c & 15);
    }
}
void rd_long(long long *x)
{
    char c;
    for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';)
    {
        *x = *x * 10 + (c & 15);
    }
}
void rd_ulong(unsigned long long *x)
{
    char c;
    for (*x = *ibuf++ & 15; (c = *ibuf++) >= '0';)
    {
        *x = *x * 10 + (c & 15);
    }
}
void wt_char(char c)
{
    obuf[obufi++] = c;
}
void wt_ulong(unsigned long long x)
{
    if (obufi > BUFFER_SIZE - 32)
    {
        flush();
    }
    if (x >= 1e16)
    {
        long long q0 = x / 100000000;
        int r0 = x % 100000000;
        int q1 = q0 / 100000000, r1 = q0 % 100000000;
        if (x >= 1e18)
        {
            memcpy(obuf + obufi, num[q1] + 1, 3);
            memcpy(obuf + obufi + 3, num[r1 / 10000], 4);
            memcpy(obuf + obufi + 7, num[r1 % 10000], 4);
            memcpy(obuf + obufi + 11, num[r0 / 10000], 4);
            memcpy(obuf + obufi + 15, num[r0 % 10000], 4);
            obufi += 19;
        }
        else if (x >= 1e17)
        {
            int q2 = (q1 * 103) >> 10;
            obuf[obufi] = q2 | '0';
            obuf[obufi + 1] = (q1 - q2 * 10) | '0';
            memcpy(obuf + obufi + 2, num[r1 / 10000], 4);
            memcpy(obuf + obufi + 6, num[r1 % 10000], 4);
            memcpy(obuf + obufi + 10, num[r0 / 10000], 4);
            memcpy(obuf + obufi + 14, num[r0 % 10000], 4);
            obufi += 18;
        }
        else
        {
            obuf[obufi] = q1 | '0';
            memcpy(obuf + obufi + 1, num[r1 / 10000], 4);
            memcpy(obuf + obufi + 5, num[r1 % 10000], 4);
            memcpy(obuf + obufi + 9, num[r0 / 10000], 4);
            memcpy(obuf + obufi + 13, num[r0 % 10000], 4);
            obufi += 17;
        }
    }
    else
    {
        for (outi = 8; x >= 10000; outi -= 4)
        {
            memcpy(out + outi, num[x % 10000], 4);
            x /= 10000;
        }
        if (x >= 1000)
        {
            memcpy(obuf + obufi, num[x], 4);
            obufi += 4;
        }
        else if (x >= 100)
        {
            memcpy(obuf + obufi, num[x] + 1, 3);
            obufi += 3;
        }
        else if (x >= 10)
        {
            int q = (x * 103) >> 10;
            obuf[obufi] = q | '0';
            obuf[obufi + 1] = (x - q * 10) | '0';
            obufi += 2;
        }
        else
        {
            obuf[obufi++] = x | '0';
        }
        memcpy(obuf + obufi, out + outi + 4, 8 - outi);
        obufi += 8 - outi;
    }
}
void wt_long(long long x)
{
    if (obufi > BUFFER_SIZE - 32)
    {
        flush();
    }
    if (x < 0)
    {
        wt_char('-');
        x = -x;
    }
    wt_ulong(x);
}
void wt_uint(unsigned x)
{
    if (obufi > BUFFER_SIZE - 32)
    {
        flush();
    }
    for (outi = 8; x >= 10000; outi -= 4)
    {
        memcpy(out + outi, num[x % 10000], 4);
        x /= 10000;
    }
    if (x >= 1000)
    {
        memcpy(obuf + obufi, num[x], 4);
        obufi += 4;
    }
    else if (x >= 100)
    {
        memcpy(obuf + obufi, num[x] + 1, 3);
        obufi += 3;
    }
    else if (x >= 10)
    {
        int q = (x * 103) >> 10;
        obuf[obufi] = q | '0';
        obuf[obufi + 1] = (x - q * 10) | '0';
        obufi += 2;
    }
    else
    {
        obuf[obufi++] = x | '0';
    }
    memcpy(obuf + obufi, out + outi + 4, 8 - outi);
    obufi += 8 - outi;
}
void wt_int(int x)
{
    if (obufi > BUFFER_SIZE - 32)
    {
        flush();
    }
    if (x < 0)
    {
        wt_char('-');
        x = -x;
    }
    wt_uint(x);
}

__attribute__((destructor)) void _d() { flush(); }

int main(void)
{
    unsigned long long A[500010];
    unsigned long long a;
    unsigned N, Q, l, r;
    rd_uint(&N);
    rd_uint(&Q);
    rd_ulong(&a);
    A[0] = a;
    for (unsigned i = 1; i < N; i++)
    {
        rd_ulong(&a);
        A[i] = a + A[i - 1];
    }
    while (Q--)
    {
        rd_uint(&l);
        rd_uint(&r);
        if (l == 0)
        {
            wt_ulong(A[r - 1]);
        }
        else
        {
            wt_ulong(A[r - 1] - A[l - 1]);
        }
        wt_char('\n');
    }
    return 0;
}