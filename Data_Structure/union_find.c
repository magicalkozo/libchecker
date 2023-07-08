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

int root(int x, size_t uf_size, int *union_find_tree)
{
    if (union_find_tree[x] < 0)
    {
        return x;
    }
    else
    {
        union_find_tree[x] = root(union_find_tree[x], uf_size, union_find_tree);
        return union_find_tree[x];
    }
}
bool is_same(int x, int y, size_t uf_size, int *union_find_tree)
{
    return root(x, uf_size, union_find_tree) == root(y, uf_size, union_find_tree);
}
bool merge(int x, int y, size_t uf_size, int *union_find_tree)
{
    x = root(x, uf_size, union_find_tree);
    y = root(y, uf_size, union_find_tree);
    if (x == y)
    {
        return false;
    }
    if (union_find_tree[x] > union_find_tree[y])
    {
        int t;
        t = x;
        x = y;
        y = t;
    }
    union_find_tree[x] += union_find_tree[y];
    union_find_tree[y] = x;
    return true;
}
int size(int x, size_t uf_size, int *union_find_tree)
{
    return -union_find_tree[root(x, uf_size, union_find_tree)];
}

int main(void)
{
    unsigned N, Q, t, u, v;
    rd_uint(&N);
    rd_uint(&Q);
    int *uf = (int *)malloc(sizeof(int) * N);
    for (unsigned i = 0; i < N; i++)
    {
        uf[i] = -1;
    }
    while (Q--)
    {
        rd_uint(&t);
        rd_uint(&u);
        rd_uint(&v);
        if (t == 0)
        {
            merge(u, v, N, uf);
        }
        else
        {
            wt_char(is_same(u, v, N, uf) ? '1' : '0');
            wt_char('\n');
        }
    }
    free(uf);
    return 0;
}