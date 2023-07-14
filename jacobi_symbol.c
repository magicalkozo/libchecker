int jacobi_symbol(long long a, long long n)
{
    int j = 1;
    long long t;
    while (a)
    {
        if (a < 0)
        {
            a = -a;
            if ((n & 3) == 3)
            {
                j = -j;
            }
        }
        int s = ctz64(a);
        a >>= s;
        if ((((n & 7) == 3) || ((n & 7) == 5)) && (s & 1))
        {
            j = -j;
        }
        if ((a & n & 3) == 3)
        {
            j = -j;
        }
        t = a, a = n, n = t;
        a %= n;
        if ((a << 1) > n)
        {
            a -= n;
        }
    }
    return n == 1 ? j : 0;
}