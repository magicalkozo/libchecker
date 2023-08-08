#include "template.h"

bool miller_rabin_u32(u32 n) {
    set_m32(n);
    u32 s = ctz32(n - 1);
    u32 d = (n - 1) >> s;
    u32 bases[3] = {2u, 7u, 61u};
    for (int i = 0; i < 3; i++) {
        if (n <= bases[i])
            return true;
        u32 a = pow_m32(to_m32(bases[i]), d);
        if (a == r1_m32)
            continue;
        u32 r = 1;
        while (a != n - r1_m32) {
            if (r == s)
                return false;
            a = squ_m32(a);
            r++;
        }
    }
    return true;
}
bool miller_rabin_u64(u64 n) {
    set_m64(n);
    u64 s = ctz64(n - 1);
    u64 d = (n - 1) >> s;
    u64 bases[7] = {2ull, 325ull, 9375ull, 28178ull, 450775ull, 9780504ull, 1795265022ull};
    for (int i = 0; i < 7; i++) {
        if (n <= bases[i])
            return true;
        u64 a = pow_m64(to_m64(bases[i]), d);
        if (a == r1_m64)
            continue;
        u64 r = 1;
        while (a != n - r1_m64) {
            if (r == s)
                return false;
            a = squ_m64(a);
            r++;
        }
    }
    return true;
}
bool is_prime(u64 n) {
    if (n < 64ull)
        return (1ull << n) & 2891462833508853932ull;
    if (!(n & 1))
        return false;
    if (gcd64(n, 15) != 1ull)
        return false;
    if (n < 1073741824ull)
        return miller_rabin_u32((u32)n);
    return miller_rabin_u64(n);
}
u32 pollard_brent32(u32 n, u32 c) {
    set_m32(n);
    const u32 one = r1_m32;
    const u32 two = to_m32(2u);
    const u32 cc = to_m32(c);
    const u32 m = 1ull << ((31 - clz32(n)) / 5);
    u32 x = one, y = two, z = one, q = one;
    u32 g = 1u;
    for (u32 r = 1; g == 1; r <<= 1) {
        x = y;
        for (u32 i = 0; i < r; ++i) {
            y = add_m32(squ_m32(y), cc);
        }
        for (u32 k = 0; k < r && g == 1u; k += m) {
            z = y;
            for (u32 _ = 0; _ < ((m < (r - k)) ? m : (r - k)); _++) {
                y = add_m32(squ_m32(y), cc);
                q = mul_m32(q, sub_m32(x, y));
            }
            g = gcd32(from_m32(q), n);
        }
    }
    if (g == n) {
        do {
            z = add_m32(squ_m32(z), cc);
            g = gcd32(from_m32(sub_m32(x, z)), n);
        } while (g == 1);
    }
    return g;
}
u64 pollard_brent64(u64 n, u64 c) {
    set_m64(n);
    const u64 one = r1_m64;
    const u64 two = to_m64(2ull);
    const u64 cc = to_m64(c);
    const u64 m = 1ull << ((63 - clz64(n)) / 5);
    u64 x = one, y = two, z = one, q = one;
    u64 g = 1ull;
    for (u64 r = 1; g == 1; r <<= 1) {
        x = y;
        for (u64 i = 0; i < r; ++i) {
            y = add_m64(squ_m64(y), cc);
        }
        for (u64 k = 0; k < r && g == 1ull; k += m) {
            z = y;
            for (u64 _ = 0; _ < ((m < (r - k)) ? m : (r - k)); _++) {
                y = add_m64(squ_m64(y), cc);
                q = mul_m64(q, sub_m64(x, y));
            }
            g = gcd64(from_m64(q), n);
        }
    }
    if (g == n) {
        do {
            z = add_m64(squ_m64(z), cc);
            g = gcd64(from_m64(sub_m64(x, z)), n);
        } while (g == 1);
    }
    return g;
}
typedef Pair_u64 MC;
u64 check(MC p) { return gcd64(from_m64(p.b), n_m64); }
static u64 a24;
MC create_curve_and_point(void) {
    while (true) {
        u64 a = randrange_64(0, n_m64 - 1ull);
        u64 x = randrange_64(0, n_m64 - 1ull);
        u64 m1 = r1_m64;
        u64 y2 = mul_m64(x, add_m64(mul_m64(x, add_m64(x, a)), m1));
        if (jacobi_symbol(from_m64(y2), n_m64) == 1) {
            a24 = div_m64(add_m64(a, to_m64(2ull)), to_m64(4ull));
            return (MC){x, m1};
        }
    }
}
MC add_MC(MC p, MC q, MC diff) {
    u64 u = mul_m64(sub_m64(p.a, p.b), add_m64(q.a, q.b));
    u64 v = mul_m64(add_m64(p.a, p.b), sub_m64(q.a, q.b));
    u64 upv = add_m64(u, v);
    u64 umv = sub_m64(u, v);
    u64 new_x = mul_m64(mul_m64(diff.b, upv), upv);
    u64 new_z = mul_m64(mul_m64(diff.a, umv), umv);
    return (MC){new_x, new_z};
}
MC dbl_MC(MC p) {
    u64 s = add_m64(p.a, p.b);
    u64 d = sub_m64(p.a, p.b);
    u64 s2 = squ_m64(s);
    u64 d2 = squ_m64(d);
    u64 t = sub_m64(s2, d2);
    u64 new_x = mul_m64(s2, d2);
    u64 new_z = mul_m64(t, add_m64(d2, mul_m64(a24, t)));
    return (MC){new_x, new_z};
}
MC mul_MC(MC p, size_t k) {
    MC p0 = p;
    MC p1 = dbl_MC(p);
    for (int b = bit_width64(k) - 2; b >= 0; --b) {
        MC t = add_MC(p1, p0, p);
        if ((k >> b) & 1) {
            p1 = dbl_MC(p1);
            p0 = t;
        } else {
            p0 = dbl_MC(p0);
            p1 = t;
        }
    }
    return p0;
}
u64 ecm(u64 n) {
    for (int k = 2; k < bit_width64(n); ++k) {
        u64 r = floor_kth_root_integer(n, k);
        u64 pw = r;
        for (int i = 1; i < k; ++i)
            pw *= r;
        if (pw == n)
            return r;
    }
    u64 t;
    set_m64(n);
    u64 ecm_blocks[10] = {5690199479092128000ull, 810162134158954261ull, 326580695497527083ull, 13784092967194631821ull, 1107997261359193637ull, 6532397423431938467ull, 96265407405451883ull, 260006624961107813ull, 707992818804600227ull, 22417030981ull};
    while (true) {
        MC point = create_curve_and_point();
        u64 f = 1ull;
        for (size_t block = 0; block < 10; ++block) {
            MC new_point = mul_MC(point, ecm_blocks[block]);
            f = check(new_point);
            if (f != 1) {
                if (f != n_m64)
                    return f;
                else
                    break;
            }
            point = new_point;
        }
        if (f == n_m64)
            continue;
        MC six = dbl_MC(add_MC(dbl_MC(point), point, point));
        MC q0 = six;
        MC q1 = dbl_MC(six);
        for (int _ = 6; _ < 400; _ += 6) {
            q0 = add_MC(q1, six, q0);
            t = q0.a;
            q0.a = q1.a;
            q1.a = t;
            t = q0.b;
            q0.b = q1.b;
            q1.b = t;
        }
        u64 xprod = r1_m64;
        u64 x_norm = div_m64(point.a, point.b);
        for (int i = 396; i < 3000; i += 6) {
            xprod = mul_m64(xprod, sub_m64(q0.a, mul_m64(q0.b, x_norm)));
            if (i % 300 == 0) {
                f = gcd64(from_m64(xprod), n_m64);
                if (f != 1) {
                    if (f != n_m64)
                        return f;
                    else
                        break;
                }
            }
            q0 = add_MC(q1, six, q0);
            t = q0.a;
            q0.a = q1.a;
            q1.a = t;
            t = q0.b;
            q0.b = q1.b;
            q1.b = t;
        }
        if (f == 1) {
            f = gcd64(from_m64(xprod), n_m64);
            if (f != 1 && f != n_m64)
                return f;
        }
    }
}
u64 find_prime_factor(u64 n) {
    if (is_prime(n))
        return n;
    for (int _ = 0; _ < 200; ++_) {
        u64 m;
        if (n < 1073741824ull) {
            m = pollard_brent32((u32)n, randrange_32(1u, (u32)n - 1));
            if (is_prime(m))
                return m;
        } else if (n < 10000000000000000ull) {
            m = pollard_brent64(n, randrange_64(1ull, n - 1));
            if (is_prime(m))
                return m;
        } else {
            m = ecm(n);
            if (is_prime(m))
                return m;
        }
        n = m;
    }
    return -1;
}
u64 *factorize(u64 n) {
    u64 *ret = (u64 *)calloc(65, sizeof(u64));
    if (ret == NULL)
        exit(EXIT_FAILURE);
    ret[64] = 0;
    for (int i = 0; i < 64; ++i)
        ret[i] = (u64)(i64)-1;
    u64 s = ctz64(n);
    n >>= s;
    ret[64] += s;
    for (u64 i = 0; i < s; ++i)
        ret[i] = 2;
    for (u64 i = 3ull; i <= 100ull && i * i <= n; i += 2ull) {
        if (n % i == 0) {
            do {
                n /= i;
                ret[ret[64]++] = i;
            } while (n % i == 0);
        }
    }
    while (n > 1) {
        u64 p = find_prime_factor(n);
        do {
            n /= p;
            ret[ret[64]++] = p;
        } while (n % p == 0);
    }
    combsort11_64(64, ret);
    return ret;
}
typedef Pair_u64 Factor;
void combsort11_factors(int a_len, Factor *a) {
    int g = a_len;
    Factor t;
    while (true) {
        bool flag = 1;
        g = (((g * 10) / 13) > 1) ? ((g * 10) / 13) : 1;
        if (g == 9 || g == 10)
            g = 11;
        for (int i = 0; i + g < a_len; i++) {
            if (a[i].a > a[i + g].a) {
                t = a[i], a[i] = a[i + g], a[i + g] = t;
                flag = 0;
            }
        }
        if (g == 1 && flag)
            break;
    }
}
Factor *factors(u64 n) {
    Factor *ret = (Factor *)calloc(17, sizeof(Factor));
    if (ret == NULL)
        exit(EXIT_FAILURE);
    for (int i = 0; i < 17; i++)
        ret[i] = (Factor){UINT64_MAX, UINT64_MAX};
    ret[16].b = 0;
    u64 s = ctz64(n);
    n >>= s;
    if (s != 0) {
        ret[0] = (Factor){2, s};
        ret[16].b += 1;
    }
    for (u64 i = 3ull; i <= 100ull && i * i <= n; i += 2ull) {
        if (n % i == 0) {
            u64 cnt = 0;
            do {
                n /= i;
                cnt++;
            } while (n % i == 0);
            ret[ret[16].b++] = (Factor){i, cnt};
        }
    }
    while (n > 1) {
        u64 p = find_prime_factor(n);
        u64 cnt = 0;
        do {
            n /= p;
            cnt++;
        } while (n % p == 0);
        ret[ret[16].b++] = (Factor){p, cnt};
    }
    combsort11_factors(16, ret);
    return ret;
}

void solve_factorize(void) {
    int Q;
    rd_int(&Q);
    while (Q--) {
        u64 x;
        rd_u64(&x);
        Factor *pf = factors(x);
        u64 cnt = 0;
        for (u64 u = 0; u < pf[16].b; u++)
            cnt += pf[u].b;
        wt_u64(cnt);
        for (u64 i = 0; i < pf[16].b; i++) {
            for (u64 u = 0; u < pf[i].b; u++) {
                wt_sp();
                wt_u64(pf[i].a);
            }
        }
        wt_nl();
    }
}
