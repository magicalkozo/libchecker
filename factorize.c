#include "template.h"

static const u32 base32[256] = {
    1216, 1836, 8885, 4564, 10978, 5228, 15613, 13941,
    1553, 173, 3615, 3144, 10065, 9259, 233, 2362,
    6244, 6431, 10863, 5920, 6408, 6841, 22124, 2290,
    45597, 6935, 4835, 7652, 1051, 445, 5807, 842,
    1534, 22140, 1282, 1733, 347, 6311, 14081, 11157,
    186, 703, 9862, 15490, 1720, 17816, 10433, 49185,
    2535, 9158, 2143, 2840, 664, 29074, 24924, 1035,
    41482, 1065, 10189, 8417, 130, 4551, 5159, 48886,
    786, 1938, 1013, 2139, 7171, 2143, 16873, 188,
    5555, 42007, 1045, 3891, 2853, 23642, 148, 3585,
    3027, 280, 3101, 9918, 6452, 2716, 855, 990,
    1925, 13557, 1063, 6916, 4965, 4380, 587, 3214,
    1808, 1036, 6356, 8191, 6783, 14424, 6929, 1002,
    840, 422, 44215, 7753, 5799, 3415, 231, 2013,
    8895, 2081, 883, 3855, 5577, 876, 3574, 1925,
    1192, 865, 7376, 12254, 5952, 2516, 20463, 186,
    5411, 35353, 50898, 1084, 2127, 4305, 115, 7821,
    1265, 16169, 1705, 1857, 24938, 220, 3650, 1057,
    482, 1690, 2718, 4309, 7496, 1515, 7972, 3763,
    10954, 2817, 3430, 1423, 714, 6734, 328, 2581,
    2580, 10047, 2797, 155, 5951, 3817, 54850, 2173,
    1318, 246, 1807, 2958, 2697, 337, 4871, 2439,
    736, 37112, 1226, 527, 7531, 5418, 7242, 2421,
    16135, 7015, 8432, 2605, 5638, 5161, 11515, 14949,
    748, 5003, 9048, 4679, 1915, 7652, 9657, 660,
    3054, 15469, 2910, 775, 14106, 1749, 136, 2673,
    61814, 5633, 1244, 2567, 4989, 1637, 1273, 11423,
    7974, 7509, 6061, 531, 6608, 1088, 1627, 160,
    6416, 11350, 921, 306, 18117, 1238, 463, 1722,
    996, 3866, 6576, 6055, 130, 24080, 7331, 3922,
    8632, 2706, 24108, 32374, 4237, 15302, 287, 2296,
    1220, 20922, 3350, 2089, 562, 11745, 163, 11951};


bool is_prime_m32(u32 n) {
    set_m32(n);
    u32 b = base32[(u32)(n * 2908904329) >> 24];
    u32 s = ctz32(n - 1);
    u32 d = (n - 1) >> s;
    u32 a = to_m32(b);
    u32 c = pow_m32(a, d);
    if (c == r1_m32)
        return true;
    for (u32 r = 0; r < s; r++) {
        if (c == n - r1_m32)
            return true;
        c = squ_m32(c);
    }
    return false;
}
bool is_prime_b32(u32 n) {
    set_b32(n);
    u32 b = base32[(u32)(n * 2908904329) >> 24];
    u32 s = ctz32(n - 1);
    u32 d = (n - 1) >> s;
    u32 c = pow_b32(b, d);
    if (c == 1)
        return true;
    for (u32 r = 0; r < s; r++) {
        if (c == n - 1)
            return true;
        c = squ_b32(c);
    }
    return false;
}
#define to_b64
#define BPSW_MILLER_TEST(one, shl, squ)                                                         \
    {                                                                                           \
        u64 d = (n - 1) << clz64(n - 1);                                                        \
        u64 t = shl(one);                                                                       \
        for (d <<= 1; d; d <<= 1) {                                                             \
            t = squ(t);                                                                         \
            if (d >> 63)                                                                        \
                t = shl(t);                                                                     \
        }                                                                                       \
        if (t != one) {                                                                         \
            u64 x = only_lsb(n - 1);                                                            \
            u64 r = n - one;                                                                    \
            for (x >>= 1; t != r; x >>= 1) {                                                    \
                if (x == 0)                                                                     \
                    return false;                                                               \
                t = squ(t);                                                                     \
            }                                                                                   \
        }                                                                                       \
    }
#define BPSW_SET_VAR(one, to_64)                                                                \
    i64 D = 5;                                                                                  \
    {                                                                                           \
        for (int i = 0; jacobi_symbol(D, n) != -1 && i < 64; ++i) {                             \
            if (i == 32) {                                                                      \
                u64 sqrt_n = (u64)sqrtl((long double)n);                                        \
                if (n == sqrt_n * sqrt_n)                                                       \
                    return false;                                                               \
            }                                                                                   \
            D = i & 1 ? 2 - D : -2 - D;                                                         \
        }                                                                                       \
    }                                                                                           \
    u64 Q, Qn, u, v, k;                                                                         \
    {                                                                                           \
        Q = to_64((D < 0) ? ((1 - D) / 4 % n) : (n - (D - 1) / 4 % n));                         \
        Qn = Q;                                                                                 \
        u = one;                                                                                \
        v = one;                                                                                \
        k = (n + 1) << clz64(n + 1);                                                            \
        D %= (i64)n;                                                                            \
        D = to_64((D < 0) ? (D + n) : D);                                                       \
    }
#define BPSW_REMOVE_S_LUCAS_PP(add, sub, mul, squ, shr)                                         \
    {                                                                                           \
        for (k <<= 1; k; k <<= 1) {                                                             \
            u = mul(u, v);                                                                      \
            v = sub(squ(v), add(Qn, Qn));                                                       \
            Qn = squ(Qn);                                                                       \
            if (k >> 63) {                                                                      \
                u64 uu = add(u, v);                                                             \
                uu = shr(uu);                                                                   \
                v = shr(add(mul(D, u), v));                                                     \
                u = uu;                                                                         \
                Qn = mul(Qn, Q);                                                                \
            }                                                                                   \
        }                                                                                       \
        if (u == 0 || v == 0)                                                                   \
            return true;                                                                        \
        u64 x = (n + 1) & ~n;                                                                   \
        for (x >>= 1; x; x >>= 1) {                                                             \
            u = mul(u, v);                                                                      \
            v = sub(squ(v), add(Qn, Qn));                                                       \
            if (v == 0)                                                                         \
                return true;                                                                    \
            Qn = squ(Qn);                                                                       \
        }                                                                                       \
    }

bool is_prime_m64(u64 n) {
    set_m64(n);
    BPSW_MILLER_TEST(r1_m64, shl_m64, squ_m64);
    BPSW_SET_VAR(r1_m64, to_m64);
    BPSW_REMOVE_S_LUCAS_PP(add_m64, sub_m64, mul_m64, squ_m64, shr_m64);
    return false;
}
bool is_prime_b64(u64 n) {
    set_b64(n);
    BPSW_MILLER_TEST(1ull, shl_b64, squ_b64);
    BPSW_SET_VAR(1, to_b64);
    BPSW_REMOVE_S_LUCAS_PP(add_b64, sub_b64, mul_b64, squ_b64, shr_b64);
    return false;
}
#undef BPSW_REMOVE_S_LUCAS_PP
#undef BPSW_SET_VAR
#undef BPSW_MILLER_TEST
#undef to_b64
bool is_prime(u64 n) {
    if (n < 64ull)
        return (1ull << n) & 2891462833508853932ull;
    if (!(n & 1))
        return false;
    if (gcd64(n, 15) != 1ull)
        return false;
    if (gcd64(n, 77) != 1ull)
        return false;
    if (n < 121ull)
        return true;
    if (n < 1073741824)
        return is_prime_m32((u32)n);
    if (n < 4294967296ull)
        return is_prime_b32((u32)n);
    if (n < 4611686018427387904ull)
        return is_prime_m64(n);
    return is_prime_b64(n);
}
#define POLLARD_BRENT_RHO(bit)                                                                  \
    set_m##bit(n);                                                                              \
    const u##bit one = r1_m##bit;                                                               \
    const u##bit two = to_m##bit(2);                                                            \
    const u##bit cc = to_m##bit(c);                                                             \
    const u##bit m = 1ull << ((((bit) - 1) - clz##bit(n)) / 5);                                 \
    u##bit x = one, y = two, z = one, q = one;                                                  \
    u##bit g = 1;                                                                               \
    for (u##bit r = 1; g == 1; r <<= 1) {                                                       \
        x = y;                                                                                  \
        for (u##bit i = 0; i < r; ++i) {                                                        \
            y = add_m##bit(squ_m##bit(y), cc);                                                  \
        }                                                                                       \
        for (u##bit k = 0; k < r && g == 1; k += m) {                                           \
            z = y;                                                                              \
            for (u##bit _ = 0; _ < ((m < (r - k)) ? m : (r - k)); _++) {                        \
                y = add_m##bit(squ_m##bit(y), cc);                                              \
                q = mul_m##bit(q, sub_m##bit(x, y));                                            \
            }                                                                                   \
            g = gcd##bit(from_m##bit(q), n);                                                    \
        }                                                                                       \
    }                                                                                           \
    if (g == n) {                                                                               \
        do {                                                                                    \
            z = add_m##bit(squ_m##bit(z), cc);                                                  \
            g = gcd##bit(from_m##bit(sub_m##bit(x, z)), n);                                     \
        } while (g == 1);                                                                       \
    }                                                                                           \
    return g;

u32 pollard_brent_rho32(u32 n, u32 c) { POLLARD_BRENT_RHO(32); }
u64 pollard_brent_rho64(u64 n, u64 c) { POLLARD_BRENT_RHO(64); }
#undef POLLARD_BRENT_RHO

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
            m = pollard_brent_rho32((u32)n, randrange_32(1u, (u32)n - 1));
            if (is_prime(m))
                return m;
        } else if (n < 10000000000000000ull) {
            m = pollard_brent_rho64(n, randrange_64(1ull, n - 1));
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
    combsort11_u64(min_sort_u64, 64, ret);
    return ret;
}
typedef Pair_u64 Factor;
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
    combsort11_Pair_u64(min_sort_fst_Pair_u64, 16, ret);
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
