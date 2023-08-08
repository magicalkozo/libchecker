#include "template.h"

const u32 n1 = 998244353u;
const u32 ni_1 = 3296722945u;
const u32 n2_1 = 1996488706u;
const u32 r1_1 = 301989884u;
const u32 r2_1 = 932051910u;
const u32 r3_1 = 679058953u;
const u32 ga1[] = {691295370, 307583142, 566821959, 878217029, 375146819, 138254384, 500602490, 79119218, 790898700, 978335284, 651424567, 308706579, 723000027, 474797508, 683394121, 44141573, 536892010, 945865189, 175417726, 536169764, 831722880, 721458245};
const u32 iga1[] = {306948983, 888603487, 138723248, 65668869, 842568658, 953245971, 195169681, 118717521, 792052763, 828450244, 908724728, 218560432, 628507989, 248210924, 566568154, 6285593, 82571768, 49985074, 225413092, 349167278, 61514562, 763211248};

const u32 n2 = 985661441u;
const u32 ni_2 = 3309305857u;
const u32 n2_2 = 1971322882u;
const u32 r1_2 = 352321532u;
const u32 r2_2 = 616455619u;
const u32 r3_2 = 862208618u;
const u32 ga2[] = {969414155, 240156868, 716651500, 728800531, 977177032, 47314842, 240475723, 876076444, 626710676, 365360170, 808202916, 560909592, 755542104, 303317332, 75348256, 259192271, 882296372, 620044766, 876870197, 256206930, 761331788};
const u32 iga2[] = {16247286, 67104299, 325946810, 44505332, 582782266, 729124870, 724673072, 173952869, 594582867, 76943556, 66752559, 892797276, 469283465, 123325105, 933929770, 911329874, 741246559, 905930185, 828158135, 9523962, 198022420};

const u32 n3 = 943718401u;
const u32 ni_3 = 3351248897u;
const u32 n2_3 = 1887436802u;
const u32 r1_3 = 520093692u;
const u32 r2_3 = 917135855u;
const u32 r3_3 = 628490905u;
const u32 ga3[] = {125689310, 401270432, 193546243, 204233475, 765072983, 793690592, 598110941, 560814539, 323055569, 635997590, 661263945, 671645950, 596439462, 577210208, 667936112, 172603057, 698142776, 3390265, 400541812, 419143563, 100582761};
const u32 iga3[] = {818029091, 177917178, 278610320, 675646939, 629165784, 803573782, 552038920, 685763768, 343497720, 610893888, 604907871, 366961343, 132493990, 882172703, 730481417, 529389095, 864269596, 777879390, 446333578, 468025435, 879098724};

const u32 n4 = 935329793u;
const u32 ni_4 = 3359637505u;
const u32 n2_4 = 1870659586u;
const u32 r1_4 = 553648124u;
const u32 r2_4 = 53961717u;
const u32 r3_4 = 473631956u;
const u32 ga4[] = {127754847, 332035729, 485487390, 767812069, 338537822, 908760351, 266422919, 15870799, 282786169, 336464639, 917837900, 600157874, 102745434, 499721765, 297856120, 400728174, 874566324, 135528821, 283962647, 427476450, 373468654};
const u32 iga4[] = {807574946, 290485311, 409574102, 869639944, 58356144, 27374476, 879628632, 793485894, 827164904, 191931428, 683432311, 448322044, 311701814, 503088307, 417978759, 513126193, 721540170, 517773789, 633460269, 360340898, 201778657};

const u32 n5 = 918552577u;
const u32 ni_5 = 3376414721u;
const u32 n2_5 = 1837105154u;
const u32 r1_5 = 620756988u;
const u32 r2_5 = 394187990u;
const u32 r3_5 = 593746783u;
const u32 ga5[] = {812999284, 432990384, 554323451, 509021345, 537430786, 537202372, 408282163, 905806130, 114600116, 74220439, 529612107, 493978550, 597815913, 752217987, 27017876, 659065034, 609107995, 447505193, 917415928, 859609177, 789582125};
const u32 iga5[] = {105553293, 141994416, 219417454, 67374992, 569187052, 228827437, 78128169, 611730064, 263847616, 903474251, 229028647, 181422347, 438135835, 126700759, 646908306, 338522933, 183400791, 871889148, 412716320, 221223285, 914195002};

const u32 n6 = 1045430273u;
const u32 ni_6 = 3249537025u;
const u32 n2_6 = 2090860546u;
const u32 r1_6 = 113246204u;
const u32 r2_6 = 798281873u;
const u32 r3_6 = 772096312u;
const u32 ga6[] = {885458361u, 900924659u, 962059782u, 415593171u, 1042963020u, 670679244u, 541100730u, 791232404u, 1017138063u, 1009790060u, 862479644u, 14783221u, 244867424u, 431174616u, 403070125u, 826977693u, 130509512u, 234433839u, 583350782u};
const u32 iga6[] = {159971912u, 22701698u, 360170016u, 9164511u, 967212496u, 707952636u, 928400066u, 788102537u, 248599354u, 954265785u, 539282025u, 749013941u, 66727014u, 835563346u, 909028729u, 532555753u, 488995841u, 863413679u, 575546441u};

const u32 n7 = 1051721729u;
const u32 ni_7 = 3243245569u;
const u32 n2_7 = 2103443458u;
const u32 r1_7 = 88080380u;
const u32 r2_7 = 748646691u;
const u32 r3_7 = 478513284u;
const u32 ga7[] = {467413740u, 850027207u, 242082110u, 469704687u, 798875242u, 948071666u, 405859818u, 794843055u, 771520840u, 176008702u, 644573208u, 715864305u, 237085256u, 357934180u, 654888890u, 961560509u, 740357787u, 329663983u, 679393590u};
const u32 iga7[] = {584307989u, 414819672u, 727414907u, 1042290158u, 763418088u, 582659108u, 106169741u, 486228312u, 106468254u, 164241005u, 856035592u, 117294265u, 988587910u, 117912639u, 562619904u, 608072619u, 298409394u, 523681378u, 94283612u};

const u32 n8 = 1053818881u;
const u32 ni_8 = 3241148417u;
const u32 n2_8 = 2107637762u;
const u32 r1_8 = 79691772u;
const u32 r2_8 = 159648582u;
const u32 r3_8 = 1029151498u;
const u32 ga8[] = {882957072u, 871538451u, 673579336u, 73537272u, 425574675u, 697950553u, 938595210u, 127556594u, 879998874u, 1019737367u, 353383985u, 36465421u, 901066455u, 93757427u, 979102987u, 501430327u, 40143827u, 92524426u, 698947874u};
const u32 iga8[] = {170861809u, 516220135u, 927824564u, 603681864u, 586404698u, 720550924u, 639641177u, 182451670u, 257278119u, 437349354u, 600413897u, 648372120u, 399254828u, 182276433u, 483901237u, 893329600u, 390180935u, 919465749u, 333795344u};

#define NTT(index)                                                              \
    int h = 0;                                                                  \
    while (A_len > (1 << h))                                                    \
        h++;                                                                    \
    for (int ph = 1; ph <= h; ph++) {                                           \
        int w = 1 << (ph - 1);                                                  \
        int p = 1 << (h - ph);                                                  \
        u32 now = r1_##index;                                                   \
        for (int s = 0; s < w; s++) {                                           \
            int offset = s << (h - ph + 1);                                     \
            for (int i = 0; i < p; i++) {                                       \
                u32 l = A[i + offset];                                          \
                u32 r = mul_dm32(A[i + offset + p], now, ni_##index, n##index); \
                A[i + offset] = relaxed_add_dm32(l, r, n2_##index);             \
                A[i + offset + p] = relaxed_sub_dm32(l, r, n2_##index);         \
            }                                                                   \
            now = mul_dm32(now, ga##index[ctz32(~s)], ni_##index, n##index);    \
        }                                                                       \
    }

void ntt1(int A_len, u32 *A) { NTT(1); }
void ntt2(int A_len, u32 *A) { NTT(2); }
void ntt3(int A_len, u32 *A) { NTT(3); }
void ntt4(int A_len, u32 *A) { NTT(4); }
void ntt5(int A_len, u32 *A) { NTT(5); }
void ntt6(int A_len, u32 *A) { NTT(6); }
void ntt7(int A_len, u32 *A) { NTT(7); }
void ntt8(int A_len, u32 *A) { NTT(8); }
#undef NTT

#define INTT(index)                                                                                           \
    int h = 0;                                                                                                \
    while (A_len > (1 << h))                                                                                  \
        h++;                                                                                                  \
    for (int ph = h; ph >= 1; ph--) {                                                                         \
        int w = 1 << (ph - 1);                                                                                \
        int p = 1 << (h - ph);                                                                                \
        u32 inow = r1_##index;                                                                                \
        for (int s = 0; s < w; s++) {                                                                         \
            int offset = s << (h - ph + 1);                                                                   \
            for (int i = 0; i < p; i++) {                                                                     \
                u32 l = A[i + offset];                                                                        \
                u32 r = A[i + offset + p];                                                                    \
                A[i + offset] = relaxed_add_dm32(l, r, n2_##index);                                           \
                A[i + offset + p] = mul_dm32(relaxed_sub_dm32(l, r, n2_##index), inow, ni_##index, n##index); \
            }                                                                                                 \
            inow = mul_dm32(inow, iga##index[ctz32(~s)], ni_##index, n##index);                               \
        }                                                                                                     \
    }                                                                                                         \
    u32 inv2t = inv_dm32(to_dm32(A_len, r2_##index, ni_##index, n##index), r3_##index, ni_##index, n##index); \
    for (int i = 0; i < A_len; i++)                                                                           \
        A[i] = mul_dm32(A[i], inv2t, ni_##index, n##index), ni_##index, n##index;

void intt1(int A_len, u32 *A) { INTT(1); }
void intt2(int A_len, u32 *A) { INTT(2); }
void intt3(int A_len, u32 *A) { INTT(3); }
void intt4(int A_len, u32 *A) { INTT(4); }
void intt5(int A_len, u32 *A) { INTT(5); }
void intt6(int A_len, u32 *A) { INTT(6); }
void intt7(int A_len, u32 *A) { INTT(7); }
void intt8(int A_len, u32 *A) { INTT(8); }
#undef INTT

#define CONVOLUTE(index)                                   \
    int ret_len = bit_ceil32(A_len + B_len - 1);           \
    u32 *C = (u32 *)calloc(ret_len, sizeof(u32));          \
    u32 *D = (u32 *)calloc(ret_len, sizeof(u32));          \
    if (C == NULL || D == NULL)                            \
        exit(EXIT_FAILURE);                                \
    memcpy(C, A, sizeof(u32) * A_len);                     \
    memcpy(D, B, sizeof(u32) * B_len);                     \
    ntt##index(ret_len, C);                                \
    ntt##index(ret_len, D);                                \
    for (int i = 0; i < ret_len; i++)                      \
        C[i] = mul_dm32(C[i], D[i], ni_##index, n##index); \
    free(D);                                               \
    intt##index(ret_len, C);                               \
    return C

u32 *conv1(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(1); }
u32 *conv2(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(2); }
u32 *conv3(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(3); }
u32 *conv4(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(4); }
u32 *conv5(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(5); }
u32 *conv6(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(6); }
u32 *conv7(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(7); }
u32 *conv8(int A_len, int B_len, u32 *A, u32 *B) { CONVOLUTE(8); }
#undef CONVOLUTE

u32 *conv_mod_arbitrary(u32 mod, int a_len, int b_len, u32 *a, u32 *b) {
    int ret_len = a_len + b_len - 1;

    u32 n0 = mod;
    u32 ni_0 = ni_dm32(mod);
    u32 n2_0 = n2_dm32(mod);
    u32 r1_0 = r1_dm32(mod);
    u32 r2_0 = r2_dm32(mod);
    u32 r3_0 = r3_dm32(mod, r1_0, r2_0);

    u32 *A1 = (u32 *)calloc(a_len, sizeof(u32));
    u32 *A2 = (u32 *)calloc(a_len, sizeof(u32));
    u32 *A3 = (u32 *)calloc(a_len, sizeof(u32));
    u32 *B1 = (u32 *)calloc(b_len, sizeof(u32));
    u32 *B2 = (u32 *)calloc(b_len, sizeof(u32));
    u32 *B3 = (u32 *)calloc(b_len, sizeof(u32));
    if (A1 == NULL || A2 == NULL || A3 == NULL || B1 == NULL || B2 == NULL || B3 == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < a_len; i++) {
        A1[i] = to_dm32(a[i], r2_6, ni_6, n6);
        A2[i] = to_dm32(a[i], r2_7, ni_7, n7);
        A3[i] = to_dm32(a[i], r2_8, ni_8, n8);
    }
    for (int i = 0; i < b_len; i++) {
        B1[i] = to_dm32(b[i], r2_6, ni_6, n6);
        B2[i] = to_dm32(b[i], r2_7, ni_7, n7);
        B3[i] = to_dm32(b[i], r2_8, ni_8, n8);
    }
    u32 *x = conv6(a_len, b_len, A1, B1);
    u32 *y = conv7(a_len, b_len, A2, B2);
    u32 *z = conv8(a_len, b_len, A3, B3);

    free(A1);
    free(A2);
    free(A3);
    free(B1);
    free(B2);
    free(B3);

    for (int i = 0; i < ret_len; ++i)
        x[i] = from_dm32(x[i], ni_6, n6);
    for (int i = 0; i < ret_len; ++i)
        y[i] = from_dm32(y[i], ni_7, n7);
    for (int i = 0; i < ret_len; ++i)
        z[i] = from_dm32(z[i], ni_8, n8);

    u32 m1_inv_m2 = inv_dm32(to_dm32(n6, r2_7, ni_7, n7), r3_7, ni_7, n7);
    u32 m12_inv_m3 = inv_dm32(mul_dm32(to_dm32(n6, r2_8, ni_8, n8), to_dm32(n7, r2_8, ni_8, n8), ni_8, n8), r3_8, ni_8, n8);
    u32 m1_m3 = to_dm32(n6, r2_8, ni_8, n8);
    u32 m1_m0 = to_dm32(n6, r2_0, ni_0, n0);
    u32 m12_m0 = mul_dm32(to_dm32(n6, r2_0, ni_0, n0), to_dm32(n7, r2_0, ni_0, n0), ni_0, n0);

    u32 *ret = (u32 *)calloc(ret_len, sizeof(u32));
    if (ret == NULL)
        exit(EXIT_FAILURE);

    for (u32 i = 0; i < ret_len; i++) {
        u32 xi_m2 = to_dm32(x[i], r2_7, ni_7, n7);
        u32 yi_m2 = to_dm32(y[i], r2_7, ni_7, n7);
        u32 zi_m3 = to_dm32(z[i], r2_8, ni_8, n8);
        u32 xi_m3 = to_dm32(x[i], r2_8, ni_8, n8);
        u32 v1 = from_dm32(mul_dm32(relaxed_sub_dm32(yi_m2, xi_m2, n2_7), m1_inv_m2, ni_7, n7), ni_7, n7);
        u32 v1_m3 = to_dm32(v1, r2_8, ni_8, n8);
        u32 v2 = from_dm32(mul_dm32(relaxed_sub_dm32(zi_m3, relaxed_add_dm32(xi_m3, mul_dm32(m1_m3, v1_m3, ni_8, n8), n8), n8), m12_inv_m3, ni_8, n8), ni_8, n8);
        u32 v2_m0 = to_dm32(v2, r2_0, ni_0, n0);
        u32 xi_m0 = to_dm32(x[i], r2_0, ni_0, n0);
        u32 v1_m0 = to_dm32(v1, r2_0, ni_0, n0);
        ret[i] = from_dm32(relaxed_add_dm32(relaxed_add_dm32(xi_m0, mul_dm32(m1_m0, v1_m0, ni_0, n0), n0), mul_dm32(m12_m0, v2_m0, ni_0, n0), n0), ni_0, n0);
    }

    free(x);
    free(y);
    free(z);

    return ret;
}

u64 *conv_mod_2_64(int a_len, int b_len, u64 *a, u64 *b) {
    u32 *A = (u32 *)calloc(a_len, sizeof(u32));
    u32 *B = (u32 *)calloc(b_len, sizeof(u32));

    if (A == NULL || B == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < a_len; i++)
        A[i] = to_dm32(a[i] % n1, r2_1, ni_1, n1);
    for (int i = 0; i < b_len; i++)
        B[i] = to_dm32(b[i] % n1, r2_1, ni_1, n1);

    u32 *C1 = conv1(a_len, b_len, A, B);

    for (int i = 0; i < a_len; i++)
        A[i] = to_dm32(a[i] % n2, r2_2, ni_2, n2);
    for (int i = 0; i < b_len; i++)
        B[i] = to_dm32(b[i] % n2, r2_2, ni_2, n2);

    u32 *C2 = conv2(a_len, b_len, A, B);

    for (int i = 0; i < a_len; i++)
        A[i] = to_dm32(a[i] % n3, r2_3, ni_3, n3);
    for (int i = 0; i < b_len; i++)
        B[i] = to_dm32(b[i] % n3, r2_3, ni_3, n3);

    u32 *C3 = conv3(a_len, b_len, A, B);

    for (int i = 0; i < a_len; i++)
        A[i] = to_dm32(a[i] % n4, r2_4, ni_4, n4);
    for (int i = 0; i < b_len; i++)
        B[i] = to_dm32(b[i] % n4, r2_4, ni_4, n4);

    u32 *C4 = conv4(a_len, b_len, A, B);

    for (int i = 0; i < a_len; i++)
        A[i] = to_dm32(a[i] % n5, r2_5, ni_5, n5);
    for (int i = 0; i < b_len; i++)
        B[i] = to_dm32(b[i] % n5, r2_5, ni_5, n5);

    u32 *C5 = conv5(a_len, b_len, A, B);

    free(A);
    free(B);

    u64 *ret = (u64 *)calloc(a_len + b_len - 1, sizeof(u64));
    if (ret == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < a_len + b_len - 1; i++) {
        u32 m2_b = mul_dm32(328554155u, relaxed_sub_dm32(C2[i], to_dm32(from_dm32(C1[i], ni_1, n1), r2_2, ni_2, n2), n2), ni_2, n2);
        u32 m3_b = mul_dm32(145185674u, relaxed_sub_dm32(C3[i], relaxed_add_dm32(mul_dm32(to_dm32(from_dm32(m2_b, ni_2, n2), r2_3, ni_3, n3), to_dm32(n1, r2_3, ni_3, n3), ni_3, n3), to_dm32(from_dm32(C1[i], ni_1, n1), r2_3, ni_3, n3), n2_3), n2_3), ni_3, n3);
        u32 m4_b = mul_dm32(228777623u, relaxed_sub_dm32(C4[i], relaxed_add_dm32(mul_dm32(relaxed_add_dm32(mul_dm32(to_dm32(from_dm32(m3_b, ni_3, n3), r2_4, ni_4, n4), to_dm32(n2, r2_4, ni_4, n4), ni_4, n4), to_dm32(from_dm32(m2_b, ni_2, n2), r2_4, ni_4, n4), n2_4), to_dm32(n1, r2_4, ni_4, n4), ni_4, n4), to_dm32(from_dm32(C1[i], ni_1, n1), r2_4, ni_4, n4), n2_4), n2_4), ni_4, n4);
        u32 m5_b = mul_dm32(578664300u, relaxed_sub_dm32(C5[i], relaxed_add_dm32(mul_dm32(relaxed_add_dm32(mul_dm32(relaxed_add_dm32(mul_dm32(to_dm32(from_dm32(m4_b, ni_4, n4), r2_5, ni_5, n5), to_dm32(n3, r2_5, ni_5, n5), ni_5, n5), to_dm32(from_dm32(m3_b, ni_3, n3), r2_5, ni_5, n5), n2_5), to_dm32(n2, r2_5, ni_5, n5), ni_5, n5), to_dm32(from_dm32(m2_b, ni_2, n2), r2_5, ni_5, n5), n2_5), to_dm32(n1, r2_5, ni_5, n5), ni_5, n5), to_dm32(from_dm32(C1[i], ni_1, n1), r2_5, ni_5, n5), n2_5), n2_5), ni_5, n5);
        u64 b1 = (u64)from_dm32(C1[i], ni_1, n1);
        u64 b2 = (u64)from_dm32(m2_b, ni_2, n2);
        u64 b3 = (u64)from_dm32(m3_b, ni_3, n3);
        u64 b4 = (u64)from_dm32(m4_b, ni_4, n4);
        u64 b5 = (u64)from_dm32(m5_b, ni_5, n5);
        ret[i] = (((b5 * n4 + b4) * n3 + b3) * n2 + b2) * n1 + b1;
    }

    free(C1);
    free(C2);
    free(C3);
    free(C4);
    free(C5);

    return ret;
}

void supset_zeta_transform(int A_len, u32 *A) {
    for (u32 w = 1; w != A_len; w <<= 1)
        for (u32 k = 0; k != A_len; k += w * 2)
            for (u32 i = 0; i != w; ++i)
                A[k + i] = add_m32(A[k + i], A[k + w + i]);
}
void supset_mobius_transform(int A_len, u32 *A) {
    for (u32 w = 1; w != A_len; w <<= 1)
        for (u32 k = 0; k != A_len; k += w * 2)
            for (u32 i = 0; i != w; i += 1)
                A[k + i] = sub_m32(A[k + i], A[k + w + i]);
}
u32 *conv_bitwise_and(int A_len, int B_len, u32 *A, u32 *B) {
    assert(A_len == B_len);
    supset_zeta_transform(A_len, A);
    supset_zeta_transform(B_len, B);
    u32 *ret = (u32 *)calloc(A_len, sizeof(u32));
    if (ret == NULL)
        exit(EXIT_FAILURE);
    for (int i = 0; i < A_len; i++)
        ret[i] = mul_m32(A[i], B[i]);
    supset_mobius_transform(A_len, ret);
    return ret;
}

void walsh_hadamard_transform(int A_len, u32 *A) {
    for (int i = 1; i < A_len; i <<= 1)
        for (int j = 0; j < A_len; j += i << 1)
            for (int k = 0; k < i; k++) {
                u32 s = A[j + k];
                u32 t = A[j + k + i];
                A[j + k] = add_m32(s, t);
                A[j + k + i] = sub_m32(s, t);
            }
}
u32 *conv_bitwise_xor(int A_len, int B_len, u32 *A, u32 *B) {
    u32 *ret = (u32 *)calloc(A_len, sizeof(u32));
    if (ret == NULL)
        exit(EXIT_FAILURE);
    walsh_hadamard_transform(A_len, A);
    walsh_hadamard_transform(B_len, B);
    for (int i = 0; i < A_len; i++)
        ret[i] = mul_m32(A[i], B[i]);
    walsh_hadamard_transform(A_len, ret);
    u32 inv2t = inv_m32(to_m32(A_len));
    for (int i = 0; i < A_len; i++)
        ret[i] = mul_m32(ret[i], inv2t);
    return ret;
}


void solve_conv_mod_998(void) {
    u32 A[524288], B[524288];
    int n, m;
    rd_int(&n);
    rd_int(&m);
    u32 x;
    for (int i = 0; i < n; i++) {
        rd_u32(&x);
        A[i] = to_dm32(x, r2_1, ni_1, n1);
    }
    for (int i = 0; i < m; i++) {
        rd_u32(&x);
        B[i] = to_dm32(x, r2_1, ni_1, n1);
    }
    u32 *C = conv1(n, m, A, B);
    for (int i = 0; i < n + m - 1; i++) {
        if (i)
            wt_char(' ');
        wt_u32(from_dm32(C[i], ni_1, n1));
    }
    wt_char('\n');
    free(C);
}
void solve_conv_mod_107(void) {
    u32 a[1 << 19], b[1 << 19];
    int n, m;
    rd_int(&n);
    rd_int(&m);
    for (int i = 0; i < n; i++)
        rd_u32(&a[i]);
    for (int i = 0; i < m; i++)
        rd_u32(&b[i]);
    u32 *c = conv_mod_arbitrary(1000000007, n, m, a, b);
    for (int i = 0; i < n + m - 1; i++) {
        if (i)
            wt_char(' ');
        wt_u32(c[i]);
    }
    wt_char('\n');
    free(c);
}
void solve_conv_mod_2_64(void) {
    u64 a[1 << 19], b[1 << 19];
    int n, m;
    rd_int(&n);
    rd_int(&m);
    for (int i = 0; i < n; i++)
        rd_u64(&a[i]);
    for (int i = 0; i < m; i++)
        rd_u64(&b[i]);
    u64 *c = conv_mod_2_64(n, m, a, b);
    for (int i = 0; i < n + m - 1; i++) {
        if (i)
            wt_char(' ');
        wt_u64(c[i]);
    }
    wt_char('\n');
    free(c);
}

void solve_conv_bitwise_and(void) {
    set_m32(998244353u);
    u32 A[1 << 20], B[1 << 20];
    int N;
    rd_int(&N);
    N = 1 << N;
    u32 x;
    for (int i = 0; i < N; i++) {
        rd_u32(&x);
        A[i] = to_m32(x);
    }
    for (int i = 0; i < N; i++) {
        rd_u32(&x);
        B[i] = to_m32(x);
    }
    u32 *C = conv_bitwise_and(N, N, A, B);
    for (int i = 0; i < N; i++) {
        if (i)
            wt_char(' ');
        wt_u32(from_m32(C[i]));
    }
    wt_char('\n');
    free(C);
}

void solve_conv_bitwise_or(void) {
    set_m32(998244353u);
    u32 A[1 << 20], B[1 << 20];
    int N;
    rd_int(&N);
    N = 1 << N;
    u32 x;
    for (int i = 0; i < N; i++) {
        rd_u32(&x);
        A[i] = to_m32(x);
    }
    for (int i = 0; i < N; i++) {
        rd_u32(&x);
        B[i] = to_m32(x);
    }
    u32 *C = conv_bitwise_xor(N, N, A, B);
    for (int i = 0; i < N; i++) {
        if (i)
            wt_char(' ');
        wt_u32(from_m32(C[i]));
    }
    wt_char('\n');
    free(C);
    return 0;
}
