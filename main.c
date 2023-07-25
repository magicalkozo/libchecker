#include "template.h"

int main(void)
{
    u32 mod = 998244353u;
    set_m32(mod);
    for (int _ = 0; _ < 1000; _++)
    {
        u32 a = pcg_range(0u, 998244352u);
        u32 b = pcg_range(0u, 998244352u);
        u32 B = to_m32(b);
        assert(mul_m32(a, b) == mul_m32(a, B));
    }
    return 0;
}
