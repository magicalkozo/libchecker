#include "template.h"

int main(void)
{
    set_m32(998244353u);
    u32 a, b;
    for (int _ = 0; _ < 10000000000; _++)
    {
        a = pcg_range(0, 998244352u);
        b = pcg_range(0, 998244352u);
        assert(add_m32(a, b) == (a + b) % 998244353u);
    }
    return 0;
}
