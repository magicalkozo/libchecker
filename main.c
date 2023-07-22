#include "template.h"

int main(void)
{
    set_m32(1000000007);
    u32 a, b;
    for (int _ = 0; _ < 100000; _++)
    {
        a = pcg_range(1, 1000000);
        b = pcg_range(1, 1000000);
        assert(add_m32(a, b) == shrink_add_m32(a, b));
    }
    return 0;
}