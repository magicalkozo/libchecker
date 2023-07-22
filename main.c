#include "template.h"

int main(void)
{
    u32 a = pcg_range(1, 10);
    assert(a > 100);
    return 0;
}
