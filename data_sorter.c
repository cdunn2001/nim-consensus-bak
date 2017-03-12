#include "data_sorter.h"
#include <assert.h>
#include <stdlib.h>

int compare_path_data2(const void* a, const void* b)
{
    const d_path_data2* arg1 = a;
    const d_path_data2* arg2 = b;
    if (arg1->d - arg2->d == 0) {
        return  arg1->k - arg2->k;
    } else {
        return arg1->d - arg2->d;
    }
}


void d_path_data2sort(d_path_data2* base, size_t max_idx, size_t size) {
    assert(size == sizeof(d_path_data2));
    qsort(base, max_idx, size, compare_path_data2);
}
