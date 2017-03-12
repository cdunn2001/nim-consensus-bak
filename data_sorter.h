#include <stddef.h>

typedef int seq_coor_t; 

typedef struct {
    seq_coor_t d;
    seq_coor_t k;
    seq_coor_t pre_k;
    seq_coor_t x1;
    seq_coor_t y1;
    seq_coor_t x2;
    seq_coor_t y2;
} d_path_data2;

extern void d_path_data2sort(d_path_data2* base, size_t length, size_t max_idx);
