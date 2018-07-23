//
// Created by Michel Lesoinne on 7/20/18.
//

#ifndef FEM_TYPES_H
#define FEM_TYPES_H
#define USE_INT32

#ifdef USE_INT32
using gl_num_t = int32_t;
#else
using gl_num_t = int64_t;
#endif
using lc_num_t = int;

using gl_node_idx = gl_num_t;
using lc_node_idx = lc_num_t;

#endif //FEM_TYPES_H
