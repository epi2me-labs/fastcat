// sdust.h --- taken from minimap2 source code
//
// The MIT License
// 
// Copyright (c) 2018-     Dana-Farber Cancer Institute
//               2017-2018 Broad Institute, Inc.
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#ifndef SDUST_H
#define SDUST_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct sdust_buf_s;
typedef struct sdust_buf_s sdust_buf_t;

// the simple interface
uint64_t *sdust(void *km, const uint8_t *seq, int l_seq, int T, int W, int *n);

// the following interface dramatically reduce heap allocations when sdust is frequently called.
sdust_buf_t *sdust_buf_init(void *km);
void sdust_buf_destroy(sdust_buf_t *buf);
const uint64_t *sdust_core(const uint8_t *seq, int l_seq, int T, int W, int *n, sdust_buf_t *buf);

#ifdef __cplusplus
}
#endif

#endif
