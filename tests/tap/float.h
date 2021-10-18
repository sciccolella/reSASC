/*
 * Floating point check function for the TAP protocol.
 *
 * This file is part of C TAP Harness.  The current version plus supporting
 * documentation is at <https://www.eyrie.org/~eagle/software/c-tap-harness/>.
 *
 * Copyright 2008, 2010, 2012, 2014, 2018 Russ Allbery <eagle@eyrie.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * SPDX-License-Identifier: MIT
 */

#ifndef TAP_FLOAT_H
#define TAP_FLOAT_H 1

#include <tests/tap/macros.h>

BEGIN_DECLS

/* Compare two values within epsilon. */
int is_double(double, double, double epsilon, const char *format, ...)
    __attribute__((__format__(printf, 4, 5)));

END_DECLS

#endif /* TAP_FLOAT_H */
