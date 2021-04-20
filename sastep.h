/*
    MIT License

    Copyright (c) 2017-2019 Simone Ciccolella

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to
   deal in the Software without restriction, including without limitation the
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
   IN THE SOFTWARE.
*/

#ifndef SASTEP_H
#define SASTEP_H

#include "mt19937ar.h"
#include "tree.h"
#include "vector.h"

typedef struct el_params {
  int single_alpha;
  int single_gamma;
  int single_delta;
  int M;

  double *ALPHAS;
  double *a_mu;
  double a_variance;
  double *a_xs;

  double *BETA;
  double b_mu;
  double b_variance;
  double b_x;

  double *GAMMAS;
  double *g_mu;
  double g_variance;
  double *g_xs;

  double* DELTAS;
  double* d_mu;
  double d_variance;
  double* d_xs;

  int changed;
} elpar_t;

int random_assignment(int MAX);

double greedy_tree_loglikelihood(node_t *root, vector tree_vec, int *sigma,
                                 int **inmatrix, int n, int m, double *alpha,
                                 double beta, double *gammas, double* deltas,
                                 int *k_loss, int *k_recurrent, int CORES);

node_t *anneal(node_t *root, vector tree_vec, int n, int m, int k, int r,
               double *alpha, double beta, double* delta, int* Fj, int **inmatrix, double start_temp,
               double cooling_rate, double min_temp, int MAX_LOSSES,
               int MAX_RECURRENCES, elpar_t *el_params, double *gamma, int *Cj,
               int MONOCLONAL, int CORES);

elpar_t *set_el_params(int single, int m, double *ALPHAS, double *a_mu,
                       double a_variance, double *a_xs, double *beta,
                       double b_mu, double b_variance, double *GAMMAS,
                       double *g_mu, double g_variance, double *g_xs,
                       int single_g, double* DELTAS, double* d_mu,
                       double d_variance, double* d_xs);

#endif
