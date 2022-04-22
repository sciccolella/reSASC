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
  bool test_tree;
  bool test_add_backmutation;
  bool test_add_recurrent;
} elpar_t;

/**
* @brief Function that creates a random 32 bit int with MAX as his maximun limit
*
* @param MAX int that rapresent the maximun limit
* @return int The random number that has been created
*/
int random_assignment(int MAX);

/**
* @brief Function that calculate the likelihood of a tree
*
* @param root Root of the tree
* @param tree_vec Vector of *node_t that rapresent the tree
* @param sigma Array of int
* @param inmatrix Single cell matrix of input
* @param n Number of cells in the input file
* @param m Number of mutations in the input file
* @param alpha False negative rate in the input file for each cell
* @param beta False positive rate in the input file
* @param gammas Loss rate in the input file for each cell
* @param deltas Recurrence rate in the input file for each cell
* @param k_loss Array of int that contains the numbers of deletions of each mutation
* @param k_recurrent Array of int that contains the numbers of recurrences of each muation
* @param CORES Total number of cores that can be used
* @return double Likelihood of the tree
*/
double greedy_tree_loglikelihood(node_t *root, vector tree_vec, int *sigma,
                                 int **inmatrix, int n, int m, double *alpha,
                                 double beta, double *gammas, double* deltas,
                                 int *k_loss, int *k_recurrent, int CORES);

/**
* @brief Function that search the max likelihood by making changes to the tree with
* the Simelated annealing algorithm
*
* @param root Root of the tree
* @param tree_vec Vector of *node_t that rapresent the tree
* @param n Number of cells in the input file
* @param m Number of mutations in the input file
* @param k Max number of losses for each mutation
* @param r Max number of recurrences for each mutation
* @param alpha False negative rate in the input file for each cell
* @param beta False positive rate in the input file
* @param delta Recurrence rate in the input file for each cell
* @param Fj Array that contains all the current numbers of recurrences of each mutations
* @param inmatrix Single cell matrix of input
* @param start_temp Starting temperature of the Simulated annealing algorithm
* @param cooling_rate Cooling rate of the Simelated annealing algorithm
* @param min_temp Temperature at which the Simelated annealing algorithm stops
* @param MAX_LOSSES Max number of deletions allowed in the solution
* @param MAX_RECURRENCES Max number of recurrences allowed in the solution
* @param el_params struct that contains all the parameters that are needed
* @param gamma Loss rate in the input file for each cell
* @param Cj Array that contains all the current numbers of losses of each mutations
* @param MONOCLONAL Int that is 1 if the tree is monoclonal, 0 if not
* @param CORES Total number of cores that can be used
* @return node_t* Return the root of the tree with the best likelihood found
*/
node_t *anneal(node_t *root, vector tree_vec, int n, int m, int k, int r,
               double *alpha, double beta, double* delta, int* Fj, int **inmatrix,
               double start_temp, double cooling_rate, double min_temp, int MAX_LOSSES,
               int MAX_RECURRENCES, elpar_t *el_params, double *gamma, int *Cj,
               int MONOCLONAL, int CORES);

/**
* @brief Function that creates el_params, struct that contains all the paramaters
* and set them to the values passsed
*
* @param single_a Int that is 1 if there is a single alpha, 0 if there are different
* ones for each mutation
* containing all the alphas
* @param m Number of mutation in the input file
* @param ALPHAS Array containing the current false negative probability for each mutation
* @param a_mu Array containing false negative rate for each mutation
* @param a_variance Standard variation for new false negative discovery
* @param a_xs Array containing the new false negative probability for each mutation
* @param beta Array containing the false positive probability for each mutation
* @param b_mu False positive rate
* @param b_variance Standard variation for new false positive discovery
* @param GAMMAS Array containing the current loss probability for each mutation
* @param g_mu Array containing loss rate for each mutation
* @param g_variance Standard variation for new gamma discovery
* @param g_xs Array containing the new loss probability for each mutation
* @param single_g Int that is 1 if there is a single gamma, 0 if there are different
* ones for each mutation
* @param DELTAS Array containing the current recurrent probability for each mutation
* @param d_mu Array containing recurrent rate for each mutation
* @param d_variance Standard variation for new delta discovery
* @param d_xs Array containing the new recurrent probability for each mutation
* @param single_d Int that is 1 if there is a single delta, 0 if there are different
* ones for each mutation
* @return elpar_t* Struct that contains all the paramaters
*/
elpar_t *set_el_params(int single_a, int m, double *ALPHAS, double *a_mu,
                       double a_variance, double *a_xs, double *beta,
                       double b_mu, double b_variance, double *GAMMAS,
                       double *g_mu, double g_variance, double *g_xs,
                       int single_g, double* DELTAS, double* d_mu,
                       double d_variance, double* d_xs, int single_d);

#endif
