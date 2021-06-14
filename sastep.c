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

#include "sastep.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef NDEBUG
#include <assert.h>
#else
#define assert(ignore) ((void)0)
#endif

int random_assignment(int MAX) {
  unsigned long randval = genrand_int32();
  return (int)(randval % (MAX + 1));
}

/*void check_subtree_losses(node_t *node, vector *tree_vec, vector *loss_vec,
                          int *k_loss, int *sigma, int n) {
  if (node == NULL)
    return;

  check_subtree_losses(node->first_child, tree_vec, loss_vec, k_loss, sigma, n);
  check_subtree_losses(node->next_sibling, tree_vec, loss_vec, k_loss, sigma,
                       n);

  if (node->loss == 1) {
    bool valid = is_loss_valid(node);
    bool lost = is_already_lost(node, node->mut_index);

    if (valid == false || lost == true) {
      if(vector_get(tree_vec, node->id) != NULL){
        node_delete(node, tree_vec, loss_vec, k_loss, sigma, n);
      }
    }
  }
}

void check_subtree_recurrences(node_t *node, vector *tree_vec, vector *recs_vec,
                               int *r_recs, int *og_muts_idx, int *sigma,
                               int n) {
  if (node == NULL)
    return;
  check_subtree_recurrences(node->first_child, tree_vec, recs_vec, r_recs,
                            og_muts_idx, sigma, n);
  check_subtree_recurrences(node->next_sibling, tree_vec, recs_vec, r_recs,
                            og_muts_idx, sigma, n);

  if (node->recurrent == 1) {
    // check if valid
    bool valid = is_recurrence_valid(node);
    // check if the og is valid, if it's not delete the recurrent
    node_t *og_mut_node = vector_get(tree_vec, og_muts_idx[node->mut_index]);

    if(node->mut_index != og_mut_node->mut_index){
      exit(0);
    }

    int x = og_muts_idx[node->mut_index];
    if (x == -1){
      valid = false;
    }
    bool og_valid = is_recurrence_valid(og_mut_node);

    if (valid == false || og_valid == false)
      if(vector_get(tree_vec, node->id) != NULL)
        node_delete(node, tree_vec, recs_vec, r_recs, sigma, n);
  }
}*/

void check_subtree(node_t *node, vector *tree_vec, vector *recs_vec,
                               int *r_recs, int *og_muts_idx, int *sigma,
                               int n, vector *loss_vec, int *k_loss) {
  if (node == NULL)
    return;
  check_subtree(node->first_child, tree_vec, recs_vec, r_recs,
                            og_muts_idx, sigma, n, loss_vec, k_loss);
  check_subtree(node->next_sibling, tree_vec, recs_vec, r_recs,
                            og_muts_idx, sigma, n, loss_vec, k_loss);

  if (node->recurrent == 1) {
    // check if valid
    bool valid = is_recurrence_valid(node);
    // check if the og is valid, if it's not delete the recurrent
    node_t *og_mut_node = vector_get(tree_vec, og_muts_idx[node->mut_index]);

    if(node->mut_index != og_mut_node->mut_index){
      exit(0);
    }

    int x = og_muts_idx[node->mut_index];
    if (x == -1){
      valid = false;
    }
    bool og_valid = is_recurrence_valid(og_mut_node);

    if (valid == false || og_valid == false) {
      if(vector_get(tree_vec, node->id) != NULL){
        node_delete(node, tree_vec, recs_vec, r_recs, sigma, n);
        check_subtree(node->parent, tree_vec, recs_vec, r_recs,
                                  og_muts_idx, sigma, n, loss_vec, k_loss);
      }
    }
  } else if (node->loss == 1) {

    bool valid = is_loss_valid(node);
    bool lost = is_already_lost(node, node->mut_index);

    if (valid == false || lost == true) {
      if(vector_get(tree_vec, node->id) != NULL){
        node_delete(node, tree_vec, loss_vec, k_loss, sigma, n);
        check_subtree(node->parent, tree_vec, recs_vec, r_recs,
                                  og_muts_idx, sigma, n, loss_vec, k_loss);
      }
    }
  }
}

node_t *climb_tree(node_t *node, vector *tree_vec, bool *back, node_t *res){

  if(node == NULL){
    *back = true;
    return NULL;
  }

  node_t* par = node->parent;

  if(*back == false)
    res = climb_tree(par, tree_vec, back, res);

 if(*back == true && node->recurrent == 1){
    *back = false;
    return node;
  }

  return res;
}

void check_tree(node_t *node, vector *tree_vec, vector *recs_vec,
                               int *r_recs, int *og_muts_idx, int *sigma,
                               int n, vector *loss_vec, int *k_loss) {

  if (node == NULL)
    return;

  node_t* test;
  node_t* test2;

  bool gino = false;

  test2 = climb_tree(node, tree_vec, &gino, test);

  if(test2 != NULL){
    if(test2->recurrent == 1 && is_ancestor(node, test2))
    check_subtree(test2, tree_vec, recs_vec, r_recs, og_muts_idx,
                            sigma, n, loss_vec, k_loss);
    return;
  }
  check_subtree(node, tree_vec, recs_vec, r_recs, og_muts_idx,
                          sigma, n, loss_vec, k_loss);
}

int prune_regraft(node_t *prune, node_t *regraft, node_t *root,
                  int monoclonal) {
  if (is_ancestor(regraft, prune) == true)
    return 1;
  if (regraft->parent == prune)
    return 1;
  if (prune->parent == NULL)
    return 1;
  if (prune == regraft)
    return 1;

  if (monoclonal == 1 && regraft->parent == NULL)
    return 1;

  node_detach(prune);
  node_append(regraft, prune);

  return 0;
}

int add_back_mutation(node_t *node, vector *tree_vec, int m, int k, int *k_loss,
                      vector *losses_vec, int MAX_LOSSES) {
  node_t *par = node->parent;
  if (par == NULL || par->parent == NULL)
    return 1;

  if (vector_total(losses_vec) >= MAX_LOSSES)
    return 1;

  // Walk to root to select possible candidates for deletion
  node_t *candidates[m * k + 1];
  for (int i = 0; i < k * m + 1; i++) {
    candidates[i] = NULL;
  }

  int tot = 0;
  while (par != NULL) {
    if (par->loss == 0) {
      candidates[tot] = par;
      tot++;
    }
    par = par->parent;
  }

  int rand_del = random_assignment(tot - 2) + 1;
  if (candidates[rand_del] == NULL)
    return 1;
  if (candidates[rand_del]->mut_index == -1)
    return 1;

  char label[255];
  strcpy(label, candidates[rand_del]->label);
  if (k_loss[candidates[rand_del]->mut_index] >= k)
    return 1;
  if (is_already_lost(node, candidates[rand_del]->mut_index) == true)
    return 1;

  node_t *node_del =
      node_new(label, candidates[rand_del]->mut_index, vector_total(tree_vec));
  node_del->loss = 1;

  vector_add(tree_vec, node_del);
  vector_add(losses_vec, node_del);

  k_loss[node_del->mut_index] += 1;

  par = node->parent;
  node_detach(node);
  node_append(par, node_del);
  node_append(node_del, node);

  return 0;
}

int add_recurrent_mutation(node_t *node, vector *tree_vec, int m, int r,
                           int *r_rec, vector *recurrents_vec,
                           int MAX_RECURRENCES, int *original_mut_idx) {
  node_t *par = node->parent;
  if (par == NULL || par->parent == NULL)
    return 1;
  if (vector_total(recurrents_vec) >= MAX_RECURRENCES)
    return 1;

  int genotype_candidate[m];
  for (int i = 0; i < m; i++) {
    genotype_candidate[i] = 0;
  }
  // Walk to root to exclude nodes
  while (par != NULL) {
    if (par->mut_index != -1) {
      if (par->loss == 1)
        genotype_candidate[par->mut_index] -= 1;
      else
        genotype_candidate[par->mut_index] += 1;
    }

    par = par->parent;
  }

  int rand_rec = 0;
  do {
    rand_rec = random_assignment(m - 1);
  } while (genotype_candidate[rand_rec] != 0);

  // Mutation already reccurred too many times
  if (r_rec[rand_rec] >= r)
    return 1;

  // Current node is ancestor of the to-be-recurred mutation
  int node_id = original_mut_idx[rand_rec];
  node_t *candidate_node = vector_get(tree_vec, node_id);

  if (is_ancestor(candidate_node, node))
    return 1;
  // and is not ancestor of any of its recurrences (if any)
  int tot_rec = vector_total(recurrents_vec);
  for (int i = 0; i < tot_rec; i++) {
    node_t *n = vector_get(recurrents_vec, i);
    if (n->recurrent == 1 && n->mut_index == rand_rec && is_ancestor(n, node))
      return 1;
  }

  // The recurrence should be valid
  char label[255];
  strcpy(label, candidate_node->label);
  node_t *node_rec = node_new(label, rand_rec, vector_total(tree_vec));
  node_rec->recurrent = 1;
  r_rec[node_rec->mut_index] += 1;

  vector_add(tree_vec, node_rec);
  vector_add(recurrents_vec, node_rec);

  par = node->parent;
  node_detach(node);
  node_append(par, node_rec);
  node_append(node_rec, node);

  return 0;
}

double greedy_tree_loglikelihood(node_t *root, vector tree_vec, int *sigma,
                                 int **inmatrix, int n, int m, double *alpha,
                                 double beta, double *gammas, double* deltas,
                                 int *k_loss, int *k_recurrent, int CORES) {
  int node_max = vector_total(&tree_vec);

  int *nodes_genotypes = calloc(node_max * m, sizeof(int));

  for (int i = 0; i < node_max; i++) {
    node_t *node = vector_get(&tree_vec, i);

    if (node == NULL) {
      for (int j = 0; j < m; j++) {
        nodes_genotypes[i * m + j] = 3;
      }
    } else {
      get_genotype_profile(node, nodes_genotypes + i * m);
    }
  }

  double loss_weight = 0;
  for (int j = 0; j < m; j++) {
    loss_weight += k_loss[j] * log(gammas[j]);
  }

  double maximum_likelihood = 0 - loss_weight;

  double double_weigth = 0;
  for (int j = 0; j < m; j++) {
    double_weigth += (k_recurrent[j]) * log(deltas[j]);
  }

  maximum_likelihood = maximum_likelihood - double_weigth;

  double like_00 = log(1 - beta);
  double like_10 = log(beta);
  double *like_01 = malloc(m * sizeof(double));
  double *like_11 = malloc(m * sizeof(double));
  double like_2 = 0;

  for (int j = 0; j < m; j++) {
    like_01[j] = log(alpha[j]);
    like_11[j] = log(1 - alpha[j]);
  }

#pragma omp parallel for reduction(+ : maximum_likelihood) num_threads(CORES)
  for (int i = 0; i < n; i++) {
    int best_sigma = -1;
    double best_lh = -DBL_MAX;

    for (int node = 0; node < node_max; node++) {
      if (nodes_genotypes[node * m + 0] != 3) {
        double lh = 0;

        for (int j = 0; j < m; j++) {
          double p = 1;
          assert(nodes_genotypes[node * m + j] == 0 ||
                 nodes_genotypes[node * m + j] == 1);
          if (inmatrix[i][j] == 0 && nodes_genotypes[node * m + j] == 0) {
            p = like_00;
          } else if (inmatrix[i][j] == 0 &&
                     nodes_genotypes[node * m + j] == 1) {
            p = like_01[j];
          } else if (inmatrix[i][j] == 1 &&
                     nodes_genotypes[node * m + j] == 0) {
            p = like_10;
          } else if (inmatrix[i][j] == 1 &&
                     nodes_genotypes[node * m + j] == 1) {
            p = like_11[j];
          } else if (inmatrix[i][j] == 2) {
            p = like_2;
          }
          lh += p;
        }

        if (lh > best_lh) {
          best_sigma = node;
          best_lh = lh;
        }
      }
    }
    sigma[i] = best_sigma;
    maximum_likelihood += best_lh;
  }

  free(nodes_genotypes);
  free(like_01);
  free(like_11);
  return maximum_likelihood;
}

double gen_gaussian(double mu, double sigma) {
  static const double two_pi = 2.0 * 3.14159265358979323846;
  double u1 = genrand_real3();
  double u2 = genrand_real3();
  double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  return fabs(z0 * sigma + mu);
}

void el_commit(elpar_t *el_params, double *beta) {
  for (int i = 0; i < el_params->M; i++) {
    el_params->ALPHAS[i] = el_params->a_xs[i];
    el_params->GAMMAS[i] = el_params->g_xs[i];
    el_params->DELTAS[i] = el_params->d_xs[i];
  }
  *beta = el_params->b_x;
  el_params->changed = 0;
}

void el_discard(elpar_t *el_params, double beta) {
  for (int i = 0; i < el_params->M; i++) {
    el_params->a_xs[i] = el_params->ALPHAS[i];
    el_params->g_xs[i] = el_params->GAMMAS[i];
    el_params->d_xs[i] = el_params->DELTAS[i];
  }
  el_params->b_x = beta;
  el_params->changed = 0;
}

void learn_a(elpar_t *el_params) {
  double x = 0;
  if (el_params->single_alpha == 1) {
    x = gen_gaussian(el_params->a_mu[0], el_params->a_variance);
    for (int i = 0; i < el_params->M; i++) {
      el_params->a_xs[i] = x;
    }
  } else {
    for (int i = 0; i < el_params->M; i++) {
      x = gen_gaussian(el_params->a_mu[i], el_params->a_variance);
      el_params->a_xs[i] = x;
    }
  }
}

void learn_b(elpar_t *el_params) {
  double x = 0;
  x = gen_gaussian(el_params->b_mu, el_params->b_variance);
  el_params->b_x = x;
}

void learn_g(elpar_t *el_params) {
  double x = 0;
  if (el_params->single_gamma == 1) {
    x = gen_gaussian(el_params->g_mu[0], el_params->g_variance);
    for (int i = 0; i < el_params->M; i++) {
      el_params->g_xs[i] = x;
    }
  } else {
    for (int i = 0; i < el_params->M; i++) {
      x = gen_gaussian(el_params->g_mu[i], el_params->g_variance);
      el_params->g_xs[i] = x;
    }
  }
}

void learn_d(elpar_t *el_params) {
  double x = 0;
  if (el_params->single_delta == 1) {
    x = gen_gaussian(el_params->d_mu[0], el_params->d_variance);
    for (int i = 0; i < el_params->M; i++) {
      el_params->d_xs[i] = x;
    }
  } else {
    for (int i = 0; i < el_params->M; i++) {
      x = gen_gaussian(el_params->d_mu[i], el_params->d_variance);
      el_params->d_xs[i] = x;
    }
  }
}

void params_learning(elpar_t *el_params) {
  double choice = genrand_real1();

  if (el_params->a_variance > 0 && el_params->b_variance > 0 &&
      el_params->g_variance > 0 && el_params->d_variance > 0) {
    if (choice < 0.25) {
      learn_a(el_params);
    } else if (choice < 0.50) {
      learn_b(el_params);
    } else if (choice < 0.75){
      learn_g(el_params);
    } else {
      learn_d(el_params);
    }
  } else if(el_params->a_variance > 0 && el_params->b_variance > 0 &&
       el_params->g_variance > 0){
      if(choice < 0.33) {
        learn_a(el_params);
      } else if (choice < 0.66) {
        learn_b(el_params);
      } else {
        learn_g(el_params);
      }
  } else if(el_params->a_variance > 0 && el_params->b_variance > 0 &&
       el_params->d_variance > 0){
      if(choice < 0.33) {
        learn_a(el_params);
      } else if (choice < 0.66) {
        learn_b(el_params);
      } else {
        learn_d(el_params);
      }
  } else if(el_params->a_variance > 0 && el_params->g_variance > 0 &&
       el_params->d_variance > 0){
      if(choice < 0.33) {
        learn_a(el_params);
      } else if (choice < 0.66) {
        learn_g(el_params);
      } else {
        learn_d(el_params);
      }
  } else if(el_params->g_variance > 0 && el_params->b_variance > 0 &&
       el_params->d_variance > 0){
      if(choice < 0.33) {
        learn_g(el_params);
      } else if (choice < 0.66) {
        learn_b(el_params);
      } else {
        learn_d(el_params);
      }
  } else if (el_params->a_variance > 0 && el_params->b_variance > 0) {
    if (choice < 0.5) {
      learn_a(el_params);
    } else {
      learn_b(el_params);
    }
  } else if (el_params->a_variance > 0 && el_params->g_variance > 0) {
    if (choice < 0.5) {
      learn_a(el_params);
    } else {
      learn_g(el_params);
    }
  } else if (el_params->a_variance > 0 && el_params->d_variance > 0) {
    if (choice < 0.5) {
      learn_a(el_params);
    } else {
      learn_d(el_params);
    }
  } else if (el_params->b_variance > 0 && el_params->g_variance > 0) {
    if (choice < 0.5) {
      learn_b(el_params);
    } else {
      learn_g(el_params);
    }
  } else if (el_params->b_variance > 0 && el_params->d_variance > 0) {
    if (choice < 0.5) {
      learn_a(el_params);
    } else {
      learn_d(el_params);
    }
  } else if (el_params->d_variance > 0 && el_params->g_variance > 0) {
    if (choice < 0.5) {
      learn_d(el_params);
    } else {
      learn_g(el_params);
    }
  } else if (el_params->a_variance > 0) {
    learn_a(el_params);
  } else if (el_params->b_variance > 0) {
    learn_b(el_params);
  } else if (el_params->g_variance > 0) {
    learn_g(el_params);
  } else if (el_params->d_variance > 0) {
    learn_d(el_params);
  }

  el_params->changed = 1;
}

void rec_mut_idx(node_t *node, int *original_muts_idx){
  node_t* child = node->first_child;
  if (child == NULL)
    return;

  while(child != NULL){
    if (original_muts_idx != NULL && child->recurrent == 0 && child->loss == 0) {
      assert(original_muts_idx[child->mut_index] == -1);

      original_muts_idx[child->mut_index] = child->id;
    }

    rec_mut_idx(child, original_muts_idx);

    child = child->next_sibling;
  }
}

void neighbor(node_t *root, vector *tree_vec, int *sigma, int m, int n, int k,
              int r, vector *loss_vec, vector *rec_vec, int *k_loss,
              int *r_recs, int MAX_LOSSES, int MAX_RECURRENCES,
              int *original_muts_idx, elpar_t *el_params, int monoclonal) {
  double el = genrand_real1();
  if (el < 0.1 && (el_params->a_variance > 0 || el_params->b_variance > 0 ||
                   el_params->g_variance > 0 || el_params->d_variance > 0)) {
    params_learning(el_params);
  } else {
    double move = genrand_real1();
    if (move < 0.25 && (k > 0 || r > 0)) {
      int choice = -1;
      if (k > 0 && r == 0) {
        // add back-mutation
        choice = 0;
      } else if (r > 0 && k == 0) {
        // add recurrent-mutation
        choice = 1;
      } else {
        // TODO: change probabilities?
        choice = (genrand_real1() < 0.5) ? 0 : 1;
      }
      assert(choice != -1);

      if (choice == 0) {
        // Add back-mutation
        int bm_res = 1;
        node_t *node_res = NULL;

        int ip = random_assignment(vector_total(tree_vec) - 1);
        node_res = vector_get(tree_vec, ip);

        bm_res = add_back_mutation(node_res, tree_vec, m, k, k_loss, loss_vec,
                                   MAX_LOSSES);

        if (bm_res == 0) {
          check_tree(node_res, tree_vec, rec_vec, r_recs,
                                    original_muts_idx, sigma, n, loss_vec, k_loss);
        }
      } else {
        // Add recurrent-mutation
        int rec_res = 1;
        node_t *node_res = NULL;

        int ip = random_assignment(vector_total(tree_vec) - 1);
        node_res = vector_get(tree_vec, ip);

        rec_res = add_recurrent_mutation(node_res, tree_vec, m, r, r_recs, rec_vec,
                                          MAX_RECURRENCES, original_muts_idx);
        // // This shouldn't be necessary, since recurrences should be added
        // only if correct
        if (rec_res == 0) {
        //  print_tree(root, 0.0);
          check_tree(node_res, tree_vec, rec_vec, r_recs,
                                      original_muts_idx, sigma, n, loss_vec, k_loss);
        }
      }
    } else if (move < 0.50 && (k > 0 || r > 0)) {
      int choice = -1;
      if (k > 0 && r == 0) {
        // delete back-mutation
        choice = 0;
      } else if (r > 0 && k == 0) {
        // delete recurrent-mutation
        choice = 1;
      } else {
        // TODO: change probabilities?
        choice = (genrand_real1() < 0.5) ? 0 : 1;
      }
      assert(choice != -1);
      if (choice == 0) {
        // Delete a backmutation
        node_t *node_res = NULL;
        if (vector_total(loss_vec) == 0)
          return;

        int node_max = vector_total(loss_vec) - 1;
        assert(node_max >= 0);
        int ip = random_assignment(node_max);
        node_res = vector_get(loss_vec, ip);

        node_delete(node_res, tree_vec, loss_vec, k_loss, sigma, n);

        check_tree(node_res->parent, tree_vec, rec_vec, r_recs, original_muts_idx,
                        sigma, n, loss_vec, k_loss);
      } else {
        // Delete a recurrent mutation
        node_t *node_res = NULL;
        if (vector_total(rec_vec) == 0)
          return;

        int node_max = vector_total(rec_vec) - 1;
        assert(node_max >= 0);
        int ip = random_assignment(node_max);
        node_res = vector_get(rec_vec, ip);

        node_delete(node_res, tree_vec, rec_vec, r_recs, sigma, n);

        check_tree(node_res->parent, tree_vec, rec_vec, r_recs, original_muts_idx,
                        sigma, n, loss_vec, k_loss);
      }
    } else if ((move < 0.75 && (k > 0 || r > 0)) ||
               (move < 0.50 && (k == 0 && r == 0))) {
      // switch nodes
      node_t *u = NULL;
      while (u == NULL || u->parent == NULL || u->loss == 1 ||
            u->recurrent == 1) {
        int node_max = vector_total(tree_vec) - 1;
        assert(node_max > 0);
        int ip = random_assignment(node_max);
        u = vector_get(tree_vec, ip);
      }

      node_t *v = NULL;
      while (v == NULL || v->parent == NULL || v->loss == 1 ||
            v->recurrent == 1 || v->id == u->id) {
        int node_max = vector_total(tree_vec) - 1;
        assert(node_max > 0);
        int ig = random_assignment(node_max);
        v = vector_get(tree_vec, ig);
      }

      int mut_tmp;
      int ori_mut_tmp;
      char label_tmp[255];

      mut_tmp = u->mut_index;
      strcpy(label_tmp, u->label);

      u->mut_index = v->mut_index;
      strcpy(u->label, v->label);

      v->mut_index = mut_tmp;
      strcpy(v->label, label_tmp);

      ori_mut_tmp = original_muts_idx[u->mut_index];
      original_muts_idx[u->mut_index] = original_muts_idx[v->mut_index];
      original_muts_idx[v->mut_index] = ori_mut_tmp;

      check_tree(u, tree_vec, rec_vec, r_recs, original_muts_idx,
                      sigma, n, loss_vec, k_loss);
      check_tree(v, tree_vec, rec_vec, r_recs, original_muts_idx,
                      sigma, n, loss_vec, k_loss);


    } else {
      // Prune-regraft two random nodes
      int pr_res = 1;
      node_t *prune_res = NULL;
      while (pr_res != 0) {
        node_t *prune = NULL;
        while (prune == NULL || prune->parent == NULL) {
          int node_max = vector_total(tree_vec) - 1;
          assert(node_max > 0);
          int ip = random_assignment(node_max);
          prune = vector_get(tree_vec, ip);
        }

        node_t *graft = NULL;
        while (graft == NULL) {
          int node_max = vector_total(tree_vec) - 1;
          assert(node_max > 0);
          int ig = random_assignment(node_max);
          graft = vector_get(tree_vec, ig);
        }
        pr_res = prune_regraft(prune, graft, root, monoclonal);
        prune_res = prune;
      }

      for (int i = 0; i < m; i++) {
        original_muts_idx[i] = -1;
      }

      node_t *child = root->first_child;

      while(child != NULL){

        if (original_muts_idx != NULL && child->recurrent == 0 &&
           child->loss == 0) {
          assert(original_mut_idx[child->mut_index] == -1);

          original_muts_idx[child->mut_index] = child->id;
        }

        rec_mut_idx(child, original_muts_idx);

        child = child->next_sibling;
      }
      check_tree(prune_res, tree_vec, rec_vec, r_recs,
                      original_muts_idx, sigma, n, loss_vec, k_loss);

    }
  }
}

elpar_t *set_el_params(int single_a, int m, double *ALPHAS, double *a_mu,
                       double a_variance, double *a_xs, double *beta,
                       double b_mu, double b_variance, double *GAMMAS,
                       double *g_mu, double g_variance, double *g_xs,
                       int single_g, double* DELTAS, double* d_mu,
                       double d_variance, double* d_xs, int single_d) {
  elpar_t *params = malloc(sizeof(elpar_t));
  params->single_alpha = single_a;
  params->M = m;
  params->changed = 0;

  params->ALPHAS = ALPHAS;
  params->a_mu = a_mu;
  params->a_variance = a_variance;
  params->a_xs = a_xs;

  params->BETA = beta;
  params->b_mu = b_mu;
  params->b_variance = b_variance;
  params->b_x = b_mu;

  params->single_gamma = single_g;

  params->GAMMAS = GAMMAS;
  params->g_mu = g_mu;
  params->g_variance = g_variance;
  params->g_xs = g_xs;

  params->DELTAS = DELTAS;
  params->d_mu = d_mu;
  params->d_variance = d_variance;
  params->d_xs = d_xs;

  params->single_delta = single_d;

  return params;
}

double accept_prob(double old_lh, double new_lh, double current_temp) {
  if (new_lh > old_lh) {
    return 1.0;
  } else {
    double a = exp((new_lh - old_lh) / current_temp);
    return a;
  }
}

node_t *anneal(node_t *root, vector tree_vec, int n, int m, int k, int r,
               double *alpha, double beta,  double* delta, int* Fj,
               int **inmatrix, double start_temp, double cooling_rate,
               double min_temp, int MAX_LOSSES, int MAX_RECURRENCES,
               elpar_t *el_params, double *gamma, int *Cj, int MONOCLONAL,
               int CORES) {
  double current_temp = start_temp;
  double current_cooling_rate = cooling_rate;

  // Vector of node indices
  vector current_tree_vec;
  vector_init(&current_tree_vec);

  // Vector of back-mutation node indices
  vector current_losses_vec;
  vector_init(&current_losses_vec);

  // Vector of recurrent-mutations node indices
  vector current_recs_vec;
  vector_init(&current_recs_vec);

  // Current count of losses per mutation to ensure them to be < k
  int current_kloss[m];
  for (int i = 0; i < m; i++) {
    current_kloss[i] = 0;
  }

  // Current assignment
  int current_sigma[n];
  for (int i = 0; i < n; i++) {
    current_sigma[i] = 0;
  }

  // array of mutation_index->node_id only for the original M mutations
  int current_og_mut_idx[m];
  for (int i = 0; i < m; i++) {
    current_og_mut_idx[i] = -1;
  }
  int current_rrecs[m];
  for (int i = 0; i < m; i++) {
    current_rrecs[i] = 0;
  }

  node_t *current_root = treecpy(root, &current_tree_vec, &current_losses_vec,
                                 &current_recs_vec, n, current_og_mut_idx);

  unsigned int step = 0;

  double current_lh =
      greedy_tree_loglikelihood(current_root, tree_vec, current_sigma, inmatrix,
                                n, m, alpha, beta, gamma, delta, current_kloss,
                                current_rrecs, CORES);

  printf("Step\t\t\tLog-likelihood\t\t\tTemperature\n");
  while (current_temp > min_temp) {
    // Create a modifiable copy
    vector copy_tree_vec;
    vector_init(&copy_tree_vec);

    vector copy_losses_vec;
    vector_init(&copy_losses_vec);

    vector copy_recs_vec;
    vector_init(&copy_recs_vec);

    int copy_kloss[m];
    for (int i = 0; i < m; i++) {
      copy_kloss[i] = current_kloss[i];
    }

    int copy_sigma[n];
    for (int i = 0; i < n; i++) {
      copy_sigma[i] = current_sigma[i];
    }

    int copy_og_mut_idx[m];
    for (int i = 0; i < m; i++) {
      copy_og_mut_idx[i] = -1;
    }
    int copy_rrces[m];
    for (int i = 0; i < m; i++) {
      copy_rrces[i] = current_rrecs[i];
    }

    node_t *copy_root = treecpy(current_root, &copy_tree_vec, &copy_losses_vec,
                                &copy_recs_vec, n, copy_og_mut_idx);

    for (int i = 0; i < vector_total(&copy_tree_vec); i++) {
      node_t *n = vector_get(&copy_tree_vec, i);
      assert(n->id == i);
    }
    assert(vector_total(&copy_tree_vec) == vector_total(&current_tree_vec));
    assert(vector_total(&copy_losses_vec) == vector_total(&current_losses_vec));
    assert(vector_total(&copy_recs_vec) == vector_total(&current_recs_vec));
    print_tree(copy_root, 10);

    neighbor(copy_root, &copy_tree_vec, copy_sigma, m, n, k, r,
             &copy_losses_vec, &copy_recs_vec, copy_kloss, copy_rrces,
             MAX_LOSSES, MAX_RECURRENCES, copy_og_mut_idx, el_params,
             MONOCLONAL);

    double new_lh = 0;
    if (el_params->changed == 1) {
      new_lh = greedy_tree_loglikelihood(
          copy_root, copy_tree_vec, copy_sigma, inmatrix, n, m, el_params->a_xs,
          el_params->b_x, el_params->g_xs, el_params->d_xs, copy_kloss, copy_rrces, CORES);
    } else {
      new_lh = greedy_tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma,
                                         inmatrix, n, m, alpha, beta, gamma,
                                         delta, copy_kloss, copy_rrces, CORES);
    }

    double acceptance = accept_prob(current_lh, new_lh, current_temp);
    double p = genrand_real1();

    if (acceptance > p) {
      if (el_params->changed == 1) {
        el_commit(el_params, &beta);
      }
      // Destroy current solution
      destroy_tree(current_root);

      vector_free(&current_tree_vec);
      vector_init(&current_tree_vec);

      vector_free(&current_losses_vec);
      vector_init(&current_losses_vec);

      vector_free(&current_recs_vec);
      vector_init(&current_recs_vec);

      // Copy new solution to current

      assert(greedy_tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma,
                                       inmatrix, n, m, alpha, beta, gamma,
                                       delta, copy_kloss, copy_rrces, CORES) == new_lh);

      current_lh = new_lh;

      for (int i = 0; i < m; i++) {
        current_kloss[i] = copy_kloss[i];
        current_rrecs[i] = copy_rrces[i];
        current_og_mut_idx[i] = 0;
      }

      current_root = treecpy(copy_root, &current_tree_vec, &current_losses_vec,
                             &current_recs_vec, n, current_og_mut_idx);

      assert(greedy_tree_loglikelihood(
                 current_root, current_tree_vec, current_sigma, inmatrix, n, m,
                 alpha, beta, gamma, delta, current_kloss, current_rrecs, CORES) == test_lh);
    } else {
      if (el_params->changed == 1) {
        el_discard(el_params, beta);
      }
    }

    // Destroy neighbor
    destroy_tree(copy_root);
    vector_free(&copy_tree_vec);
    vector_free(&copy_losses_vec);
    vector_free(&copy_recs_vec);

    if (current_temp <= 0.25 * start_temp) {
      current_temp *= (1 - (current_cooling_rate * 0.1));
    } else {
      current_temp *= (1 - current_cooling_rate);
    }

    ++step;
    if (step % 1000 == 0 || step == 1) {
      printf("%d\t\t\t%lf\t\t\t%lf\n", step, current_lh, current_temp);
    }
  }

  // Change SIGMA to best solution
  int bs[n];
  for (int i = 0; i < n; i++) {
    bs[i] = 0;
  }

  vector_init(&tree_vec);
  node_t *best_root = treecpy(current_root, &tree_vec, NULL, NULL, n, NULL);

  for (int i = 0; i < m; i++) {
    Cj[i] = current_kloss[i];
  }

  for (int i = 0; i < m; i++) {
    Fj[i] = current_rrecs[i];
  }

  // Return best TREE
  return best_root;
}
