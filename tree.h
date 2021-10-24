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

#ifndef TREE_H
#define TREE_H

#include "vector.h"
#include <stdbool.h>

typedef struct Node {
  int id;
  char label[255];
  int loss;
  int recurrent;
  int mut_index;

  struct Node *first_child;
  struct Node *next_sibling;
  struct Node *previous_sibling;
  struct Node *parent;
} node_t;

/**
 * @brief Create a new node.
 *
 * @param label Name of the node
 * @param mut_index Index of mutation encoded by the node
 * @param id Node's ID, this should be used with a global incrementing counter
 * @param tree Array of tree[ID] = *node_t, used to access nodes directly
 * @return node_t* Return pointer to new allocated node
 */
node_t *node_new(char *label, int mut_index, int id);

/**
 * @brief Append node as child of node
 *
 * @param parent Pointer to parent
 * @param node Pointer to node
 */
void node_append(node_t *parent, node_t *node);

/**
 * @brief Deallocate a tree
 *
 * @param node Root of the (sub)tree
 */
void destroy_tree(node_t *node);

/**
 * @brief Print tree in DOT code on stdout
 *
 * @param root Root of the (sub)tree
 */
void print_tree(node_t *root, double score);

/**
 * @brief Check if cand_ancestor is an ancestor of node
 *
 * @param node Current node
 * @param cand_ancestor Candidate ancestor
 * @return true cand_ancestor is ancestor of node
 * @return false cand_ancestor is NOT ancestor of node
 */
bool is_ancestor(node_t *node, node_t *cand_ancestor);

/**
 * @brief Print tree in DOT code on stdout and attach the leaves to the nodes
 *
 * @param root Root of the tree
 * @param tree Array of *node_t
 * @param leaves Leaves assignment (SIGMA)
 * @param MAX Number of cells
 */
void print_tree_leaves(node_t *root, node_t *tree[], int leaves[], int MAX,
                       double score);

/**
 * @brief Print tree in DOT code on outpath
 *
 * @param root Root of the (sub)tree
 * @param outpath Path to output file
 */
void fprint_tree(node_t *root, char *outpath, double score);

/**
 * @brief Print tree in DOT code on stdout and attach the leaves to the nodes
 *
 * @param root Root of the tree
 * @param tree Array of *node_t
 * @param leaves Leaves assignment (SIGMA)
 * @param MAX Number of cells
 * @param outpath Path to output file
 */
void fprint_tree_leaves(node_t *root, vector *tree_vec, int leaves[], int MAX,
                        char *outpath, double score, char cell_names[][255]);

/**
 * @brief Get the genotype profile of the selected node
 *
 * @param node Node to get the genotype profile
 * @param genotype Array int[MUTATIONS] initialized to all 0s
 */
void get_genotype_profile(node_t *node, int genotype[]);

/**
 * @brief Return a copy of the tree rooted in root. It also modifies the values
 * of tree[] and max_id
 *
 * @param root Root of the tree
 * @param tree Array of *node_t -- this will be modified accordingly
 * @return node_t* Pointer to the root of the new tree
 */
node_t *treecpy(node_t *root, vector *tree_vec, vector *losses_vec,
                vector *recs_vec, int n, int *original_mut_idx);

void node_detach(node_t *node);

/**
* @brief Check if the backmutation can be added
*
* @param loss Red node (backmutation)
* @return bool Boolean that rapresent if loss can be added to the tree
*/
bool is_loss_valid(node_t *loss);

/**
* @brief Check if all the nodes in the subtree, that lead to node,
* which have mutation index mut_index are valid
*
* @param node Red node (backmutation) that is being checked
* @param mut_index Mutation index of the nodes that are being checked
* @return bool Boolean that rapresent if the nodes checked are valid
*/
bool is_already_lost(node_t *node, int mut_index);

/**
* @brief Check if all the nodes in the subtree, that lead to recurrence,
* which have mutation index recurrece->mut_index are valid
*
* @param recurrence Green node (recurrent node)
* @return bool Boolean that rapresent if the nodes checked are valid
*/
bool is_recurrence_valid(node_t *recurrence);

/**
* @brief Check if all the nodes in the subtree under node are valid and if all
* the limits have been respected
*
* @param node Root of the tree or subtree that is being checked
* @param tree_vec Vector of *node_t
* @param og_muts_idx Array of int that save the index of all the original mutations
* @param m Number of mutations
* @param k Max number of losses for each mutation
* @param r Max number of recerrences for each mutation
* @param MAX_LOSSES Max number of deletions allowed in the solution
* @param MAX_RECURRENCES Max number of recurrences allowed in the solution
* @param loss Register that contain the number of deletions counted so far
* @param rec Register that contain the number of recerrences counted so far
* @param k_loss Array of int that contains the numbers of deletions of each mutation
* @param r_recs Array of int that contains the numbers of recurrences of each muation
* @return bool Boolean that rapresent if the subtree under node is valid
*/
bool is_tree_valid(node_t *node, vector *tree_vec, int *og_muts_idx, int m, int k,
                  int r, int MAX_LOSSES, int MAX_RECURRENCES, int* loss, int* rec,
                  int *k_loss, int *r_recs);

/**
* @brief Delete a node from the tree, reattach the parent and the childrens and
* delete the node from the vector
*
* @param node Node that is going to be deleted
* @param tree_vec Vector of *node_t that rapresent the tree
* @param vec Vector of *node_t
* @param mut Array of int
* @param sigma Array of int
* @param n Number of cells
*/
void node_delete(node_t *node, vector *tree_vec, vector *node_vec, int *mut,
                 int *sigma, int n);

#endif
