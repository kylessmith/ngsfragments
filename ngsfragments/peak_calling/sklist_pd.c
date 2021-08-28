
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sklist_pd.h"

float __skiplist_nanf(void) {
    const union {
        int __i;
        float __f;
    } __bint = {0x7fc00000UL};
    return __bint.__f;
}
#define PANDAS_NAN ((double)__skiplist_nanf())

double Log2(double val) { return log(val) / log(2.); }

double urand(void) {
    return ((double)rand() + 1) / ((double)RAND_MAX + 2);
}

int int_min(int a, int b) { return a < b ? a : b; }

node_t *node_init(double value, int levels) {
    node_t *result;
    result = (node_t *)malloc(sizeof(node_t));
    if (result) {
        result->value = value;
        result->levels = levels;
        result->is_nil = 0;
        result->ref_count = 0;
        result->next = (node_t **)malloc(levels * sizeof(node_t *));
        result->width = (int *)malloc(levels * sizeof(int));
        if (!(result->next && result->width) && (levels != 0)) {
            free(result->next);
            free(result->width);
            free(result);
            return NULL;
        }
    }
    return result;
}

// do this ourselves
void node_incref(node_t *node) { ++(node->ref_count); }

void node_decref(node_t *node) { --(node->ref_count); }

static void node_destroy(node_t *node) {
    int i;
    if (node) {
        if (node->ref_count <= 1) {
            for (i = 0; i < node->levels; ++i) {
                node_destroy(node->next[i]);
            }
            free(node->next);
            free(node->width);
            // printf("Reference count was 1, freeing\n");
            free(node);
        } else {
            node_decref(node);
        }
        // pretty sure that freeing the struct above will be enough
    }
}

void skiplist_destroy(skiplist_t *skp) {
    if (skp) {
        node_destroy(skp->head);
        free(skp->tmp_steps);
        free(skp->tmp_chain);
        free(skp);
    }
}

skiplist_t *skiplist_init(int expected_size) {
    skiplist_t *result;
    node_t *NIL, *head;
    int maxlevels, i;

    maxlevels = 1 + Log2((double)expected_size);
    result = (skiplist_t *)malloc(sizeof(skiplist_t));
    if (!result) {
        return NULL;
    }
    result->tmp_chain = (node_t **)malloc(maxlevels * sizeof(node_t *));
    result->tmp_steps = (int *)malloc(maxlevels * sizeof(int));
    result->maxlevels = maxlevels;
    result->size = 0;

    head = result->head = node_init(PANDAS_NAN, maxlevels);
    NIL = node_init(0.0, 0);

    if (!(result->tmp_chain && result->tmp_steps && result->head && NIL)) {
        skiplist_destroy(result);
        node_destroy(NIL);
        return NULL;
    }

    node_incref(head);

    NIL->is_nil = 1;

    for (i = 0; i < maxlevels; ++i) {
        head->next[i] = NIL;
        head->width[i] = 1;
        node_incref(NIL);
    }

    return result;
}

// 1 if left < right, 0 if left == right, -1 if left > right
int _node_cmp(node_t *node, double value) {
    if (node->is_nil || node->value > value) {
        return -1;
    } else if (node->value < value) {
        return 1;
    } else {
        return 0;
    }
}

double skiplist_get(skiplist_t *skp, int i) {
    node_t *node;
    int level;

    if (i < 0 || i >= skp->size) {
        return 0;
    }

    node = skp->head;
    ++i;
    for (level = skp->maxlevels - 1; level >= 0; --level) {
        while (node->width[level] <= i) {
            i -= node->width[level];
            node = node->next[level];
        }
    }

    return node->value;
}

int skiplist_insert(skiplist_t *skp, double value) {
    node_t *node, *prevnode, *newnode, *next_at_level;
    int *steps_at_level;
    int size, steps, level;
    node_t **chain;

    chain = skp->tmp_chain;

    steps_at_level = skp->tmp_steps;
    memset(steps_at_level, 0, skp->maxlevels * sizeof(int));

    node = skp->head;

    for (level = skp->maxlevels - 1; level >= 0; --level) {
        next_at_level = node->next[level];
        while (_node_cmp(next_at_level, value) >= 0) {
            steps_at_level[level] += node->width[level];
            node = next_at_level;
            next_at_level = node->next[level];
        }
        chain[level] = node;
    }

    size = int_min(skp->maxlevels, 1 - ((int)Log2(urand())));

    newnode = node_init(value, size);
    if (!newnode) {
        return -1;
    }
    steps = 0;

    for (level = 0; level < size; ++level) {
        prevnode = chain[level];
        newnode->next[level] = prevnode->next[level];

        prevnode->next[level] = newnode;
        node_incref(newnode);  // increment the reference count

        newnode->width[level] = prevnode->width[level] - steps;
        prevnode->width[level] = steps + 1;

        steps += steps_at_level[level];
    }

    for (level = size; level < skp->maxlevels; ++level) {
        chain[level]->width[level] += 1;
    }

    ++(skp->size);

    return 1;
}

int skiplist_remove(skiplist_t *skp, double value) {
    int level, size;
    node_t *node, *prevnode, *tmpnode, *next_at_level;
    node_t **chain;

    chain = skp->tmp_chain;
    node = skp->head;

    for (level = skp->maxlevels - 1; level >= 0; --level) {
        next_at_level = node->next[level];
        while (_node_cmp(next_at_level, value) > 0) {
            node = next_at_level;
            next_at_level = node->next[level];
        }
        chain[level] = node;
    }

    if (value != chain[0]->next[0]->value) {
        return 0;
    }

    size = chain[0]->next[0]->levels;

    for (level = 0; level < size; ++level) {
        prevnode = chain[level];

        tmpnode = prevnode->next[level];

        prevnode->width[level] += tmpnode->width[level] - 1;
        prevnode->next[level] = tmpnode->next[level];

        tmpnode->next[level] = NULL;
        node_destroy(tmpnode);  // decrement refcount or free
    }

    for (level = size; level < skp->maxlevels; ++level) {
        --(chain[level]->width[level]);
    }

    --(skp->size);
    return 1;
}

// Calculate median
double find_median(skiplist_t *sklist)
{
    int median_index = sklist->size / 2;
	if(sklist->size % 2 == 0)
	{
        double median1 = skiplist_get(sklist, median_index);
        double median2 = skiplist_get(sklist, median_index-1);
		return (median1 + median2) / 2.0;
	}
	else
    {
        double median = skiplist_get(sklist, median_index);
		return median/1.0;
    }
};

// Calculate rolling median for entire array
void rolling_median(const double values[], double medians[], int n, int window)
{
	skiplist_t *sklist = skiplist_init(window);
	
	int half_window = window / 2;
	int i;
	// Initial skip list insertions
	for(i = 0; i < window; i++) 
	{
		skiplist_insert(sklist, values[i]);
	}
	
	int k;
	// Iterate over values for rolling median
	for(k = half_window; k < n - half_window; k++)
	{
		medians[k] = find_median(sklist);
		skiplist_remove(sklist, values[k - half_window]);
		skiplist_insert(sklist, values[k + half_window]);
	}

	skiplist_destroy(sklist);
};

// Calculate rolling median from a non-zero start
void rolling_median_shifted(const double values[], double medians[], int n, int window, int start)
{
	skiplist_t *sklist = skiplist_init(window);
	
	int half_window = window / 2;
	int i;
	// Initial skip list insertions
	for(i = 0+start; i < start+window; i++) 
	{
		skiplist_insert(sklist, values[i]);
	}
	
	int k;
	// Iterate over values for rolling median
	for(k = half_window+start; k < n - half_window; k++)
	{
		medians[k] = find_median(sklist);
		skiplist_remove(sklist, values[k - half_window]);
		skiplist_insert(sklist, values[k + half_window]);
	}
	
	skiplist_destroy(sklist);
};

// Display skip list level wise 
void displayList(skiplist_t *sklist) 
{ 
    printf("\n*****Skip List values*****\n");
	int i;
	node_t *node;
    for(i=0; i<sklist->maxlevels; i++) 
    { 
        node = sklist->head->next[i]; 
        printf("Level %d: ", i); 
        while(node != NULL) 
        { 
            printf("%f ", node->value); 
            node = node->next[i]; 
        } 
        printf("\n"); 
    } 
}; 

// Display skip list widths level wise 
void displayWidths(skiplist_t *sklist) 
{ 
    printf("\n*****Skip List widths*****\n");
	int i;
	node_t *node;
    for(i=0; i<sklist->maxlevels; i++) 
    { 
        node = sklist->head->next[i]; 
        printf("Level %d: ", i); 
        while(node != NULL) 
        { 
            printf("%d ", node->width[i]); 
            node = node->next[i]; 
        } 
        printf("\n"); 
    } 
}; 
