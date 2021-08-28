#ifndef SKLIST_PD_H
#define SKLIST_PD_H

typedef struct node_t node_t;

struct node_t {
    node_t **next;
    int *width;
    double value;
    int is_nil;
    int levels;
    int ref_count;
};

typedef struct {
    node_t *head;
    node_t **tmp_chain;
    int *tmp_steps;
    int size;
    int maxlevels;
} skiplist_t;

float __skiplist_nanf(void);
double Log2(double val);
double urand(void);
int int_min(int a, int b);
node_t *node_init(double value, int levels);
void node_incref(node_t *node);
void node_decref(node_t *node);
static void node_destroy(node_t *node);
void skiplist_destroy(skiplist_t *skp);
skiplist_t *skiplist_init(int expected_size);
int _node_cmp(node_t *node, double value);
double skiplist_get(skiplist_t *skp, int i);
int skiplist_insert(skiplist_t *skp, double value);
int skiplist_remove(skiplist_t *skp, double value);
double find_median(skiplist_t *sklist);
void rolling_median(const double values[], double medians[], int n, int window);
void rolling_median_shifted(const double values[], double medians[], int n, int window, int start);
void displayList(skiplist_t *sklist);
void displayWidths(skiplist_t *sklist);


#endif