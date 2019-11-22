#include <Python.h>


typedef struct weightstruct {
    int gap;
    int gapstart;
    int match;
    int mismatch;
    } Weight;

typedef struct pathstruct {
    int score;
    int matrix;
    int direction;
    } Path;

typedef struct familystruct {
    char *seq;
    int allowed;
    char *seq2;
    int allowed22;
    PyObject *pylist;
    bool merged;
    } Family;



Path *initialise_matrix(char *testseq, size_t len_seq, char *refseq, size_t len_ref, int reverse, int ignore_gaps_before_seq, int ignore_gaps_after_seq,
                        int ignore_gaps_before_ref, int ignore_gaps_after_ref, Weight *weight, int *initial_i, int *initial_j, int *initial_m);

void reversecomplement(char *start);

