#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include "pipeline.h"



typedef struct readpairstruct {
    char *name;
    char *seq;
    char *qual;
    char *seq2;
    char *qual2;
    char *umi;
    size_t len;
    size_t len2;
    int fragment_size;
    int family;
    int copy_number;
    } ReadPair;



static char cpipeline_docstring[] = "This module provides a number of optimised functions for NGS pipelines including the Needleman-Wunch alignment algorithm.";
static char dedup_docstring[] = "Deduplicate and error correct FASTQ reads.";



static PyObject *dedup(PyObject *self, PyObject *args, PyObject *kwargs);



static PyMethodDef cpipeline_methods[] = {
    {"dedup", (PyCFunction)dedup, METH_VARARGS | METH_KEYWORDS, dedup_docstring},
    {NULL, NULL, 0, NULL}
    };



static struct PyModuleDef cpipeline_Module = {
    PyModuleDef_HEAD_INIT,
    "cpipeline",
    cpipeline_docstring,
    -1,
    cpipeline_methods
};



PyMODINIT_FUNC PyInit_cpipeline(void) {
    return PyModule_Create(&cpipeline_Module);
}



void reversecomplement(char *start);
void reverse(char *start);
static size_t strip_newlines(char *line, size_t length);
static int do_they_overlap(ReadPair *read, int min_overlap, int allowed);
static int compare_by_fragment_size_descending(const void *one, const void *two);
static int compare_by_family(const void *one, const void *two);
static int compare_by_sequence(const void *one, const void *two);
static int are_they_duplicates(ReadPair *read_a, ReadPair *read_b, int allowed);
static void print_sequences(ReadPair *read);
static int read_fasqs(char *read1, char *read2, int r1_umi_len, int r2_umi_len, ReadPair **_readpairs, char **_contents);



static void print_sequences(ReadPair *read) {
    printf("%s    %s\n", read->seq, read->seq2);
    }
    
    

static int compare_by_sequence(const void *first, const void *second) {
    ReadPair *item1 = (ReadPair *)first;
    ReadPair *item2 = (ReadPair *)second;
    int comp;

    comp = strcmp(item1->seq, item2->seq);
    if (comp == 0) {
        comp = strcmp(item1->seq2, item2->seq2);
        }
    return comp;
    }

    

static int compare_by_fragment_size_descending(const void *first, const void *second) {
    ReadPair *item1 = (ReadPair*)first;
    ReadPair *item2 = (ReadPair*)second;
    
    if (item2->fragment_size > item1->fragment_size)
        return 1;
    else if (item2->fragment_size < item1->fragment_size)
        return -1;
    else
        return 0;
    }
    


static int compare_by_family(const void *first, const void *second) {
    ReadPair *item1 = (ReadPair*)first;
    ReadPair *item2 = (ReadPair*)second;
    
    if (item2->family < item1->family)
        return 1;
    else if (item2->family > item1->family)
        return -1;
    else
        return 0;
    }
    


static size_t strip_newlines(char *line, size_t length) {
    if (length && line[length - 1] == '\n') {
        length -= 1;
        line[length] = '\0';
        if (length && line[length - 1] == '\r') {
            length -= 1;
            line[length] = '\0';
            }
        }
    return length;
    }


    
static int are_they_duplicates(ReadPair *read_a, ReadPair *read_b, int allowed) {
    int mismatches = 0;
    char *seq1, *seq2;
    
    seq1 = read_a->seq;
    seq2 = read_b->seq;
    while (*seq1 != '\0' && *seq2 != '\0') {
        if (*seq1 != *seq2  && *seq2 != 'N' && *seq2 != 'N') {
            mismatches += 1;
            if (mismatches > allowed)
                return 0;
            }
        seq1++;
        seq2++;
        }
    
    seq1 = read_a->seq2 + read_a->len2 - 1;
    seq2 = read_b->seq2 + read_b->len2 - 1;
    while (seq1 >= read_a->seq2 && seq2 >= read_b->seq2) {
        if (*seq1 != *seq2  && *seq2 != 'N' && *seq2 != 'N') {
            mismatches += 1;
            if (mismatches > allowed)
                return 0;
            }
        seq1--;
        seq2--;
        }
    
    return 1;
    }
    
    
// Do we need to make this handle 'N's in the input data?
static int do_they_overlap(ReadPair *read, int min_overlap, int allowed) {
    int read1_start, mismatches = 0;
    char *seq1, *seq2, *qual1, *qual2;
    
    if (min_overlap > read->len || min_overlap > read->len2)
        return 0;
    
    if (read->len <= read->len2) 
        read1_start = 0;
    else 
        read1_start = read->len - read->len2;
    
    // read1_start = index position in read1 to start comparing with read2[0]
    for(; read1_start < read->len - min_overlap + 1; read1_start++) {
        seq1 = read->seq + read1_start;
        seq2 = read->seq2;
        mismatches = 0;
        while (*seq1 != '\0' && *seq2 != '\0') {
            if (*seq1 != *seq2) {
                mismatches += 1;
                if (mismatches > allowed)
                    break;
                }
            seq1++;
            seq2++;
            }
        
        // It is a overlapping pair so 'N' out the mismatches and reduce quality to minimum.
        if (mismatches <= allowed) {
            seq1 = read->seq + read1_start;
            seq2 = read->seq2;
            qual1 = read->qual + read1_start;
            qual2 = read->qual2;
            while (*seq1 != '\0' && *seq2 != '\0') {
                if (*seq1 != *seq2) {
                    *seq1 = 'N';
                    *seq2 = 'N';
                    *qual1 = '!';
                    *qual2 = '!';
                    }
                seq1++;
                seq2++;
                qual1++;
                qual2++;
                }
            return read->len2 + read1_start;
            }
        }
    return 0;
    }
    

    
static PyObject *dedup(PyObject *self, PyObject *args, PyObject *kwargs) {
    // Expect two tuples of (name, seq, qual) as positional arguments
    //
    //
    
    PyObject *ret = NULL;

    ReadPair *readpairs = NULL;
    
    char *read1 = NULL, *read2 = NULL, *contents = NULL;
    int total_reads = 0, i = 0, j = 0, k = 0, n = 0, size = 0, matches = 0, next_family = 0, joined_family = 0, start_of_bin = 0, fragment_size = 0, last_fragment_size = -1;
    int unsized_start = 0, r1_umi_len = 0, r2_umi_len = 0;
    
    int64_t v_umi = 0, last_v_umi = 0;
    
    // Parse the input tuple
    static char *kwlist[] = {"read1", "read2", "r1_umi_len", "r2_umi_len", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|ii", kwlist, &read1, &read2, &r1_umi_len, &r2_umi_len))
        return NULL;

    if ((total_reads = read_fasqs(read1, read2, r1_umi_len, r2_umi_len, &readpairs, &contents)) == -1)
        goto error;
     
    
    // If the read is empty or consists entirely of 'N's then remove it as it is junk and will mess up our comparisons
    j = 0;
    n = 0;
    for (i = 0; i <  total_reads; i++) {
        for (k = 0; k < readpairs[i].len; k++) {
            if (readpairs[i].seq[k] != 'N')
                break;
            }
        if (k == readpairs[i].len) {
            n++;
            }
        else {
            j++;
            if (i > j)
                readpairs[j] = readpairs[i];
            }
        }
    printf("Removed %i N only reads.\n", n);
    total_reads -= n;
    

    // Remove exact duplicates.
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_sequence);
    j = 0;
    n = 0;
    for (i = 1; i < total_reads; i++) {
        if (are_they_duplicates(&readpairs[i], &readpairs[j], 0)) {            
            readpairs[j].copy_number += 1;
            n++;
            // Merge quality scores, keeping the highest.
            for (k =0; readpairs[i].qual[k] != '\0'; k++) {
                if (readpairs[j].qual[k] < readpairs[i].qual[k])
                    readpairs[j].qual[k] = readpairs[i].qual[k];
                }
            for (k =0; readpairs[i].qual2[k] != '\0'; k++) {
                if (readpairs[j].qual2[k] < readpairs[i].qual2[k])
                    readpairs[j].qual2[k] = readpairs[i].qual2[k];
                }
            }
        else {
            j++;
            if (i > j)
                readpairs[j] = readpairs[i];
            }
        }
    printf("Exact duplicates = %i of %i (%i%%)\n", n, total_reads, n * 100 / total_reads);
    total_reads -= n;
    
    
    
    // Fragment size all paired reads, 'N' out mismatches and sort largest first and unpaired (size=0) last.
    for (i = 0; i < total_reads; i++) {
        size = do_they_overlap(&readpairs[i], 100, 3);
        readpairs[i].fragment_size = size;
        if (size) {
            matches += 1;
            //printf("%i\n", size);
            }
        }
    printf("Sized = %i of %i (%i%\%)\n", matches, total_reads, matches * 100 / total_reads);
    goto finished;
    
    
    // Group sized fragments into families
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_fragment_size_descending);
    for (i = 0; i < total_reads && readpairs[i].fragment_size != 0; i++) {
        if (readpairs[i].family == 0)
            readpairs[i].family = ++next_family;
        fragment_size = readpairs[i].fragment_size;
        
        if (fragment_size != last_fragment_size) {
            start_of_bin = i;
            last_fragment_size = fragment_size;
            }
        
        for(j = i + 1; j < total_reads && readpairs[j].fragment_size == fragment_size; j++) {
            if (are_they_duplicates(&readpairs[i], &readpairs[j], 3)) {
                if (readpairs[j].family == 0)
                    readpairs[j].family = readpairs[i].family;
                else if (readpairs[j].family != readpairs[i].family) {
                    joined_family = readpairs[j].family;
                    for(n = start_of_bin; n < total_reads && readpairs[n].fragment_size == fragment_size; n++) { // Probably an inefficient loop
                        if (readpairs[n].family == joined_family)
                            readpairs[n].family = readpairs[i].family;
                        }
                    }
                }
            }
        }
    
    printf("i=%i, total=%i\n", i, total_reads);
    
    
    
    // Now work on unsized fragments
    v_umi = !((int64_t *)(readpairs[i].seq))[0];
    unsized_start = i;
    for (k = 0; k < 2; k++) {
        
        //qsort(readpairs + unsized_start, total_reads - unsized_start, sizeof(ReadPair), compare_umi1);
        for (i = unsized_start; i < total_reads; i++) {
            if (readpairs[i].family == 0)
                readpairs[i].family = ++next_family;
            v_umi = ((int64_t *)(readpairs[i].seq))[k];
            
            if (v_umi != last_v_umi) {
                start_of_bin = i;
                last_v_umi = v_umi;
                printf("%i\n", last_fragment_size);
                }
            
            for(j = i + 1; j < total_reads && ((int64_t *)(readpairs[j].seq))[k] == v_umi; j++) {
                if (are_they_duplicates(&readpairs[i], &readpairs[j], 3)) {
                    if (readpairs[j].family == 0)
                        readpairs[j].family = readpairs[i].family;
                    else if (readpairs[j].family != readpairs[i].family) {
                        joined_family = readpairs[j].family;
                        for(n = start_of_bin; n < total_reads && ((int64_t *)(readpairs[n].seq))[k] == v_umi; n++) { // Probably an inefficient loop
                            if (readpairs[n].family == joined_family)
                                readpairs[n].family = readpairs[i].family;
                            }
                        }
                    }
                }
            }
        }
    
    
    
    
    
    
    
    
    
    
    
    
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_family);
    
 /*   
    n = 0;
    family = 1;
    for (i = 0; i < total_reads; i++) {
        if (readpairs[i].family == 0) // only needed if above function not allowed to complete
            continue;

        if (readpairs[i].family == family)
            n += 1;
        else {
            printf(".\n.\n");//family = %i, members = %i\n", family, n);
            family = readpairs[i].family;
            n = 1;
            }
        print_sequences(&readpairs[i]);
        }*/
        
    finished:
    ret = Py_BuildValue("");
    error:
    free(readpairs);
    free(contents);
    return ret;
    }
    
    

static int read_fasqs(char *read1, char *read2, int r1_umi_len, int r2_umi_len, ReadPair **_readpairs, char **_contents) { 

    FILE *fp1 = NULL, *fp2 = NULL;
    
    ReadPair *readpairs;
    char *line = NULL, *contents = NULL;
    int total_reads = 0, i = 0, n = 0, retval = -1;
    ssize_t bytes_read = 0;
    size_t len = 0;
    long filesize = 0;

    
    // Open both fastqs, calculate file size and number of lines for memory allocation then reset file pointers to the start.
    fp1 = fopen(read1, "r");
    if (fp1 == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unable to open first file.");
        goto cleanup;
        }
    
    while ((bytes_read = getline(&line, &len, fp1)) != -1) {
        total_reads += 1;
        }
    filesize = ftell(fp1);
    total_reads = total_reads / 4;
    fseek(fp1, 0L, SEEK_SET);

    fp2 = fopen(read2, "r");
    if (fp2 == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unable to open second file.");
        goto cleanup;
        }
    fseek(fp2, 0L, SEEK_END);
    filesize += ftell(fp2);
    fseek(fp2, 0L, SEEK_SET);
    
    //printf("MB = %i, reads = %i\n", (int)filesize / 1024 / 1024, total_reads);

    
    // Allocate memory
    readpairs = (ReadPair *) malloc(sizeof(ReadPair) * total_reads);
    if (readpairs == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    *_readpairs = readpairs;
    contents = (char *) malloc(filesize + (total_reads * 8)); // only works with ascii files but hey fastqs are always ascii.
    if (contents == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    *_contents = contents;
    
    // Read paired fasqs into memory and reverse complement read2
    for (i = 0; i < total_reads; i++) {
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        strcpy(contents, line);
        readpairs[i].name = contents;
        contents += bytes_read + 1;
    
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        
        strcpy(contents, line);
        readpairs[i].seq = contents;
        contents += bytes_read + 1;
        readpairs[i].len = bytes_read;
    
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        strcpy(contents, line);
        readpairs[i].qual = contents;
        contents += bytes_read + 1;
        if (bytes_read != readpairs[i].len) {
            printf("r1, s=%i, q=%i\n",(int)bytes_read, (int)readpairs[i].len);
            PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
            goto cleanup;
            }
 
        
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        for (n = 0; n < bytes_read; n++) {
            if (line[n] != readpairs[i].name[n] && (line[n] != '2' || readpairs[i].name[n] != '1')) {
                PyErr_SetString(PyExc_TypeError, "Reads 1 and 2 names dont match.");
                goto cleanup;
                }
            }
    
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        strcpy(contents, line);
        readpairs[i].seq2 = contents;
        reversecomplement(contents);
        contents += bytes_read + 1;
        readpairs[i].len2 = bytes_read;
    
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        strcpy(contents, line);
        readpairs[i].qual2 = contents;
        reverse(contents);
        contents += bytes_read + 1;
        if (bytes_read != readpairs[i].len2) {
            PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
            goto cleanup;
            }
            
        readpairs[i].family = 0;
        readpairs[i].copy_number = 1;
        }

    retval = total_reads;
    cleanup:
    fclose(fp1);
    fclose(fp2);
    free(line);
    return retval;
    }



void reversecomplement(char *start) {
    if (start == NULL || *start == '\0')
        return;

    char *end = start;
    while (*end != '\0' && *end != '\n')
        end++;
    end -= 1;
    //size_t length = strlen(start);
    //char *end = start + length - 1;
    char temp ='\0';

    while (end > start) {
        switch (*start) {
            case 'A':
                temp = 'T'; break;
            case 'T':
                temp = 'A'; break;
            case 'C':
                temp = 'G'; break;
            case 'G':
                temp = 'C'; break;
            default:
                temp = *start; break;
            }
        switch (*end) {
            case 'A':
                *start = 'T'; break;
            case 'T':
                *start = 'A'; break;
            case 'C':
                *start = 'G'; break;
            case 'G':
                *start = 'C'; break;
            default:    
                *start = *end; break;
            }
        *end = temp;
        start++;
        end--;
        }
    if (start == end) {
        temp = *end;
        switch (temp) {
            case 'A':
                *end = 'T'; break;
            case 'T':
                *end = 'A'; break;
            case 'C':
                *end = 'G'; break;
            case 'G':
                *end = 'C'; break;
            default:
                *end = temp; break;
            }
        }
    }

    
    
void reverse(char *start) {
    if (start == NULL || *start == '\0')
        return;

    char *end = start;
    while (*end != '\0' && *end != '\n')
        end++;
    end -= 1;
    
    char temp ='\0';

    while (end > start) {
        temp = *start;
        *start = *end;
        *end = temp;
        start++;
        end--;
        }
    }

    
    
    
    
    
    