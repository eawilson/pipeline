#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>

#define THRUPLEX_UMT 6
#define THRUPLEX_STEM 11
#define SIGNIFICANT_PHRED_DIFFERENCE 10
#define MINIMUM_NON_N_BASES 50
#define R1 0
#define R2 1



typedef struct readstruct {
    char *name;
    char *seq;
    char *qual;
    char *umi;
    unsigned short len;
    unsigned short nonoverlapping_len; // Only meaningful for R2
    } Read;

    
    
typedef struct readpairstruct {
    Read read[2];
    unsigned short fragment_size;
    unsigned short copy_number;
    int family;
    int prevfamily;
    } ReadPair;

    
    
typedef struct pairedfastqstruct {
    ReadPair *readpairs;
    int total_reads;
    } PairedFastq;

    
    
typedef struct mergematrixstruct {
    int first_match;
    int second_match;
    int swap;
    } MergeMatrix;

    
    
static char cpipeline_docstring[] = "This module provides a number of optimised functions for NGS pipelines.";
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



static void reversecomplement(char *start, int len);
static void reverse(char *start, int len);
static size_t strip_newlines(char *line, size_t length);
static int do_they_overlap(ReadPair *read, int min_overlap, int allowed);
static int compare_by_fragment_size_descending(const void *one, const void *two);
static int compare_by_family(const void *one, const void *two);
static int compare_by_sequence(const void *one, const void *two);
static int compare_by_short_sequence(const void *first, const void *second, const void *offset);
static int are_they_duplicates(ReadPair *readpair1, ReadPair *readpair2, int allowed);
static PairedFastq read_fastqs(char *read1, char *read2);
static void debug_print_read(ReadPair *readpair);
static void debug_print_fastq(PairedFastq fq);
static void brute_assign_families(ReadPair *bin_start, ReadPair *bin_end, int *current_family, int allowed);
static int remove_exact_duplicates(PairedFastq *fq);
static int remove_n_only_reads(PairedFastq *fq);
static int size_and_remove_umis(PairedFastq *fq, int allowed, int thruplex);
static int assign_families(PairedFastq *fq, int allowed);
static void print_family_sizes(PairedFastq *fq);
static int collapse_families(PairedFastq *fq);
static void merge_reads(ReadPair *bin_start, ReadPair * bin_end);
static int remove_unconfirmed_reads(PairedFastq *fq);
static int write_fastqs(PairedFastq fq, char *read1, char *read2);
static void print_family_sequences(PairedFastq *fq);



static void brute_assign_families(ReadPair *bin_start, ReadPair *bin_end, int *current_family, int allowed) {
    ReadPair *readpair, *readpair2, *readpair3;
    int joined_family;
    //printf("bin size = %i, allowed = %i\n", bin_end - bin_start, allowed);
    for (readpair = bin_start; readpair < bin_end; readpair++) {
        if (readpair->family == 0)
            readpair->family = ++(*current_family);
        
        for(readpair2 = readpair + 1; readpair2 < bin_end; readpair2++) {
            //printf("comparing\n%s\n%s\n", readpair->read[R1].seq, readpair2->read[R1].seq);
            if (are_they_duplicates(readpair, readpair2, allowed)) {
                //printf("match\n");
                if (readpair2->family == 0)
                    readpair2->family = readpair->family;
                else if (readpair2->family != readpair->family) {
                    joined_family = readpair2->family;
                    for(readpair3 = bin_start; readpair3 < bin_end; readpair3++) {
                        if (readpair3->family == joined_family)
                            readpair3->family = readpair->family;
                        }
                    }
                }
            }
        }
    }
    
    
    
// DO WE NEED TO WORRY ABOUT OVERSHOOTING ON SHORT READS?
static int compare_by_short_sequence(const void *first, const void *second, const void *offset) {
    ReadPair *item1 = (ReadPair *)first;
    ReadPair *item2 = (ReadPair *)second;
    int index = 10 + (*((int *)offset) * 6);
    int comp;

    comp = memcmp(item1->read[R1].seq + index, item2->read[R1].seq + index, 6);
    if (comp == 0) {
        comp = memcmp(item1->read[R2].seq + index, item2->read[R2].seq + index, 6);
        }
    return comp;
    }

    

static int compare_by_sequence(const void *first, const void *second) {
    ReadPair *item1 = (ReadPair *)first;
    ReadPair *item2 = (ReadPair *)second;
    int comp;

    comp = strcmp(item1->read[R1].seq, item2->read[R1].seq);
    if (comp == 0) {
        comp = strcmp(item1->read[R2].seq, item2->read[R2].seq);
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


    
static int are_they_duplicates(ReadPair *readpair1, ReadPair *readpair2, int allowed) {
    int mismatches = 0, i = 0, len = 0;
    
    if (readpair1->read[R1].nonoverlapping_len > readpair2->read[R1].nonoverlapping_len)
        len = readpair1->read[R1].nonoverlapping_len;
    else
        len = readpair2->read[R1].nonoverlapping_len;
    for (i = 0; i < len; i++) {
        if (readpair1->read[R1].seq[i] != readpair2->read[R1].seq[i] && readpair1->read[R1].seq[i] != 'N' && readpair2->read[R1].seq[i] != 'N') {
            mismatches++;
            if (mismatches > allowed)
                return 0;
            }
        }
        
    if (readpair1->read[R2].len < readpair2->read[R2].len)
        len = readpair1->read[R2].len;
    else
        len = readpair2->read[R2].len;
    for (i = 0; i < len; i++) {
        if (readpair1->read[R2].seq[i] != readpair2->read[R2].seq[i] && readpair1->read[R2].seq[i] != 'N' && readpair2->read[R2].seq[i] != 'N') {
            mismatches++;
            if (mismatches > allowed)
                return 0;
            }
        }
    return 1;
    }

    
    
static int do_they_overlap(ReadPair *readpair, int min_overlap, int allowed) {
    int read1_start = 0, best_read1_start = 0, mismatches = 0, best = allowed + 1;
    char *seq1 = NULL, *seq2 = NULL, *qual1 = NULL, *qual2 = NULL;
    if (min_overlap > readpair->read[R1].len || min_overlap > readpair->read[R2].len)
        return 0;
    
    if (readpair->read[R1].len <= readpair->read[R2].len) 
        read1_start = 0;
    else 
        read1_start = readpair->read[R1].len - readpair->read[R2].len;
    
    // read1_start = index position in read1 to start comparing with read2[0]
    for(; read1_start < readpair->read[R1].len - min_overlap + 1; read1_start++) {
        seq1 = readpair->read[R1].seq + read1_start;
        seq2 = readpair->read[R2].seq;
        mismatches = 0;
        while (*seq1 != '\0' && *seq2 != '\0') {
            if (*seq1 != *seq2 && *seq1 != 'N' && *seq2 != 'N') {
                mismatches += 1;
                if (mismatches > allowed)
                    break;
                }
            seq1++;
            seq2++;
            }
        // Don't need to store best as we quit after first match, but just in case we ever want to search to the end...
        if (mismatches < best) {
            best = mismatches;
            best_read1_start = read1_start;
            break;
            }
        }
        
    // It is a overlapping pair so 'N' out the mismatches and reduce quality to minimum.
    if (best <= allowed) {
        seq1 = readpair->read[R1].seq + best_read1_start;
        seq2 = readpair->read[R2].seq;
        qual1 = readpair->read[R1].qual + best_read1_start;
        qual2 = readpair->read[R2].qual;
        while (*seq1 != '\0' && *seq2 != '\0') {
            if (*seq1 != *seq2) {
                if (*qual1 > *qual2 + SIGNIFICANT_PHRED_DIFFERENCE) {
                    *seq2 = *seq1;
                    *qual2 = *qual1;
                    }
                else if (*qual2 > *qual1 + SIGNIFICANT_PHRED_DIFFERENCE) {
                    *seq1 = *seq2;
                    *qual1 = *qual2;
                    }
                else {
                    *seq1 = 'N';
                    *seq2 = 'N';
                    *qual1 = '!';
                    *qual2 = '!';
                    }
                }
            seq1++;
            seq2++;
            qual1++;
            qual2++;
            }
            
        readpair->read[R1].nonoverlapping_len = best_read1_start;
        readpair->fragment_size = readpair->read[R2].len + best_read1_start;
        return 1;
        }
    return 0;
    }


    
static PyObject *dedup(PyObject *self, PyObject *args, PyObject *kwargs) {
    /*
     * 
     */
    
    PyObject *ret = NULL;
    PairedFastq fq;
    
    clock_t time_start = clock();    
    char *read1 = NULL, *read2 = NULL;
    int thruplex = 0, allowed = 3,  n = 0;

    // Parse the input tuple
    static char *kwlist[] = {"read1", "read2", "allowed", "thruplex", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|ip", kwlist, &read1, &read2, &allowed, &thruplex))
        return NULL;

    // Read the paired fastqs
    fq = read_fastqs(read1, read2);
    if (fq.total_reads == 0)
        goto error;

    // Remove N only and too short reads
    n = remove_n_only_reads(&fq);
    fprintf(stderr, "Removed %i N only reads (%i%% of %i)\n", n, n * 100 / (fq.total_reads + n), fq.total_reads + n);

    // Remove exact duplicates
    n = remove_exact_duplicates(&fq);
    fprintf(stderr, "Removed %i exact duplicates (%i%% of %i)\n", n, n * 100 / (fq.total_reads + n), fq.total_reads + n);
    
    // Size fragments and remove thruplex umis
    n = size_and_remove_umis(&fq, allowed, thruplex);
    fprintf(stderr, "Sized %i reads (%i%% of %i)\n", n, n * 100 / fq.total_reads, fq.total_reads);
    
    // Group reads into families
    if (assign_families(&fq, allowed) == -1)
        goto error;

    // Print family sizes
    print_family_sizes(&fq);
    //print_family_sequences(&fq);
    
    // Collapse families into single consensus reads
    n = collapse_families(&fq);
    fprintf(stderr, "Removed %i inexact duplicates (%i%% of %i)\n", n, n * 100 / (fq.total_reads + n), fq.total_reads + n);
    
    // Remove reads with no size or family
    n = remove_unconfirmed_reads(&fq);
    fprintf(stderr, "Removed %i reads without size or sibligs (%i%% of %i)\n", n, n * 100 / (fq.total_reads + n), fq.total_reads + n);

    // Write remaining reads to file
    fprintf(stderr, "Writing %i reads to file\n", fq.total_reads);
    if (write_fastqs(fq, read1, read2) == -1)
        goto error;
    
    finished:
    fprintf(stderr, "Completed in %f seconds\n", (double)(clock() - time_start) / CLOCKS_PER_SEC);
    ret = Py_BuildValue("");
    
    error:
    free(fq.readpairs);
    return ret;
    }


    
static int write_fastqs(PairedFastq fq, char *read1, char *read2) {
    FILE *fp = NULL;
    int read = R1;
    ReadPair *readpair = NULL;
    char *filenames[2] = {read1, read2}, *buffer = NULL;
    
    buffer = malloc(sizeof(char) * (strlen(filenames[R1]) + strlen(filenames[R2]) + 9));
    if (buffer == NULL) {
        PyErr_NoMemory();
        return -1;
        }

    for (read = R1; read <= R2; read++) {
        strcpy(buffer, filenames[read]);
        strcpy(buffer + strlen(filenames[read]) - 6, ".deduped.fastq");
        
        fp = fopen(buffer, "w");
        if (fp == NULL) {
            PyErr_SetString(PyExc_TypeError, "Unable to open output file.");
            free(buffer);
            return -1;
            }
        for (readpair = fq.readpairs; readpair < fq.readpairs + fq.total_reads; readpair++)
            fprintf(fp, "%s\n%s\n+\n%s\n", readpair->read[read].name, readpair->read[read].seq, readpair->read[read].qual);
        fclose(fp);
        }
        
    free(buffer);
    return 0;
    }

    
    
static PairedFastq read_fastqs(char *read1, char *read2) { 
    /* Read paired fastq files into an array of ReadPair structs. 
     * Strip newlines (\n or \r\n) from all strings.
     * Ensure R1 sequence and quality are the same length and that R2 sequence and quality are the same length.
     * Ensure that R1 and R2 names only differ where there is a '1' in the R1 name and a '2' in the R2 name.
     */
    PairedFastq fq = {.readpairs = NULL, .total_reads = 0};
    ReadPair *readpair;
    FILE *fp[2] = {NULL, NULL};
    char *contents = NULL, *line = NULL, *filenames[2] = {read1, read2};
    int j = 0, read = R1;
    ssize_t bytes_read = 0;
    size_t line_len = 0;
    long filesize = 0;    
    
    // Open both fastqs, calculate file size and number of lines for memory allocation then reset file pointers to the start.
    for (read = R1; read <= R2; read++) {
        if (strcmp(filenames[read] + strlen(filenames[read]) - 6 , ".fastq") != 0) {
            PyErr_SetString(PyExc_TypeError, "Not a fastq.");
            goto cleanup;
            }        
        fp[read] = fopen(filenames[read], "r");
        if (fp[read] == NULL) {
            PyErr_SetString(PyExc_TypeError, "Unable to open fastq.");
            goto cleanup;
            }
        fseek(fp[read], 0L, SEEK_END);
        filesize += ftell(fp[read]);
        fseek(fp[read], 0L, SEEK_SET);
        }
           
    while ((bytes_read = getline(&line, &line_len, fp[R1])) != -1)
        fq.total_reads += 1;
    fq.total_reads = fq.total_reads / 4;
    fseek(fp[R1], 0L, SEEK_SET);
   
    fprintf(stderr, "Filesize = %ld MB, Reads = %i\n", filesize / 1024 / 1024, fq.total_reads);

    
    // Allocate memory
    fq.readpairs = malloc((fq.total_reads * sizeof(ReadPair)) + filesize);
    if (fq.readpairs == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    memset(fq.readpairs, 0, fq.total_reads * sizeof(ReadPair));
    contents = (char *)(fq.readpairs + fq.total_reads);
    
    // Read paired fasqs into memory
    for (readpair = fq.readpairs; readpair < fq.readpairs + fq.total_reads; readpair++) {
        for (read = R1; read <= R2; read++) {
            // Name
            if ((bytes_read = getline(&line, &line_len, fp[read])) == -1) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }
            strip_newlines(line, bytes_read);
            memcpy(contents, line, bytes_read + 1);
            readpair->read[read].name = contents;
            contents += bytes_read + 1;
            
            // Sequence
            if ((bytes_read = getline(&line, &line_len, fp[read])) == -1) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }
            bytes_read = strip_newlines(line, bytes_read);
            memcpy(contents, line, bytes_read + 1);
            readpair->read[read].seq = contents;
            readpair->read[read].len = bytes_read;
            contents += bytes_read + 1;
        
            // Plus
            if ((bytes_read = getline(&line, &line_len, fp[read])) == -1) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }
        
            // Quality
            if ((bytes_read = getline(&line, &line_len, fp[read])) == -1) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }
            bytes_read = strip_newlines(line, bytes_read);
            memcpy(contents, line, bytes_read + 1);
            readpair->read[read].qual = contents;
            contents += bytes_read + 1;
            if (bytes_read != readpair->read[read].len) {
                PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
                goto cleanup;
                }
            }
            
        for (j = 0; !(readpair->read[R1].name[j] == '\0' && readpair->read[R2].name[j] == '\0'); j++) {
            if (readpair->read[R1].name[j] != readpair->read[R2].name[j] && !(readpair->read[R1].name[j] == '1' && readpair->read[R2].name[j] == '2')) {
                PyErr_SetString(PyExc_TypeError, "Reads 1 and 2 names dont match.");
                goto cleanup;
                }
            }
        
        readpair->read[R1].nonoverlapping_len = readpair->read[R1].len;
        readpair->copy_number = 1;
        }
    
    goto finished;
    
    cleanup:
    fq.total_reads = 0;
    
    finished:
    if (fp[R1] != NULL)
        fclose(fp[R1]);
    if (fp[R2] != NULL)
        fclose(fp[R2]);
    free(line);
    return fq;
    }



static void reversecomplement(char *start, int len) {
    if (len == 0)
        return;

    char *end = start + len - 1;
    char temp = '\0';

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

    
    
static void reverse(char *start, int len) {
    if (len == 0)
        return;

    char *end = start + len - 1;
    char temp ='\0';

    while (end > start) {
        temp = *start;
        *start = *end;
        *end = temp;
        start++;
        end--;
        }
    }

    

static void debug_print_read(ReadPair *readpair) {
    /* Print out all the members of a readpair struct, for use in debugging only.
     */ 
    fprintf(stderr, "Name   = %s\n", readpair->read[R1].name);
    fprintf(stderr, "Name2  = %s\n", readpair->read[R2].name);
    fprintf(stderr, "Seq    = %s\n", readpair->read[R1].seq);
    fprintf(stderr, "Seq2   = %s\n", readpair->read[R2].seq);
    fprintf(stderr, "Qual   = %s\n", readpair->read[R1].qual);
    fprintf(stderr, "Qual2  = %s\n", readpair->read[R2].qual);
    fprintf(stderr, "Umi1   = %s\n", readpair->read[R1].umi);
    fprintf(stderr, "Umi2   = %s\n", readpair->read[R2].umi);
    fprintf(stderr, "Len    = %i\n", (int)(readpair->read[R1].len));
    fprintf(stderr, "Len2   = %i\n", (int)(readpair->read[R2].len));
    fprintf(stderr, "NO Len = %i\n", (int)(readpair->read[R1].nonoverlapping_len));
    fprintf(stderr, "Size   = %i\n", readpair->fragment_size);
    fprintf(stderr, "Family = %i\n", readpair->family);
    fprintf(stderr, "Copies = %i\n", readpair->copy_number);
    }


    
static void debug_print_fastq(PairedFastq fq) {
    /* Print out every readpair, for use in debugging only.
     * WARNING - This will print out millions of lines if used with normal input.
     */ 
    ReadPair *readpair = NULL;
    for (readpair = fq.readpairs; readpair < fq.readpairs + fq.total_reads; readpair++) {
        debug_print_read(readpair);
        fprintf(stderr, ".\n");
        }
    }
    
    

static int remove_exact_duplicates(PairedFastq *fq) {
    /* Remove exact duplicates. Ns don't count in the comparison but the reads do have to be the same length. 
     * Therefore they may be duplicates even if they are not identical. Merge duplicate reads, correcting Ns 
     * where possible and keeping the higher quality score.
     * WILL MISS SOME EXACT DUPLICATES IF AN N INTEREFERES WITH THE SORT however these should be caught later
     * and the purpose of this function is to reduce the burden on later functionns rather than catch everything.
     */
    int n = 0, k = 0, read = R1;
    ReadPair *readpair = NULL, *readpair2 = NULL;
    
    qsort(fq->readpairs, fq->total_reads, sizeof(ReadPair), (void *)compare_by_sequence);
    readpair2 = fq->readpairs;
    for (readpair = fq->readpairs + 1; readpair < fq->readpairs + fq->total_reads; readpair++) {
        if (readpair->read[R1].len == readpair2->read[R1].len && readpair->read[R2].len == readpair2->read[R2].len && are_they_duplicates(readpair, readpair2, 0)) {
            readpair2->copy_number += 1;
            n++;
            
            // Merge qualities, keeping the highest and try to correct 'N's.
            for (read = R1; read <= R2; read++) {
                for (k = 0; k < readpair2->read[read].len; k++) {
                    if (readpair2->read[read].seq[k] == 'N') {
                        readpair2->read[read].seq[k] = readpair->read[read].seq[k];
                        readpair2->read[read].qual[k] = readpair->read[read].qual[k];                
                        }
                    else if (readpair2->read[read].qual[k] < readpair->read[read].qual[k])
                        readpair2->read[read].qual[k] = readpair->read[read].qual[k];                
                    }
                }
            }
        else {
            readpair2++;
            if (readpair > readpair2)
                *readpair2 = *readpair;
            }
        }
    fq->total_reads -= n;
    return n;
    }

    

static int remove_n_only_reads(PairedFastq *fq) {
    /* Remove reads consisting of only Ns in either the forward or reverse read as these will match everything.
     * Also removes reads with less than MINIMUM_NON_N_BASES of non N bases.
     */ 
    int n = 0, i = 0, bases = 0, read = R1;
    ReadPair *readpair = NULL, *readpair2 = NULL;

    readpair2 = fq->readpairs;
    for (readpair = fq->readpairs; readpair <  fq->readpairs + fq->total_reads; readpair++) {
        for (read = R1; read <= R2; read++) {
            bases = 0;
            for (i = 0; i < readpair->read[read].len; i++) {
                if (readpair->read[read].seq[i] != 'N') {
                    bases++;
                    if (bases == MINIMUM_NON_N_BASES)
                        break;
                    }
                }
            if (bases < MINIMUM_NON_N_BASES)
                break;
            }
        
        if (bases < MINIMUM_NON_N_BASES) {
            n++;
            }
        else {
            if (readpair > readpair2)
                *readpair2 = *readpair;
            readpair2++;
            }
        }
    fq->total_reads -= n;
    return n;
    }
    
    
    
static int size_and_remove_umis(PairedFastq *fq, int allowed, int thruplex) {
    /* Fragment size paired reads, Correct sequencing errors if phred difference large enough or 'N' them out otherwise.
     * If thruplex then trim the umis including if there has been runthrough sequencing into the stem / umi at the other end.
     * Makes sense to deal with sizing and umis together as they both need a reverse complement read 2.
     */
    int n = 0, r1_readthrough = 0, r2_readthrough = 0;
    ReadPair *readpair = NULL;

    for (readpair = fq->readpairs; readpair < fq->readpairs + fq->total_reads; readpair++) {
        reversecomplement(readpair->read[R2].seq, readpair->read[R2].len);
        reverse(readpair->read[R2].qual, readpair->read[R2].len);
        
        n += do_they_overlap(readpair, 70, allowed);

        if (thruplex) {
            if (readpair->fragment_size > 0) {
                r2_readthrough = THRUPLEX_UMT + THRUPLEX_STEM - readpair->read[R1].nonoverlapping_len;
                r1_readthrough = THRUPLEX_UMT + THRUPLEX_STEM - (readpair->read[R2].len - (readpair->read[R1].len - readpair->read[R1].nonoverlapping_len));
                if (r2_readthrough > 0) {
                    readpair->read[R2].seq += r2_readthrough;
                    readpair->read[R2].qual += r2_readthrough;
                    readpair->read[R2].len -= r2_readthrough;
                    readpair->read[R1].nonoverlapping_len += r2_readthrough; // The overhang at the start of R1 has increased but will be zeroed when the umis are subtracted
                    }
                if (r1_readthrough > 0) {
                    readpair->read[R1].len -= r1_readthrough;
                    readpair->read[R1].seq[readpair->read[R1].len] = '\0';
                    readpair->read[R1].qual[readpair->read[R1].len] = '\0';
                    }
                }
            
            readpair->read[R1].umi = readpair->read[R1].seq;
            readpair->read[R1].umi[THRUPLEX_UMT] = '\0';
            readpair->read[R1].len -= THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->read[R1].nonoverlapping_len -= THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->read[R1].seq += THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->read[R1].qual += THRUPLEX_UMT + THRUPLEX_STEM;
            
            readpair->read[R2].umi = readpair->read[R2].seq + readpair->read[R2].len - THRUPLEX_UMT;
            readpair->read[R2].len -= THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->read[R2].seq[readpair->read[R2].len] = '\0';
            readpair->read[R2].qual[readpair->read[R2].len] = '\0'; // only needed if quality is not lazy loaded
            
            if (readpair->fragment_size > 0) // Can't subtract size from unsized pairs
                readpair->fragment_size -= 2 * (THRUPLEX_UMT + THRUPLEX_STEM);
            }
            
        reversecomplement(readpair->read[R2].seq, readpair->read[R2].len);
        reverse(readpair->read[R2].qual, readpair->read[R2].len);
        }
    return n;
    }

    

static int assign_families(PairedFastq *fq, int allowed) {
    int biggest_bin = 0, current_family = 0, temp_family = 0, family = 0, max_family = 0, n = 0, i = 0, incomplete = 0, merge_required = 0;

    ReadPair *bin_start = NULL, *bin_end = NULL, *subbin_start = NULL, *subbin_end = NULL;
    MergeMatrix *mergematrix = NULL;
    
    qsort(fq->readpairs, fq->total_reads, sizeof(ReadPair), (void *)compare_by_fragment_size_descending);

    // Calculate size of largest bin to allocate memory.
    bin_start = fq->readpairs;
    while (bin_start < fq->readpairs + fq->total_reads) {
        for (bin_end = bin_start + 1; bin_end < fq->readpairs + fq->total_reads; bin_end++) {
            if (bin_end->fragment_size != bin_start->fragment_size)
                break;
            }
        if (bin_end - bin_start > biggest_bin)
            biggest_bin = bin_end - bin_start;
        bin_start = bin_end;
        }
    mergematrix = malloc(biggest_bin * sizeof(MergeMatrix));
    if (mergematrix == NULL) {
        PyErr_NoMemory();
        return -1;
        }
        
    // Do the actual family assessment.
    bin_start = fq->readpairs;
    while (bin_start < fq->readpairs + fq->total_reads) {
        for (bin_end = bin_start + 1; bin_end < fq->readpairs + fq->total_reads; bin_end++) {
            if (bin_end->fragment_size != bin_start->fragment_size)
                break;
            }
            
        if (bin_end - bin_start < 2000) {
            brute_assign_families(bin_start, bin_end, &current_family, allowed);
            }
            
        else {
            for (n = 0; n < allowed + 4; n++) {
                temp_family = 0;
                qsort_r(bin_start, bin_end - bin_start, sizeof(ReadPair), (void *)compare_by_short_sequence, &n);
                subbin_start = bin_start;
                while (subbin_start < bin_end) {
                    for (subbin_end = subbin_start + 1; subbin_end < bin_end; subbin_end++) {
                        if (compare_by_short_sequence(subbin_start, subbin_end, &n) != 0)
                            break;
                        }
                    brute_assign_families(subbin_start, subbin_end, &temp_family, allowed);
                    subbin_start = subbin_end;
                    }
                
                if (n == 0) { // Don't need to merge families on first iteration.
                    for (i = 0; bin_start + i < bin_end; i++) {
                        bin_start[i].prevfamily = bin_start[i].family;
                        bin_start[i].family = 0;
                        }
                    }
                    
                else {
                    do {
                        memset(mergematrix, 0, biggest_bin * sizeof(MergeMatrix));
                        incomplete = 0;
                        merge_required = 0;

                        for (i = 0; bin_start + i < bin_end; i++) {
                            family = bin_start[i].family;
                            if (mergematrix[family].first_match == 0)
                                mergematrix[family].first_match = bin_start[i].prevfamily;
                            else if (mergematrix[family].first_match == bin_start[i].prevfamily)
                                ;
                            else if (mergematrix[family].second_match == 0) {
                                mergematrix[family].second_match = bin_start[i].prevfamily;
                                if (mergematrix[mergematrix[family].first_match].swap == 0) {
                                    mergematrix[mergematrix[family].second_match].swap = mergematrix[family].first_match;
                                    merge_required = 1;
                                    }
                                else
                                    incomplete = 1;
                                }
                            else if (mergematrix[family].second_match == bin_start[i].prevfamily)
                                ;
                            else
                                incomplete = 1;
                            }

                        if (merge_required) {
                            //printf("merging\n");
                            for (i = 0; bin_start + i < bin_end; i++) {
                                family = bin_start[i].prevfamily;
                                if (mergematrix[family].swap != 0) {
                                    bin_start[i].prevfamily = mergematrix[family].swap;
                                    }
                                }
                            }
                            
                        } while (incomplete);

                    for (i = 0; bin_start + i < bin_end; i++) {
                        bin_start[i].family = 0;
                        }
                    } 
                }
            max_family = 0;
            for (i = 0; bin_start + i < bin_end; i++) {
                family = bin_start[i].prevfamily + current_family;
                bin_start[i].family = family;
                if (family > max_family)
                    max_family = family;
                }
            current_family = max_family;
            }
            
        //fprintf(stderr, "Fragment size = %i, bin size = %i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
        //printf("%i,%i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
        bin_start = bin_end;
        }
        
    free(mergematrix);
    return 0;
    }

    

static void print_family_sizes(PairedFastq *fq) {
    int counts[10] = { 0 }, members = 0, n = 0;
    ReadPair *bin_start = NULL, *bin_end = NULL;

    qsort(fq->readpairs, fq->total_reads, sizeof(ReadPair), (void *)compare_by_family);    
    n = 0;
    
    for (bin_start = fq->readpairs; bin_start < fq->readpairs + fq->total_reads;) {
        members = 0;
        for (bin_end = bin_start; bin_end < fq->readpairs + fq->total_reads; bin_end++) {
            if (bin_end->family != bin_start->family)
                break;
            members += bin_end->copy_number;
            }
        if (members < 10)
            counts[members]++;
        else
            counts[0]++;
        bin_start = bin_end;
        }
        
    for (n = 1; n < 10; n++) {
        if (counts[n] > 0) 
            fprintf(stderr, "%i x family size %i\n", counts[n], n);
        }
    if (counts[0] > 0)
        fprintf(stderr, "%i x family size 10+\n", counts[0]);
    
    }


    
static void print_family_sequences(PairedFastq *fq) {
    ReadPair *bin_start = NULL, *bin_end = NULL, *readpair = NULL, *readpair2 = NULL;
    
    qsort(fq->readpairs, fq->total_reads, sizeof(ReadPair), (void *)compare_by_family);    
    bin_start = fq->readpairs;
    readpair = fq->readpairs;
    while (bin_start < fq->readpairs + fq->total_reads) {
        for (bin_end = bin_start + 1; bin_end < fq->readpairs + fq->total_reads; bin_end++) {
            if (bin_end->family != bin_start->family)
                break;
            }
        if (bin_end->fragment_size == 0 && bin_end > bin_start + 1) {
            for (readpair2 = bin_start; readpair2 < bin_end; readpair2++)
                printf("%s %s\n", readpair2->read[R1].umi, readpair2->read[R2].umi);
            printf(".\n");
            }
        readpair++;
        bin_start = bin_end;
        }
    }


    
static int collapse_families(PairedFastq *fq) {
    ReadPair *bin_start = NULL, *bin_end = NULL, *readpair = NULL;
    int n = 0;
    
    qsort(fq->readpairs, fq->total_reads, sizeof(ReadPair), (void *)compare_by_family);    
    bin_start = fq->readpairs;
    readpair = fq->readpairs;
    while (bin_start < fq->readpairs + fq->total_reads) {
        for (bin_end = bin_start + 1; bin_end < fq->readpairs + fq->total_reads; bin_end++) {
            if (bin_end->family != bin_start->family)
                break;
            }
        if (bin_end > bin_start + 1) {
            merge_reads(bin_start, bin_end);
            n += bin_end - bin_start - 1;
            }
            
        if (bin_start > readpair)
            *readpair = *bin_start;
        readpair++;
        bin_start = bin_end;
        }
    fq->total_reads -= n;
    return n;
    }



static void merge_reads(ReadPair *bin_start, ReadPair *bin_end) {
    /* Merge all readpairs from bin_start+1 to bin_end-1 into bin_start.
     * 
     */
    int j = 0, i = 0, swap_int = 0, required = 0, tot = 0, counts[5] = {0, 0, 0, 0, INT_MAX}, read = R1;
    char *swap_charstar = NULL, bases[5] = {'A', 'T', 'C', 'G', 'N'};
    ReadPair *readpair = NULL;
    
    // Make sure that the first member of the family has the biggest seq and qual buffers.
    for (readpair = bin_start + 1; readpair < bin_end; readpair++) {
        for (read = R1; read <= R2; read++) {
            if (readpair->read[read].len > bin_start->read[read].len) {
                swap_charstar = bin_start->read[read].seq;
                bin_start->read[read].seq = readpair->read[read].seq;
                readpair->read[read].seq = swap_charstar;
                swap_charstar = bin_start->read[read].qual;
                bin_start->read[read].qual = readpair->read[read].qual;
                readpair->read[read].qual = swap_charstar;
                swap_int = bin_start->read[read].len;
                bin_start->read[read].len = readpair->read[read].len;
                readpair->read[read].len = swap_int;
                }
            }
        }
    
    for (read = R1; read <= R2; read++) {
        for (j = 0; j < bin_start->read[read].len; j++) {
            for (i = 0; i < 4; i++)
                counts[i] = 0;
            for (readpair = bin_start; readpair < bin_end; readpair++) {
                if (j < readpair->read[read].len) {
                    for (i = 0; i < 4; i++) {
                        if (readpair->read[read].seq[j] == bases[i]) {
                            counts[i] += readpair->copy_number;
                            break;
                            }
                        }
                    }
                }
            tot = counts[0] + counts[1] + counts[2] + counts[3];
            required = ((6 * tot) + 9) / 10;
            for (i = 0; i < 5; i++) {
                if (counts[i] >= required) {
                    bin_start->read[read].seq[j] = bases[i];
                    break;
                    }
                }
            
            if (bin_start->read[read].seq[j] == 'N')
                bin_start->read[read].qual[j] = '!';
            else {
                for (readpair = bin_start + 1; readpair < bin_end; readpair++) {
                    if (readpair->read[read].seq[j] == bin_start->read[read].seq[j] && readpair->read[read].qual[j] > bin_start->read[read].qual[j])
                        bin_start->read[read].qual[j] = readpair->read[read].qual[j];
                    }
                }
            }
        }
            
    for (readpair = bin_start + 1; readpair < bin_end; readpair++)
        bin_start->copy_number += readpair->copy_number;
    }



static int remove_unconfirmed_reads(PairedFastq *fq) {
    ReadPair *readpair = NULL, *readpair2 = NULL;
    int n = 0;
    
    readpair2 = fq->readpairs;
    for (readpair = fq->readpairs; readpair < fq->readpairs + fq->total_reads; readpair++) {
        if (readpair->fragment_size > 0 || readpair->copy_number > 1) {
            if (readpair > readpair2)
                *readpair2 = *readpair;
            readpair2++;
            }
        else
            n++;
        }
    fq->total_reads -= n;
    return n;
    }
    
    