#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>

int THRUPLEX_UMT = 6;
int THRUPLEX_STEM = 11;
int SIGNIFICANT_PHRED_DIFFERENCE = 10;



typedef struct readpairstruct {
    char *seq;
    char *seq2;
    char *umi;
    char *umi2;
    long name;
    long qual;
    long qual2;
    unsigned short len;
    unsigned short len2;
    unsigned short fragment_size;
    unsigned short copy_number;
    int family;
    } ReadPair;

    
    
typedef struct pairedfastqstruct {
    ReadPair *readpairs;
    int total_reads;
    char *contents;
    char *line;
    size_t line_len;
    FILE *fp;
    } PairedFastq;

    
    
typedef struct mergefamilystruct {
    int first_match;
    int second_match;
    int swap;
    } MergeFamily;

    
    
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



static void reversecomplement(char *start, int len);
static void reverse(char *start, int len);
static size_t strip_newlines(char *line, size_t length);
static int do_they_overlap(ReadPair *read, int min_overlap, int allowed, int thruplex);
static int compare_by_fragment_size_descending(const void *one, const void *two);
static int compare_by_family(const void *one, const void *two);
static int compare_by_sequence(const void *one, const void *two);
static int are_they_duplicates(ReadPair *read_a, ReadPair *read_b, int allowed);
static void print_sequences(ReadPair *read);
static PairedFastq read_fastqs(char *read1, char *read2);
static int debug_print_fastq(ReadPair *readpair, PairedFastq fq);
static void assign_families(ReadPair *bin_start, ReadPair *bin_end, int *current_family, int allowed);



static void assign_families(ReadPair *bin_start, ReadPair *bin_end, int *current_family, int allowed) {
    ReadPair *readpair, *readpair2, *readpair3;
    int joined_family;
    
    for (readpair = bin_start; readpair < bin_end; readpair++) {
        if (readpair->family == 0)
            readpair->family = ++(*current_family);
        
        for(readpair2 = readpair + 1; readpair2 < bin_end; readpair2++) {
            if (are_they_duplicates(readpair, readpair2, allowed)) {
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
static int do_they_overlap(ReadPair *read, int min_overlap, int allowed, int thruplex) {
    int read1_start, read_through, mismatches = 0;
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
            
            if (thruplex) {
                read_through = THRUPLEX_UMT + THRUPLEX_STEM - read1_start;
                if (read_through > 0) {
                    read->seq2 += read_through;
                    read->qual2 += read_through;
                    read->len2 -= read_through;
                    }
                read_through = THRUPLEX_UMT + THRUPLEX_STEM - (read->len2 - (read->len - read1_start));
                if (read_through > 0) {
                    read->len -= read_through;
                    read->seq[read->len] = '\0';
                    //read->qual[read->len] = '\0';
                    }
                }
            return read->len2 + read1_start;
            }
        }
    return 0;
    }


    
static PyObject *dedup(PyObject *self, PyObject *args, PyObject *kwargs) {
    //
    //
    
    PyObject *ret = NULL;

    ReadPair *readpairs = NULL, *bin_start = NULL, *bin_end = NULL, *subbin_start = NULL, *subbin_end = NULL, *readpair = NULL, *readpair2 = NULL;
    PairedFastq fq;
    
    MergeFamily *merge = NULL;
    int *merge_family = NULL, biggest_bin = 0, bin_size = 0;
    
    char *read1 = NULL, *read2 = NULL;
    int i = 0, j = 0, k = 0, n = 0, matches = 0, current_family = 0, temp_family = 0 ,thruplex = 0, family = 0, incomplete = 0, merge_required = 0, max_family = 0;
    
    
    // Parse the input tuple
    static char *kwlist[] = {"read1", "read2", "thruplex", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|p", kwlist, &read1, &read2, &thruplex))
        return NULL;

    fq = read_fastqs(read1, read2);
    if (fq.readpairs == NULL)
        goto error;

    
    // DO WE NEED A MIN NUMBER OF VIABLE BASES READ?
    // If the read is empty or consists entirely of 'N's then remove it as it is junk and will mess up our comparisons
    n = 0;
    readpair2 = fq.readpairs;
    for (readpair = fq.readpairs; readpair <  fq.readpairs + fq.total_reads; readpair++) {
        for (k = 0; k < readpair->len; k++) {
            if (readpair->seq[k] != 'N')
                break;
            }
        if (k == readpair->len) {
            n++;
            }
        else {
            readpair2++;
            if (readpair > readpair2)
                readpair2 = readpair;
            }
        }
    printf("Removed %i N only reads.\n", n);
    fq.total_reads -= n;

    
    
//     // Fragment size paired reads, Correct sequencing errors if phred difference large enough or 'N' them out otherwise.
//     // If thruplex then trim the end of the reads if there has been runthrough sequencing into the stem / umi at the other end.
//     for (i = 0; i < fq.total_reads; i++) {
//         if ((readpairs[i].fragment_size = do_they_overlap(&readpairs[i], 70, 3, thruplex)) > 0)
//             matches += 1;
//         }
//     printf("Sized = %i of %i (%i%\%)\n", matches, fq.total_reads, matches * 100 / fq.total_reads);
//     
//     
    
    // Trim thruplex umis
    if (thruplex) {
        for (readpair = fq.readpairs; readpair < fq.readpairs + fq.total_reads; readpair++) {
            readpair->umi = readpair->seq;
            readpair->umi[THRUPLEX_UMT] = '\0';
            readpair->seq += THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->qual += THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->len -= THRUPLEX_UMT + THRUPLEX_STEM;
            
            readpair->umi2 = readpair->seq2 + readpair->len2 - THRUPLEX_UMT;
            readpair->len2 -= THRUPLEX_UMT + THRUPLEX_STEM;
            readpair->seq2[readpair->len2] = '\0';
            //readpair->qual2[readpair->len2] = '\0'; // only needed if quality is not lazy loaded
            
            if (readpair->fragment_size > 0) // Can't subtract size from unsized pairs
                readpair->fragment_size -= 2 * (THRUPLEX_UMT + THRUPLEX_STEM);
            }
            printf("Trimmed thruplex umis and stems.\n");
        }
    
    fseek(fq.fp, fq.readpairs[0].qual, SEEK_SET);
    fputs("ZUB!!!!", fq.fp);
    debug_print_fastq(&fq.readpairs[0], fq);
    
//     // Group sized fragments into families
//     //printf("Size,Number\n");
//     qsort(readpairs, fq.total_reads, sizeof(ReadPair), compare_by_fragment_size_descending);
// 
// //     bin_start = readpairs;
// //     biggest_bin = 0;
// //     bin_size = 0;
// //     for (bin_end = bin_start + 1; bin_end < readpairs + fq.total_reads; bin_end++) {
// //         if (bin_end->fragment_size != bin_start->fragment_size) {
// //             bin_size = bin_end - bin_start;
// //             if (bin_size > biggest_bin) {
// //                 biggest_bin = bin_end - bin_start;
// //                 }
// //             }
// //         printf("bin size = %i\n");
// //         bin_start = bin_end;
// //         }
// //     printf("Largest size bin = %i\n", biggest_bin);    
// //     merge_family = (int *) calloc(fq.total_reads, (sizeof(MergeFamily) + sizeof(int)));
// //     if (merge_family == NULL) {
// //         PyErr_NoMemory();
// //         goto error;
// //         }
// //     merge = (ReadPair *)(merge_family + biggest_bin);
//     
//     goto finished;
//     
//     
//     
//     while (bin_start < readpairs + fq.total_reads) {
//         for (bin_end = bin_start + 1; bin_end < readpairs + fq.total_reads; bin_end++) {
//             if (bin_end->fragment_size != bin_start->fragment_size)
//                 break;
//             }
//         
//         if (bin_end - bin_start < 25000) {
//             assign_families(bin_start, bin_end, &current_family, 3);
//             }
//         else {
//             printf("woo hoo %i\n", (int)(bin_end - bin_start));
//             for (n = 0; n < 4; n++) {
//                 temp_family = 0;
//                 qsort(bin_start, (bin_end - bin_start) / sizeof(ReadPair), sizeof(ReadPair), compare_by_fragment_size_descending); // NEED TO FIX
//                 subbin_start = bin_start;
//                 while (subbin_start < bin_end) {
//                     for (subbin_end = subbin_start + 1; subbin_end < bin_end; subbin_end++) {
//                         if (memcmp(subbin_start->seq, subbin_end->seq, 8) != 0) // NEED TO FIX
//                             break;
//                         }
//                     assign_families(subbin_start, subbin_end, &temp_family, 3);
//                     
//                     subbin_start = subbin_end;
//                     }
//                 
//                 if (n == 0) { // Don't need to merge families on first iteration.
//                     for (i = 0; bin_start + i < bin_end; i++) {
//                         merge_family[i] = bin_start[i].family;
//                         bin_start[i].family = 0;
//                         }
//                     }
//                 else {
//                     memset(merge, 0, fq.total_reads * sizeof(ReadPair));
//                     
//                     do {
//                         incomplete = 0;
//                         merge_required = 0;
//                         
//                         for (i = 0; bin_start + i < bin_end; i++) {
//                             family = bin_start[i].family;
//                             if (merge[family].first_match == 0)
//                                 merge[family].first_match = merge_family[i];
//                             else if (merge[family].second_match == 0) {
//                                 merge[family].second_match = merge_family[i];
//                                 merge[merge[family].second_match].swap = merge[family].first_match;
//                                 merge_required = 1;
//                                 }
//                             else
//                                 incomplete = 1;
//                             }
// 
//                         if (merge_required) {
//                             for (i = 0; bin_start + i < bin_end; i++) {
//                                 family = merge_family[i];
//                                 if (merge[family].swap != 0)
//                                     merge_family[i] = merge[family].swap;
//                                 }
//                             }
//                             
//                         } while (incomplete);
// 
//                     for (i = 0; bin_start + readpairs[i].i < bin_end; i++)
//                         bin_start[i].family = 0;
//                     }
//                 } 
//             
//             max_family = 0;
//             for (i = 0; bin_start + 1 < bin_end; i++) {
//                 family = merge_family[i] + current_family;
//                 bin_start[i].family = family;
//                 if (family > max_family)
//                     max_family = family;
//                 }
//             current_family = max_family;
//             }
//         printf("Fragment size = %i, bin size = %i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
//         //printf("%i,%i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
//         bin_start = bin_end;
//         }
//     printf("Grouped into %i families.\n", current_family);
//  
//     
//     n = 0;
//     qsort(readpairs, fq.total_reads, sizeof(ReadPair), compare_by_family);
//     bin_start = readpairs;
//     while (bin_start < readpairs + fq.total_reads) {
//         for (bin_end = bin_start + 1; bin_end < readpairs + fq.total_reads; bin_end++)
//             if (bin_end->family != bin_start->family)
//                 break;
//             }
//         if (bin_end - bin_start > 1) {
//             n++;
//             for (readpair = bin_start; readpair < bin_end; readpair++) {
//                 ;
// //                 j = 0;
// //                 for (i = 0; i < 6; i++) {
// //                     if (readpair->umi[i] != 0)
//                         
// //                printf("umi = %s, umi2 = %s, copies = %i\n", readpair->umi, readpair->umi2, readpair->copy_number);
// //            printf(".\n.\n");
//             }
//             
//         bin_start = bin_end;
//         }
//     
//     
// 
//     
//     
//     
//     
//     
//     
//  /*   
//     n = 0;
//     family = 1;
//     for (i = 0; i < fq.total_reads; i++) {
//         if (readpairs[i].family == 0) // only needed if above function not allowed to complete
//             continue;
// 
//         if (readpairs[i].family == family)
//             n += 1;
//         else {
//             printf(".\n.\n");//family = %i, members = %i\n", family, n);
//             family = readpairs[i].family;
//             n = 1;
//             }
//         print_sequences(&readpairs[i]);
//         }*/
        
    finished:
    ret = Py_BuildValue("");
    error:
    free(fq.readpairs);
    free(fq.contents);
    free(fq.line);
    fclose(fq.fp);

    free(merge_family);
    return ret;
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
static PairedFastq read_fastqs(char *read1, char *read2) { 
// Read paired fastq files and store the results in in ReadPair struct. 
// Strip newlines (\n or \r\n) from all strings.
// Ensure R1 sequence and quality are the same length and that R2 sequence and quality are the same length.
// Ensure that R1 and R2 names only differ where there is a '1' in the R1 name and a '2' in the R2 name. 
//
//
    PairedFastq fq = {.readpairs = NULL, .total_reads = 0, .contents = NULL, .line = NULL, .line_len = 152, .fp = NULL};
    ReadPair *readpair;
    FILE *fp1, *fp2;
    char *contents = NULL, *line = NULL;
    int i = 0, j = 0, lines = 0, total_reads = 0;
    ssize_t bytes_read = 0;
    size_t line_len = 0;
    long buffer_len = 0, filesize = 0;
    
    
    // Open a new read/write file to act as a store for names and quality scores, we can't use the original fasqs as the quality score will be updated for all overlapping pairs.
    if (strlen(read1) > fq.line_len + 1)
        fq.line_len = strlen(read1) + 1;
    fq.line = calloc(fq.line_len, 1);
    if (fq.line == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    if (sprintf(fq.line, "%s.temp", read1) < 0) {
        PyErr_SetString(PyExc_TypeError, "Unable to sprintf file name!");
        goto cleanup;
        }
    fq.fp = fopen(fq.line, "w+");
    if (fq.fp == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unable to create temp buffer file.");
        goto cleanup;
        }
    
    
    // Open both fastqs, calculate file size and number of lines for memory allocation then reset file pointers to the start.
    fp1 = fopen(read1, "r");
    if (fp1 == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unable to open first file.");
        goto cleanup;
        }
    while ((bytes_read = getline(&fq.line, &fq.line_len, fp1)) != -1) {
        fq.total_reads += 1;
        if ((fq.total_reads - 2) % 4 == 0)
            buffer_len += bytes_read;
        }
    fq.total_reads = fq.total_reads / 4;
    filesize = ftell(fp1);
    fseek(fp1, 0L, SEEK_SET);
    
    fp2 = fopen(read2, "r");
    if (fp2 == NULL) {
        PyErr_SetString(PyExc_TypeError, "Unable to open second file.");
        goto cleanup;
        }
    while ((bytes_read = getline(&fq.line, &fq.line_len, fp2)) != -1) {
        total_reads += 1;
        if ((total_reads - 2) % 4 == 0)
            buffer_len += bytes_read;
        }
    total_reads = total_reads / 4;
    filesize += ftell(fp2);
    fseek(fp2, 0L, SEEK_SET);

    if (total_reads != fq.total_reads) {
        PyErr_SetString(PyExc_TypeError, "Read 1 and read 2 fastqs are of different length.");
        goto cleanup;
        }    
    printf("Filesize = %ld MB, Reads = %ld, Buffer = %d MB\n", filesize / 1024 / 1024, fq.total_reads, buffer_len / 1024 / 1024);

    
    // Allocate memory
    fq.readpairs = (ReadPair *) calloc(fq.total_reads, sizeof(ReadPair));
    if (fq.readpairs == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    fq.contents = (char *) malloc(buffer_len);
    if (fq.contents == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    contents = fq.contents;
    
    
    // Read paired fasqs into memory and reverse complement read2
    for (readpair = fq.readpairs; readpair < fq.readpairs + total_reads; readpair++) {
    
        // read 1 name
        readpair->name = ftell(fq.fp);
        if ((bytes_read = getline(&line, &line_len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        fputs(line, fq.fp);
//         strip_newlines(line, bytes_read);
//         memcpy(contents, fq.line, bytes_read + 1);
//         fq.readpairs[i].name = contents;
//         contents += bytes_read + 1;
        
        // read 1 sequence
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(fq.line, bytes_read);
        memcpy(contents, fq.line, bytes_read + 1);
        readpair->seq = contents;
        readpair->len = bytes_read;
        contents += bytes_read + 1;
    
        // read 1 plus
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        // read 1 quality
        readpair->qual = ftell(fq.fp);
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        fputs(fq.line, fq.fp);
        bytes_read = strip_newlines(fq.line, bytes_read);
//         memcpy(contents, fq.line, bytes_read + 1);
//         fq.readpairs[i].qual = contents;
//         contents += bytes_read + 1;
        if (bytes_read != readpair->len) {
            PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
            goto cleanup;
            }
 
        
        // read2 name
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(fq.line, bytes_read);
        for (j = 0; j < bytes_read; j++) {
            if (fq.line[j] != line[j] && (fq.line[j] != '2' || line[j] != '1')) {
                PyErr_SetString(PyExc_TypeError, "Reads 1 and 2 names dont match.");
                goto cleanup;
                }
            }
    
        // read 2 sequence
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(fq.line, bytes_read);
        memcpy(contents, fq.line, bytes_read + 1);
        readpair->seq2 = contents;
        readpair->len2 = bytes_read;
        reversecomplement(contents, readpair->len2);
        contents += bytes_read + 1;
    
        // read 2 plus
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        // read 2 quality
        readpair->qual2 = ftell(fq.fp);
        if ((bytes_read = getline(&fq.line, &fq.line_len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        fputs(fq.line, fq.fp);
        bytes_read = strip_newlines(fq.line, bytes_read);
//         memcpy(contents, fq.line, bytes_read + 1);
//         fq.readpairs[i].qual2 = contents;
//         reverse(contents, fq.readpairs[i].len2);
//         contents += bytes_read + 1;
        if (bytes_read != readpair->len2) {
            PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
            goto cleanup;
            }
            
        readpair->family = 0;
        readpair->copy_number = 1;
        }

        
    goto finished;
    cleanup:
    free(fq.readpairs);
    fq.readpairs = NULL;
    
    finished:
    close(fp1);
    close(fp2);
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

    

static int debug_print_fastq(ReadPair *readpair, PairedFastq fq) {
    int bytes_read = 0;
    
    fseek(fq.fp, readpair->name, SEEK_SET);
    if ((bytes_read = getline(&fq.line, &fq.line_len, fq.fp)) == -1) {
        PyErr_SetString(PyExc_TypeError, "Error reading file.");
        return -1;
        }
    strip_newlines(fq.line, bytes_read);
    printf("Name   = %s\n", fq.line);

    printf("Seq    = %s\n", readpair->seq);
    printf("Seq2   = %s\n", readpair->seq2);

    fseek(fq.fp, readpair->qual, SEEK_SET);
    if (fread(fq.line, sizeof(char), readpair->len, fq.fp) != readpair->len) {
        PyErr_SetString(PyExc_TypeError, "Error reading file.");
        return -1;
        }
    fq.line[readpair->len] = '\0';
    printf("Qual   = %s\n", fq.line);
 
    fseek(fq.fp, readpair->qual2, SEEK_SET);
    if (fread(fq.line, sizeof(char), readpair->len2, fq.fp) != readpair->len2) {
        PyErr_SetString(PyExc_TypeError, "Error reading file.");
        return -1;
        }
    fq.line[readpair->len2] = '\0';
    printf("Qual   = %s\n", fq.line);

    printf("Umi1   = %s\n", readpair->umi);
    printf("Umi2   = %s\n", readpair->umi2);
    printf("Len    = %i\n", (int)(readpair->len));
    printf("Len2   = %i\n", (int)(readpair->len2));
    printf("Size   = %i\n", readpair->fragment_size);
    printf("Family = %i\n", readpair->family);
    printf("Copies = %i\n", readpair->copy_number);
    
    return 0;
    }

    

// Remove exact duplicates. Ns don't count in the comparison but the reads do have to be the same length. 
// Therefore they may be duplicates even if they are not identical.
int remove_exact_duplicates(PairedFastq fq) {
    char *quality1 = NULL, *quality2 = NULL;    
    quality1 = malloc(fq.line_len + 1);
    
    
    qsort(fq.readpairs, fq.total_reads, sizeof(ReadPair), compare_by_sequence);
    n = 0;
    readpair2 = fq.readpairs;
    for (readpair = fq.readpairs; readpair < fq.readpairs + fq.total_reads; readpair++) {
        if (readpair->len == readpair2->len && readpair->len2 == readpair2->len2 && are_they_duplicates(readpair, readpair2, 0)) {            
            readpair2->copy_number += 1;
            n++;
            
            // Merge qualities, keeping the highest and try to correct 'N's.
            fseek(fq.fp, readpair->qual1, SEEK_SET);
            if (fread(quality1, sizeof(char), readpair->len, fq.fp) != readpair->len) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }
            fseek(fq.fp, readpair2->qual1, SEEK_SET);
            if (fread(quality2, sizeof(char), readpair2->len, fq.fp) != readpair2->len) {
                PyErr_SetString(PyExc_TypeError, "Error reading file.");
                goto cleanup;
                }                
            for (k = 0; k < readpair2->len; k++) {
                if (readpair2->seq[k] == 'N')
                    readpair2->seq[k] = readpair->seq[k];
                
//                 if (readpair2->qual[k] < readpair->qual[k])
//                     readpair2->qual[k] = readpair->qual[k];
                if (readpair2->qual[k] < readpair->qual[k])
                    readpair2->qual[k] = readpair->qual[k];
                
                }
            for (k = 0; k < readpair2->len2; k++) {
                if (readpair2->seq2[k] == 'N')
                    readpair2->seq2[k] = readpair->seq2[k];
                if (readpair2->qual2[k] < readpair->qual2[k])
                    readpair2->qual2[k] = readpair->qual2[k];
                }
             }
        else {
            j++;
            if (i > j)
                readpairs[j] = readpairs[i];
            }
        }
    printf("Exact duplicates = %i of %i (%i%%)\n", n, fq.total_reads, n * 100 / fq.total_reads);
    fq.total_reads -= n;
    return 0;
    }
    
    
    