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
    char *name;
    char *seq;
    char *qual;
    char *seq2;
    char *qual2;
    char *umi;
    char *umi2;
    int len;
    int len2;
    int fragment_size;
    int family;
    int copy_number;
    } ReadPair;

typedef struct mergefamilystruct {
    int family;
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
static int read_fastqs(char *read1, char *read2, ReadPair **_readpairs, char **_contents);
static void debug_print_fastq(ReadPair *readpair);
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

    
    
static void debug_print_fastq(ReadPair *read) {
    printf("Name   = %s\n", read->name);
    printf("Seq    = %s\n", read->seq);
    printf("Seq2   = %s\n", read->seq2);
    printf("Qual   = %s\n", read->qual);
    printf("Qual2  = %s\n", read->qual2);
    printf("Umi1   = %s\n", read->umi);
    printf("Umi2   = %s\n", read->umi2);
    printf("Len    = %i\n", (int)(read->len));
    printf("Len2   = %i\n", (int)(read->len2));
    printf("Size   = %i\n", read->fragment_size);
    printf("Family = %i\n", read->family);
    printf("Copies = %i\n", read->copy_number);
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
                    read->qual[read->len] = '\0';
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

    ReadPair *readpairs = NULL, *bin_start = NULL, *bin_end = NULL, *subbin_start = NULL, *subbin_end = NULL, *readpair = NULL;
    MergeFamily *merge = NULL;
    
    char *read1 = NULL, *read2 = NULL, *contents = NULL;
    int total_reads = 0, i = 0, j = 0, k = 0, n = 0, matches = 0, current_family = 0, temp_family = 0, joined_family = 0, start_of_bin = 0, fragment_size = 0, last_fragment_size = -1;
    int unsized_start = 0, thruplex = 0, family = 0, incomplete = 0, merge_required = 0;
    
    int64_t v_umi = 0, last_v_umi = 0;
    
    // Parse the input tuple
    static char *kwlist[] = {"read1", "read2", "thruplex", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ss|p", kwlist, &read1, &read2, &thruplex))
        return NULL;

    if ((total_reads = read_fastqs(read1, read2, &readpairs, &contents)) == -1)
        goto error;
    
    merge = (MergeFamily *) calloc(total_reads, sizeof(MergeFamily));
    if (merge == NULL) {
        PyErr_NoMemory();
        goto error;
        }

     
    // DO WE NEED A MIN NUMBER OF VIABLE BASES READ?
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
    

    // Remove exact duplicates. Ns don't count in the comparison but the reads do have to be the same length. 
    // Therefore they may be duplicates even if they are not identical.
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_sequence);
    j = 0;
    n = 0;
    for (i = 1; i < total_reads; i++) {
        if (readpairs[i].len == readpairs[j].len && readpairs[i].len2 == readpairs[j].len2 && are_they_duplicates(&readpairs[i], &readpairs[j], 0)) {            
            readpairs[j].copy_number += 1;
            n++;
            
            // Merge qualities, keeping the highest and try and correct 'N's.
            for (k = 0; k < readpairs[j].len; k++) {
                if (readpairs[j].seq[k] == 'N')
                    readpairs[j].seq[k] = readpairs[i].seq[k];
                if (readpairs[j].qual[k] < readpairs[i].qual[k])
                    readpairs[j].qual[k] = readpairs[i].qual[k];
                }
            for (k = 0; k < readpairs[j].len2; k++) {
                if (readpairs[j].seq2[k] == 'N')
                    readpairs[j].seq2[k] = readpairs[i].seq2[k];
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

    
    
    // Fragment size paired reads, Correct sequencing errors if phred difference large enough or 'N' them out otherwise.
    // If thruplex then trim the end of the reads if there has been runthrough sequencing into the stem / umi at the other end.
    for (i = 0; i < total_reads; i++) {
        if ((readpairs[i].fragment_size = do_they_overlap(&readpairs[i], 70, 3, thruplex)) > 0)
            matches += 1;
        }
    printf("Sized = %i of %i (%i%\%)\n", matches, total_reads, matches * 100 / total_reads);
    
    
    
    // Trim thruplex umis
    if (thruplex) {
        for (i = 0; i < total_reads; i++) {
            readpairs[i].umi = readpairs[i].seq;
            readpairs[i].umi[THRUPLEX_UMT] = '\0';
            readpairs[i].seq += THRUPLEX_UMT + THRUPLEX_STEM;
            readpairs[i].qual += THRUPLEX_UMT + THRUPLEX_STEM;
            readpairs[i].len -= THRUPLEX_UMT + THRUPLEX_STEM;
            
            readpairs[i].umi2 = readpairs[i].seq2 + readpairs[i].len2 - THRUPLEX_UMT;
            readpairs[i].len2 -= THRUPLEX_UMT + THRUPLEX_STEM;
            readpairs[i].seq2[readpairs[i].len2] = '\0';
            readpairs[i].qual2[readpairs[i].len2] = '\0';
            
            if (readpairs[i].fragment_size > 0)
                readpairs[i].fragment_size -= 2 * (THRUPLEX_UMT + THRUPLEX_STEM);
            }
            printf("Trimmed thruplex umis and stems.\n");
        }
        
    
    // Group sized fragments into families
    //printf("Size,Number\n");
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_fragment_size_descending);
    bin_start = readpairs;
    while (bin_start < readpairs + total_reads) {
        for (bin_end = bin_start + 1; bin_end < readpairs + total_reads; bin_end++) {
            if (bin_end->fragment_size != bin_start->fragment_size)
                break;
            }
        
        if (bin_end - bin_start < 25000)
            assign_families(bin_start, bin_end, &current_family, 3);
        
        else {
            for (n = 0; n < 4; n++) {
                qsort(bin_start, bin_end - bin_start, sizeof(ReadPair), compare_by_fragment_size_descending); // NEED TO FIX
                subbin_start = bin_start;
                temp_family = 0;
                while (subbin_start < bin_end) {
                    for (subbin_end = subbin_start + 1; subbin_end < bin_end; subbin_end++) {
                        if (memcmp(subbin_start->seq, subbin_end->seq, 8) != 0) // NEED TO FIX
                            break;
                        }
                    assign_families(subbin_start, subbin_end, &temp_family, 3);
                    
                    subbin_start = subbin_end;
                    }
                
                if (n == 0) { // Don't need to merge families on first iteration.
                    for (i = 0; bin_start + i < bin_end; i++) {
                        merge[i].family = bin_start[i].family;
                        bin_start[i].family = 0;
                        }
                    }
                else {
                    incomplete = 0;
                    do {
                        merge_required = 0;
                        
                        
                        for (i = 0; bin_start + i < bin_end; i++) {
                            family = bin_start[i].family;
                            if (merge[family].first_match == 0)
                                merge[family].first_match = merge[i].family;
                            else if (merge[family].second_match == 0) {
                                merge[family].second_match = merge[i].family;
                                merge[merge[family].second_match].swap = merge[family].first_match;
                                merge_required = 1;
                                }
                            else
                                incomplete = 1;
                            }

                        if (merge_required) {
                            for (i = 0; bin_start + i < bin_end; i++) {
                                if (merge[bin_start[i].family].swap != 0)
                                    bin_start[i].family = merge[bin_start[i].family].swap;
                                }
                            }
                            
                        } while (incomplete);
                            
                            
                        
                        
                        merge[i].family = bin_start[i].family;
                        bin_start[i].family = 0;
                        }
                    }

                for (i = 0; bin_start + i < bin_end; i++) {
                    
                    
                } 
            }
                     
                     
                     
        printf("Fragment size = %i, bin size = %i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
        //printf("%i,%i\n", bin_start->fragment_size, (int)(bin_end - bin_start));
        bin_start = bin_end;
        }
    printf("Grouped into %i families.\n", current_family);
 
    
    int mismatches = 0;
    n = 0;
    qsort(readpairs, total_reads, sizeof(ReadPair), compare_by_family);
    bin_start = readpairs;
    while (bin_start < readpairs + total_reads) {
        for (bin_end = bin_start + 1; bin_end < readpairs + total_reads; bin_end++) {
            if (bin_end->family != bin_start->family)
                break;
            }
        if (bin_end - bin_start > 1) {
            n++;
            for (readpair = bin_start; readpair < bin_end; readpair++)
//                 j = 0;
//                 for (i = 0; i < 6; i++) {
//                     if (readpair->umi[i] != 0)
                        
                printf("umi = %s, umi2 = %s, copies = %i\n", readpair->umi, readpair->umi2, readpair->copy_number);
            printf(".\n.\n");
            }
            
        bin_start = bin_end;
        }
    
    
    
    
    
    
    
    
    
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
    free(merge);
    return ret;
    }
    
    
    
static int read_fastqs(char *read1, char *read2, ReadPair **_readpairs, char **_contents) { 
// Read paired fastq files and store the results in in ReadPair struct. 
// Strip newlines (\n or \r\n) from all strings.
// Ensure R1 sequence and quality are the same length and that R2 sequence and quality are the same length.
// Ensure that R1 and R2 names only differ where there is a '1' in the R1 name and a '2' in the R2 name. 
//
//
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
    
    printf("Filesize = %l MB, Reads = %l\n", filesize / 1024 / 1024, total_reads);

    
    // Allocate memory
    readpairs = (ReadPair *) malloc(sizeof(ReadPair) * total_reads);
    if (readpairs == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    *_readpairs = readpairs;
    contents = (char *) malloc(filesize + (total_reads * 5)); // (total_reads *5) for the '\0' at the end of name, seq, qual, seq2, qual2 .
    if (contents == NULL) {
        PyErr_NoMemory();
        goto cleanup;
        }
    *_contents = contents;
    
    
    // Read paired fasqs into memory and reverse complement read2
    for (i = 0; i < total_reads; i++) {
        // read 1 name
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        memcpy(contents, line, bytes_read + 1);
        readpairs[i].name = contents;
        contents += bytes_read + 1;
        
        // read 1 sequence
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        memcpy(contents, line, bytes_read + 1);
        readpairs[i].seq = contents;
        readpairs[i].len = bytes_read;
        contents += bytes_read + 1;
    
        // read 1 plus
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        // read 1 quality
        if ((bytes_read = getline(&line, &len, fp1)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        memcpy(contents, line, bytes_read + 1);
        readpairs[i].qual = contents;
        contents += bytes_read + 1;
        if (bytes_read != readpairs[i].len) {
            PyErr_SetString(PyExc_TypeError, "Sequence and quality differ in length.");
            goto cleanup;
            }
 
        
        // readv2 name
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
    
        // read 2 sequence
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        memcpy(contents, line, bytes_read + 1);
        readpairs[i].seq2 = contents;
        readpairs[i].len2 = bytes_read;
        reversecomplement(contents, readpairs[i].len2);
        contents += bytes_read + 1;
    
        // read 2 plus
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
     
        // read 2 quality
        if ((bytes_read = getline(&line, &len, fp2)) == -1) {
            PyErr_SetString(PyExc_TypeError, "Error reading file.");
            goto cleanup;
            }
        bytes_read = strip_newlines(line, bytes_read);
        memcpy(contents, line, bytes_read + 1);
        readpairs[i].qual2 = contents;
        reverse(contents, readpairs[i].len2);
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

    
    
    
    
    
    