#ifndef MY_VM_H_INCLUDED
#define MY_VM_H_INCLUDED
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

//Assume the address space is 32 bits, so the max memory size is 4GB
//Page size is 4KB

//Add any important includes here which you may need
#include <math.h>
#include <pthread.h>

#define PGSIZE 4096

// Maximum size of virtual memory
#define MAX_MEMSIZE 4ULL*1024*1024*1024

// Size of "physcial memory"
#define MEMSIZE 1024*1024*1024

// Represents a page table entry
typedef unsigned long pte_t;

// Represents a page directory entry
typedef unsigned long pde_t;

// Represents a virtual address
typedef unsigned long v_addr;

#define TLB_ENTRIES 512

/**
 * Structure to represent a tlb entry
 */
typedef struct tlb_entry{
    unsigned int virtual_page_number; // virtual page number
    unsigned int physical_page_number; // physical page number
    char valid; // is this tlb valid should be ascii 1 or 0
}tlb_entry;

/**
 * Structure to represents TLB
 */
typedef struct tlb {
    /*Assume your TLB is a direct mapped TLB with number of entries as TLB_ENTRIES
    * Think about the size of each TLB entry that performs virtual to physical
    * address translation.
    */
    tlb_entry entries[TLB_ENTRIES];
}tlb;


// Setup functions

void set_physical_mem();

// TLB Functions

int TLB_add(void *va, void *pa);
pte_t *TLB_check(void *va);
void print_TLB_missrate();

// Page Table Functions

void* translate(pde_t *pgdir, v_addr va);
int map_page(pde_t *pgdir, v_addr va, void** pa);

// Allocation functions

void *n_malloc(unsigned int num_bytes);
void n_free(void *va, int size);

// Data operations

int put_data(void *va, void *val, int size);
void get_data(void *va, void *val, int size);
void mat_mult(void *mat1, void *mat2, int size, void *answer);

#endif
