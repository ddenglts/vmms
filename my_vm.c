#include "my_vm.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

struct tlb tlb_store;


void *get_next_avail(int num_pages); // get next set of contiguous pages in physical memory
int _clear_page_phys_bitmap(void *page_phys_addr);
int _set_page_phys_bitmap(void *page_phys_addr);
pte_t* _setup_page_table();
int _setup_page_directory(void);
pde_t *_find_first_free_pde();
pte_t *_find_first_free_pte(pte_t *page_table);
void * _create_page_and_link_to_page_table(pte_t *page_table);
pte_t* _create_page_table_and_link_to_dir(pde_t *pgdir_entry_location);
void _setupbitcounts();
static void free_program();
void _get_indices_from_va(v_addr va, unsigned long *pde_index, unsigned long *pte_index, unsigned long *offset);
v_addr _get_va_from_indicies( unsigned long pde_index, unsigned long pte_index, unsigned long offset);
v_addr _get_next_free_va(ulong size_in_pages);



//---------------------------------------------------bitmanip functions START---------------------------------------------------




/**
 * Function to set a bit at "index" of bitmap, zero-based numbering and little-endian
 * Each index is mapped to something and the bit at index states if it is valid
 * @param bitmap the bitmap to set
 * @param index the index of the bit to set
 */
static void set_bit(unsigned char *bitmap, ulong index){
    bitmap[index/8] = bitmap[index/8] | (1 << (index % 8));
}

/**
 * Function to get a bit at "index", zero-based numbering and little endian
 * Each index is mapped to something and the bit at index states if it is valid
 * @param bitmap the bitmap to get from
 * @param index the index of the bit to get
 */
static int get_bit(unsigned char *bitmap, ulong index){
    return (bitmap[index/8] >> (index % 8)) & 1;

}

/**
 * Function to clear a bit at "index", zero-based numbering and little endian
 * Each index is mapped to something and the bit at index states if it is valid
 * @param bitmap the bitmap to get from
 * @param index the index of the bit to get
 */
static void clear_bit(unsigned char *bitmap, ulong index){
    bitmap[index/8] = bitmap[index/8] & (~(1 << (index % 8)));
}

//---------------------------------------------------bitmanip functions END---------------------------------------------------

// Global Variables
unsigned int tlb_lookups = 0; // records the total amount of lookups for tlb
unsigned int tlb_misses = 0; // records the total amount of lookup misses for tlb

// Memory Management Variables
unsigned long long mem_size = MEMSIZE; // actual size of physical mem depending on what is smaller (MEMSIZE, MAX_MEMSIZE)
void* phys_mem = NULL; // physical memory pointer
pde_t* pde_start; // pointer to the start of the page directory entry (should only be one)
unsigned char *virt_bitmap = NULL; // virtual memory bitmap tells what pages have been validated for translation
unsigned char *phys_bitmap = NULL; // physical memory bitmap tells what pages have been allocated in physical memory
unsigned int total_pages; // total number of pages to keep track of given physical memory size and page size


// Bit Counts
unsigned int num_bits_offset; // stores the number of bits that offset takes
unsigned int num_bits_vpn; // stores the number of bits that physical page number takes (num_bits_pde + int num_bits_pte)
unsigned int num_bits_pde; // stores the number of bits that page directory entry takes
unsigned int num_bits_pte; // stores the number of bits for a page table entry takes

unsigned long VBITMAP_SIZE = 0; // measured in bits
unsigned long PBITMAP_SIZE = 0; // measured in bits


// All required functions in writeup ----------------------START-----------------------------

void set_physical_mem(){
    // setup bit counts
    _setupbitcounts();

    // allocate physical memory for main memory
    phys_mem = calloc(1, mem_size);
    if(phys_mem == NULL){
        fprintf(stderr, "Error: Unable to allocate physical memory.\n");
        exit(EXIT_FAILURE);
    }

    // dynamically allocate space for the bitmaps and set everything to 0
    ulong bytes_to_allocate = PBITMAP_SIZE/8;
    if (PBITMAP_SIZE % 8 != 0) {
        bytes_to_allocate++;
    }
    phys_bitmap = calloc(1, bytes_to_allocate); // 1 bit per page
    ulong bytes_to_allocate_v = VBITMAP_SIZE/8;
    if (VBITMAP_SIZE % 8 != 0) {
        bytes_to_allocate_v++;
    }
    virt_bitmap = calloc(1, bytes_to_allocate_v);
    if(phys_bitmap == NULL || virt_bitmap == NULL){
        fprintf(stderr, "Error: Unable to allocate physical memory.\n");
        free(phys_mem);
        exit(EXIT_FAILURE);
    }
    
    // edit bitmaps so that index 0 is used for storage of page directories
    // (physical memory + 0) is also the page directories page, makes it easier to manage
    pde_start = phys_mem;
    set_bit(phys_bitmap, 0);
    set_bit(virt_bitmap, 0);
    
    printf("Physical memory and bitmaps initialized successfully.\n");
    atexit(free_program);
}


/*
The function takes a virtual address and page directories starting address and
performs translation to return the physical address
*/
void *translate(pde_t *pgdir, v_addr va){

    unsigned long pde_index, pte_index, offset;
    _get_indices_from_va((v_addr) va, &pde_index, &pte_index, &offset);
    printf("translate(): va: %lu, pde_index: %lu, pte_index: %lu, offset: %lu\n", va, pde_index, pte_index, offset);

    pde_t *pde_entry = pgdir + pde_index;

    printf("translate(): pde_entry: %p, *: %p \n", pde_entry, (void*)*pde_entry);

    if ( (ulong) *pde_entry == 0){
        fprintf(stderr, "Error: Page directory entry not found in translate().\n");
        return NULL;
    }

    pte_t *pte_entry = (pte_t *) *pde_entry + pte_index;

    printf("translate(): pte_entry: %p , *: %p \n", pte_entry, (void*)*pte_entry);

    if ( (ulong) *pte_entry == 0){
        fprintf(stderr, "Error: Page table not found in translate().\n");
        return NULL;
    }

    return (void *) ((ulong) *pte_entry + offset);


    // check the tlb
    // TODO



    //If translation not successful, then return NULL
    return NULL; 
}



int map_page(pde_t *pgdir, v_addr va, void **pa){
    printf("map_page(): Mapping virtual address %lu to physical address %p\n", va, pa);

    assert(get_bit(virt_bitmap, (unsigned long) va >> num_bits_offset) == 0 && "map_page(): Virtual address already in use");

    unsigned long pde_index, pte_index, offset;
    _get_indices_from_va(va, &pde_index, &pte_index, &offset);
    printf("map_page(): va: %lu, pde_index: %lu, pte_index: %lu, offset: %lu\n", va, pde_index, pte_index, offset);

    // check if the page directory entry is valid
    pde_t *pde_entry = pgdir + pde_index;
    if ((ulong) *pde_entry == 0) {
        printf("map_page(): Creating new page table\n");
        // create a new page table and link it to the page directory
        pte_t *new_page_table = _create_page_table_and_link_to_dir(pde_entry);
        if (new_page_table == NULL) {
            fprintf(stderr, "man_page(): Error: Unable to create page table.\n");
            return -1;
        }
    }

    // check if the page table entry is valid
    pte_t *pte_entry_location = ((pte_t *) (*pde_entry)) + pte_index;
    if ((ulong) *pte_entry_location == 0) {
        printf("map_page(): Creating new page\n");
        // create a new page and link it to the page table
        if (_create_page_and_link_to_page_table(pte_entry_location) == NULL) {
            fprintf(stderr, "man_page(): Error: Unable to create page.\n");
            return -1;
        }
    }

    if (pa != NULL){
        *pa = (void *) *pte_entry_location;
    }

    return 0;

}


/*Function that gets the next available page 
*/
void *get_next_avail(int num_pages) {
    uint byte_index = 0;
    uint bit_index = 0;
    ulong current_seq = 0;

    // traverse the bitmap
    for (uint i = 0; i < PBITMAP_SIZE; ++i) {
        unsigned char byte = phys_bitmap[i];

        // skip fully occupied bytes
        if (byte == 0xFF) {
            current_seq = 0;
            continue;
        }

        for (int bit = 0; bit < 8; ++bit) {
            // check if the current bit is free
            if (get_bit(&byte, bit) == 0) {
                if (current_seq == 0) {
                    // start of a new sequence
                    byte_index = i;
                    bit_index = bit;
                }
                current_seq++;

                // if we found a valid sequence, calculate the address
                if (current_seq == num_pages) {
                    ulong start_bit = byte_index * 8 + bit_index;
                    printf("Found a sequence of %d pages starting at addr %p\n", num_pages, (void *) (((ulong) phys_mem) + start_bit * PGSIZE));
                    return (void *) (((ulong) phys_mem) + start_bit * PGSIZE);
                }
            } else {
                // reset the sequence if a used bit is found
                current_seq = 0;
            }
        }
    }

    // no suitable sequence found
    return NULL;
}


void *n_malloc(unsigned int num_bytes){


    if (phys_mem == NULL) {
        set_physical_mem();
    }



    if (pde == NULL) {
        _setup_page_directory();
    }


    // calculate the number of pages needed
    ulong num_pages = num_bytes / PGSIZE;
    if (num_bytes % PGSIZE != 0) { // round up if needed
        num_pages++;
    }

    ulong va = _get_next_free_va(num_pages); // get the next free virtual address from vbitmap ONLY
    if (va == 1) {
        fprintf(stderr, "Error: No available pages for allocation.\n");
        return NULL;
    }

    // if size is greater than 1 page, we need to map multiple pages
    for (int i = 0; i < num_pages; i++) {
        void *pa = NULL;
        if (map_page(pde, va + i * PGSIZE, pa) == -1) {
            fprintf(stderr, "Error: Unable to map page.\n");
            return NULL;
        }
    }

    // mark the virtual pages as used
    assert(va % PGSIZE == 0 && "Virtual address not a multiple of PGSIZE");
    for (int i = 0; i < num_pages; i++) {
        set_bit(virt_bitmap, (unsigned long) (va + i * PGSIZE) >> num_bits_offset);
    }


    return (void*) va;
}



void n_free(void *va, int size){

    // check if the virtual address is valid
    ulong size_in_pages = ceil((double) size / PGSIZE);
    // find all virtual pages in the range
    v_addr vpages[size_in_pages];
    for (int i = 0; i < size_in_pages; i++) {
        vpages[i] = (v_addr) va + i * PGSIZE;
        if (get_bit(virt_bitmap, (unsigned long) vpages[i] >> num_bits_offset) == 0) {
            fprintf(stderr, "Error: Virtual address not valid.\n");
            return;
        }
    }

    // keep track of unique page directory indices to check for empty page tables later
    unsigned long pde_indices[size_in_pages];
    int pde_indices_count = 0;

    for (int i = 0; i < size_in_pages; i++) {
        v_addr vpage = vpages[i];

        unsigned long pde_index, pte_index, offset;
        _get_indices_from_va(vpage, &pde_index, &pte_index, &offset);

        pde_t *pde_entry = pde + pde_index;

        if (*pde_entry == 0) {
            fprintf(stderr, "Error: Page directory entry not found in n_free().\n");
            return;
        }

        pte_t *pte_table = (pte_t *)(*pde_entry);
        pte_t *pte_entry = pte_table + pte_index;

        void *pa = (void *)*pte_entry;

        if (pa == NULL) {
            fprintf(stderr, "Error: Page table entry not found in n_free().\n");
            return;
        }

        // clear physical bitmap at physical address
        _clear_page_phys_bitmap(pa);
        printf("free: Cleared physical bitmap at %p\n", pa);

        // set page table entry to 0
        *pte_entry = 0;
        printf("free: Set page table entry to 0 at %p\n", pte_entry);

        // clear virtual bitmap at virtual page number
        clear_bit(virt_bitmap, (unsigned long)vpage >> num_bits_offset);
        printf("free: Cleared virtual bitmap at %lu\n", (unsigned long)vpage >> num_bits_offset);


        // add pde_index to the list if not already present
        int found = 0;
        for (int j = 0; j < pde_indices_count; j++) {
            if (pde_indices[j] == pde_index) {
                found = 1;
                break;
            }
        }
        if (!found) {
            pde_indices[pde_indices_count++] = pde_index;
        }
    }

    // check if page tables are empty and free them if necessary
    for (int i = 0; i < pde_indices_count; i++) {
        unsigned long pde_index = pde_indices[i];
        pde_t *pde_entry = pde + pde_index;

        if (*pde_entry == 0) {
            continue; // Already cleared
        }

        pte_t *pte_table = (pte_t *)(*pde_entry);
        int is_empty = 1;
        unsigned long num_pte_entries = 1UL << num_bits_pte;

        for (unsigned long j = 0; j < num_pte_entries; j++) {
            if (pte_table[j] != 0) {
                is_empty = 0;
                break;
            }
        }

        if (is_empty) {
            printf("free: Freeing page table at %p\n", pte_table);
            // free page table by clearing physical bitmap entries
            unsigned long pte_size_bytes = num_pte_entries * sizeof(pte_t);
            unsigned long pte_size_pages = pte_size_bytes / PGSIZE;
            if (pte_size_bytes % PGSIZE != 0) {
                pte_size_pages++;
            }

            for (unsigned long k = 0; k < pte_size_pages; k++) {
                void *addr = (void *)((ulong)pte_table + k * PGSIZE);
                _clear_page_phys_bitmap(addr);
                printf("free: Cleared physical bitmap at %p\n", addr);
            }

            // clear page directory entry
            *pde_entry = 0;
            printf("free: Set page directory entry to 0 at %p\n", pde_entry);
        }
    }
}



int put_data(void *va, void *val, int size){


    /*return -1 if put_data failed and 0 if put is successfull*/

     int bytes_written = 0;

    while (bytes_written < size) {
        // get the current virtual address
        v_addr current_va = ((unsigned long)va + bytes_written);

        // translate virtual address to physical address
        void *pa = translate(pde, current_va);
        if (pa == NULL) {
            fprintf(stderr, "Error: Invalid virtual address in put_data.\n");
            return -1;
        }

        // calculate bytes to copy in this iteration
        unsigned long offset_in_page = ((unsigned long)current_va) & ((1 << num_bits_offset) - 1);
        unsigned long bytes_left_in_page = PGSIZE - offset_in_page;
        unsigned long bytes_to_copy = size - bytes_written;
        if (bytes_to_copy > bytes_left_in_page) {
            bytes_to_copy = bytes_left_in_page;
        }

        // copy data from val to physical memory
        memcpy(pa, (void *)((unsigned long)val + bytes_written), bytes_to_copy);

        bytes_written += bytes_to_copy;
    }

    return 0;

}


void get_data(void *va, void *val, int size){


    int bytes_read = 0;

    while (bytes_read < size) {
        // get the current virtual address
        v_addr current_va = ((unsigned long)va + bytes_read);

        // translate virtual address to physical address
        void *pa = translate(pde, current_va);
        if (pa == NULL) {
            fprintf(stderr, "Error: Invalid virtual address in get_data.\n");
            return;
        }

        // calculate bytes to copy in this iteration
        unsigned long offset_in_page = ((unsigned long)current_va) & ((1 << num_bits_offset) - 1);
        unsigned long bytes_left_in_page = PGSIZE - offset_in_page;
        unsigned long bytes_to_copy = size - bytes_read;
        if (bytes_to_copy > bytes_left_in_page) {
            bytes_to_copy = bytes_left_in_page;
        }

        // copy data from physical memory to val
        memcpy((void *)((unsigned long)val + bytes_read), pa, bytes_to_copy);

        bytes_read += bytes_to_copy;
    }

}




void mat_mult(void *mat1, void *mat2, int size, void *answer){

    int x, y, val_size = sizeof(int);
    int i, j, k;
    for (i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            unsigned int a, b, c = 0;
            for (k = 0; k < size; k++){
                int address_a = (unsigned int)mat1 + ((i * size * sizeof(int))) + (k * sizeof(int));
                int address_b = (unsigned int)mat2 + ((k * size * sizeof(int))) + (j * sizeof(int));
                get_data( (void *)address_a, &a, sizeof(int));
                get_data( (void *)address_b, &b, sizeof(int));
                // printf("Values at the index: %d, %d, %d, %d, %d\n", 
                //     a, b, size, (i * size + k), (k * size + j));
                c += (a * b);
            }
            int address_c = (unsigned int)answer + ((i * size * sizeof(int))) + (j * sizeof(int));
            // printf("This is the c: %d, address: %x!\n", c, address_c);
            put_data((void *)address_c, (void *)&c, sizeof(int));
        }
    }
}


int TLB_add(void *va, void *pa){

    unsigned int vpn = (unsigned int) va >> num_bits_offset; // extract the virtaul page number
    unsigned int index = vpn % TLB_ENTRIES; // find the index of the tlb entry
    tlb_store.entries[index].virtual_page_number = vpn; // store vpn 
	tlb_store.entries[index].physical_page_number = (unsigned long) pa; // store physical frame
	tlb_store.entries[index].valid = 1; // sets tlb_entry as valid 
    return 0;
}



pte_t *TLB_check(void *va){

    unsigned int vpn = (unsigned int) va >> num_bits_offset; // extract the virtaul page number
    unsigned int index = vpn % TLB_ENTRIES; // find the index of the tlb entry

    // actual look up
    tlb_lookups++;
    // ensures that tlb is valid and the vpn matches
    if(tlb_store.entries[index].valid && tlb_store.entries[index].virtual_page_number == vpn){
        return (pte_t *) tlb_store.entries[index].physical_page_number;
    }

    // if dosen't return above then it is a miss
    tlb_misses++;
    /*This function should return a pte_t pointer*/
    return NULL;
}



void print_TLB_missrate(){
    double miss_rate = 0;	

    if(tlb_lookups == 0) miss_rate = 0.0; // makes sure we don't divide by zero
    
    miss_rate = (double)tlb_misses/tlb_lookups;
    fprintf(stderr, "TLB miss rate %lf \n", miss_rate);
}





static void free_program(){
    free(phys_mem);
    free(virt_bitmap);
    free(phys_bitmap);
}

/**
 * Function to setup the bit counts during setup
 */
void _setupbitcounts(){
    num_bits_offset = log2(PGSIZE);
    num_bits_pte = (32 - num_bits_offset) / 2;
    num_bits_pde = 32 - num_bits_pte - num_bits_offset;
    num_bits_vpn = num_bits_pte + num_bits_pde;
    total_pages = mem_size/PGSIZE;
    PBITMAP_SIZE = total_pages;
    VBITMAP_SIZE = MAX_MEMSIZE/PGSIZE;

    printf("num_bits_offset: %u\n", num_bits_offset);
    printf("num_bits_pte: %u\n", num_bits_pte);
    printf("num_bits_pde: %u\n", num_bits_pde);
    printf("num_bits_vpn: %u\n", num_bits_vpn);
    printf("total_pages: %u\n", total_pages);
    printf("PBITMAP_SIZE: %lu\n", PBITMAP_SIZE);
    printf("VBITMAP_SIZE: %lu\n", VBITMAP_SIZE);
}


// sets the value of the virtual address from the indices, notice POINTER va
v_addr _get_va_from_indicies(unsigned long pde_index, unsigned long pte_index, unsigned long offset){
    unsigned long vpn = (pde_index << num_bits_pte) + pte_index;
    return ((vpn << num_bits_offset) + offset);
}

void _get_indices_from_va(v_addr va, unsigned long *pde_index, unsigned long *pte_index, unsigned long *offset){
    unsigned long vpn = ((unsigned long) va) >> num_bits_offset;
    *pde_index = vpn >> num_bits_pte;
    *pte_index = vpn & ((1 << num_bits_pte) - 1);
    *offset = ((unsigned long) va) & ((1 << num_bits_offset) - 1);
}

// what is the VA?
// each bit represents a page
// it seems that VA will ALWAYS be a multiple of PGSIZE for the start of n_malloc
v_addr _get_next_free_va(ulong size_in_pages){
    uint byte_index = 0;
    uint bit_index = 0;
    ulong current_seq = 0;

    // traverse the bitmap
    for (uint i = 0; i < VBITMAP_SIZE; ++i) {
        unsigned char byte = virt_bitmap[i];

        // skip fully occupied bytes
        if (byte == 0xFF) {
            current_seq = 0;
            continue;
        }

        for (int bit = 0; bit < 8; ++bit) {
            // check if the current bit is free
            if (get_bit(&byte, bit) == 0) {
                if (current_seq == 0) {
                    // start of a new sequence
                    byte_index = i;
                    bit_index = bit;
                }
                current_seq++;

                // if we found a valid sequence, calculate the address
                if (current_seq == size_in_pages) {
                    ulong start_bit = byte_index * 8 + bit_index;
                    printf("Found a sequence of %lu virt pages starting at v addr %lx\n", current_seq, ( start_bit * PGSIZE));
                    return  (start_bit * PGSIZE);
                }
            } else {
                // reset the sequence if a used bit is found
                current_seq = 0;
            }
        }
    }

    // no suitable sequence found
    printf("No free virtual pages addresses.\n");
    return 1; // VAs should all be multiples of PGSIZE
}


// sets up the first level page dir (i.e. root directory) that will always be loaded
int _setup_page_directory(void){
    // CR3 >>> RECURSIVE ACCESS, we're using a agreed upon page dir address (pde)
    // how large is the page directory?
    unsigned long pde_size_bytes = pow(2, num_bits_pde) * sizeof(pde_t); // # of ptes per table, times size of pte
    unsigned long pde_size_pages = pde_size_bytes / PGSIZE;
    if (pde_size_bytes % PGSIZE != 0) { // so dir doesnt fit perfectly, add one more page for the rest
        pde_size_pages++;
    }
    printf("pde_size_bytes: %lu\n", pde_size_bytes);
    printf("pde_size_pages: %lu\n", pde_size_pages);

    // we set up the first level page directory at the start of the physical memory
    pde = (pde_t *) phys_mem;
    

    // mark p and v bitmaps; remember index 0 is for first page at phys_mem
    for (int i = 0; i < pde_size_pages; i++) {
        set_bit(phys_bitmap, i);
    }

    printf("Page directory initialized successfully.\n");
    return 0;


}

// sets up a page table and returns the address of the page table
// in detail, it finds the next available page for the page table based on bitmap
// marks the page in the bitmap
// links the page table to the page directory
// returns the address of the page table
pte_t* _create_page_table_and_link_to_dir(pde_t *pgdir_entry_location){
    // make sure the page directory is set up
    assert(pde != NULL && "Page directory not initialized");

    // how large is a page table?
    unsigned long pte_size_bytes = pow(2, num_bits_pte) * sizeof(pte_t);
    unsigned long pte_size_pages = pte_size_bytes / PGSIZE;
    if (pte_size_bytes % PGSIZE != 0) { // so dir doesnt fit perfectly, add one more page for the rest
        pte_size_pages++;
    }
    printf("pte_size_bytes: %lu\n", pte_size_bytes);
    printf("pte_size_pages: %lu\n", pte_size_pages);

    // find the next available page(s) for the page table
    void *first_page_in_seq = get_next_avail(pte_size_pages);
    if (first_page_in_seq == NULL) {
        fprintf(stderr, "Error: No available pages for page table.\n");
        return NULL;
    }
    
    
    // mark each page in the bitmap
    for (int i = 0; i < pte_size_pages; i++) {
        _set_page_phys_bitmap((void *)(((ulong) first_page_in_seq) + i * PGSIZE));
    }

    pde_t *free_pde = pgdir_entry_location;
    assert((ulong) *free_pde == 0 && "A Page table already linked to page directory");

    // link the page table to the page directory
    *free_pde = (pde_t) first_page_in_seq;

    return first_page_in_seq;

}




// creates a new page by finding the next available page in physical memory
// marks the page in the physical bitmap
// links the page to the page table specified by page_table_entry_location
// returns start address of the page
void * _create_page_and_link_to_page_table(pte_t *page_table_entry_location){
    // find the next available page
    void *page_phys_addr = get_next_avail(1);
    if (page_phys_addr == NULL) {
        fprintf(stderr, "Error: No available pages for page.\n");
        return NULL;
    }


    // mark the page in the bitmap
    _set_page_phys_bitmap(page_phys_addr);

    // link the page to the page table
    pte_t *free_pte = page_table_entry_location;
    assert((ulong) *free_pte == 0 && "A Page already linked to page table");
    printf("Setting page table entry to %p\n", page_phys_addr);

    *free_pte = (pte_t) page_phys_addr;

    return page_phys_addr;
}


int _clear_page_phys_bitmap(void *page_phys_addr){
    // get the page index
    unsigned long page_index = (unsigned long) ((unsigned long) page_phys_addr - (unsigned long) phys_mem)/PGSIZE;
    // clear the bit
    clear_bit(phys_bitmap, page_index);

    return 0;
}

int _set_page_phys_bitmap(void *page_phys_addr){
    // get the page index
    unsigned long page_index = (unsigned long) ((unsigned long) page_phys_addr - (unsigned long) phys_mem)/PGSIZE;
    // set the bit
    printf("Setting bit at index %lu\n", page_index);
    set_bit(phys_bitmap, page_index);

    return 0;
}



int main(){
    ulong va = (ulong) n_malloc(4097);
    printf("va: %lx\n", va);
    ulong va2 = (ulong) n_malloc(2048);
    printf("va2: %lx\n", va2);
    ulong va3 = (ulong) n_malloc(4096);
    printf("va3: %lx\n", va3);

    // print first 16 ulongs of pde, corresponding physical bitmap
    printf("pde location: %p\n", (void *) pde);
    for (int i = 0; i < 16; i++) {
        printf("%p pde[%d]: %p ---", pde + i, i, (void *) pde[i]);
        if (pde[i] != 0) {
            // print first 10 ulongs of the sub page table
            for (int j = 0; j < 10; j++) {
                printf(" %p ", (void *) ((pte_t *) pde[i])[j]);
            }
        }
        printf("\n");
    }
    for (int i = 0; i < 16; i++) {
        printf("%p phys_bitmap[%d]: %d, 0x%lx virt_bitmap[%d]: %d\n", phys_mem + i*PGSIZE, i, get_bit(phys_bitmap, i),((ulong) i * PGSIZE),i, get_bit(virt_bitmap, i));
    }

    printf("translate(pde, va): %p\n", translate(pde, va));
    printf("translate(pde, va2): %p\n", translate(pde, va2));
    printf("translate(pde, va3): %p\n", translate(pde, va3));

    // test put and get with a large array
    int size = 1024; // size of the array
    int *array = malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        array[i] = i;
    }

    put_data((void *) va, (void *) array, size * sizeof(int));

    int *array2 = malloc(size * sizeof(int));
    get_data((void *) va, (void *) array2, size * sizeof(int));

    for (int i = 0; i < size; i++) {
        if (array[i] != array2[i]) {
            printf("Mismatch at index %d: expected %d, got %d\n", i, array[i], array2[i]);
        }
        else{
            printf("Match at index %d: expected %d, got %d\n", i, array[i], array2[i]);
        }
    }

    free(array);
    free(array2);


    n_free((void *) va, 4097);
    n_free((void *) va2, 2048);
    n_free((void *) va3, 4096);

     // print first 16 ulongs of pde, corresponding physical bitmap
    printf("pde location: %p\n", (void *) pde);
    for (int i = 0; i < 16; i++) {
        printf("%p pde[%d]: %p ---", pde + i, i, (void *) pde[i]);
        if (pde[i] != 0) {
            // print first 10 ulongs of the sub page table
            for (int j = 0; j < 10; j++) {
                printf(" %p ", (void *) ((pte_t *) pde[i])[j]);
            }
        }
        printf("\n");
    }
    for (int i = 0; i < 16; i++) {
        printf("%p phys_bitmap[%d]: %d, 0x%lx virt_bitmap[%d]: %d\n", phys_mem + i*PGSIZE, i, get_bit(phys_bitmap, i),((ulong) i * PGSIZE),i, get_bit(virt_bitmap, i));
    }


    return 0;
}