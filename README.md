# vmms

The `vmms` project is a virtual memory management system using pages. It provides a set of library calls for memory management and data manipulation:

## Memory allocation and de-allocation calls

- `void *n_malloc(unsigned int num_bytes);`: Allocates `num_bytes` of memory and returns a pointer to the allocated memory.
- `void n_free(void *va, int size);`: Frees the memory block pointed to by `va` of the specified `size`.

## Data operation calls

- `int put_data(void *va, void *val, int size);`: Stores `size` bytes of data from `val` into the memory location pointed to by `va`.
- `void get_data(void *va, void *val, int size);`: Retrieves `size` bytes of data from the memory location pointed to by `va` and stores it in `val`.