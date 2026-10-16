/* Independent test instrumentation: counts allocation references linked into
 * the executable/static solver, excluding private calls inside shared DSOs. */
#include <stddef.h>
#include <stdint.h>
void *__real_malloc(size_t);
void *__real_calloc(size_t,size_t);
void *__real_realloc(void*,size_t);
static int enabled;
static int64_t calls, bytes;
void fsparse_refinement_alloc_begin(void) { calls=bytes=0;enabled=1; }
void fsparse_refinement_alloc_end(int64_t *count,int64_t *requested)
{ enabled=0;*count=calls;*requested=bytes; }
void *__wrap_malloc(size_t n)
{ if(enabled){++calls;bytes+=(int64_t)n;}return __real_malloc(n); }
void *__wrap_calloc(size_t n,size_t s)
{ if(enabled){++calls;bytes+=(int64_t)(n*s);}return __real_calloc(n,s); }
void *__wrap_realloc(void *p,size_t n)
{ if(enabled){++calls;bytes+=(int64_t)n;}return __real_realloc(p,n); }
