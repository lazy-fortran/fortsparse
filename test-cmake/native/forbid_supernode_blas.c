/* Independent execution oracle: the strict CPU array budget must not reach
 * vendor routines whose hidden workspace has no accounting contract. */
#include <stdio.h>
#include <stdlib.h>
static void forbidden(const char *name)
{
    fprintf(stderr, "budgeted supernode execution reached %s\n", name);
    abort();
}
void __wrap_dgemm_(void) { forbidden("DGEMM"); }
void __wrap_dsyrk_(void) { forbidden("DSYRK"); }
void __wrap_dgemv_(void) { forbidden("DGEMV"); }
void __wrap_dtrsv_(void) { forbidden("DTRSV"); }
void __wrap_dtrsm_(void) { forbidden("DTRSM"); }
