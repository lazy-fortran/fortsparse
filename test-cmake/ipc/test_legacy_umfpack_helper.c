/* Protocol adversary: reproduce the old helper's unknown-opcode reply.
 * Reuse only its transport; record every factor request independently. */
#define main current_umfpack_main
#include "../../helper/fortsparse_umfpack_helper.c"
#undef main

int main(int argc, char **argv)
{
    helper_ctx c;
    fsparse_shm_header *h;
    int factor_requests = 0;
    const char *receipt = getenv("FORTSPARSE_LEGACY_RECEIPT");
    FILE *out;
    if (argc < 5 || receipt == NULL) return 2;
    memset(&c, 0, sizeof(c));
    if (attach(&c, argv, strtoll(argv[4], NULL, 10)) != 0) return 2;
    h = (fsparse_shm_header *) c.region;
#if !defined(_WIN32)
    post_done(&c);
#endif
    for (;;) {
        wait_req(&c);
        if (h->opcode == FSPARSE_OP_SHUTDOWN) break;
        if (h->opcode == FSPARSE_OP_FACTOR_REAL ||
            h->opcode == FSPARSE_OP_FACTOR_COMPLEX) ++factor_requests;
        /* All unsupported requests, including CAPS, complete without work. */
        h->status = FSPARSE_ST_ERROR;
        post_done(&c);
    }
    out = fopen(receipt, "w");
    if (out == NULL) return 2;
    fprintf(out, "%d\n", factor_requests);
    fclose(out);
    detach(&c);
    return 0;
}
