/* GPL helper protocol fixture: retain existing factor and normal operations,
 * reject only the new operation-capability query. Header layout is unchanged. */
#define main full_helper_main
#include "../../helper/fortsparse_umfpack_helper.c"
#undef main

int main(int argc, char **argv)
{
    helper_ctx c;
    factors_t f;
    fsparse_shm_header *h;
    if (argc < 5) return 2;
    memset(&c, 0, sizeof(c));
    memset(&f, 0, sizeof(f));
    if (attach(&c, argv, strtoll(argv[4], NULL, 10)) != 0) return 2;
    h = (fsparse_shm_header *) c.region;
#if !defined(_WIN32)
    post_done(&c);
#endif
    for (;;) {
        int done = 0;
        wait_req(&c);
        if (h->opcode == FSPARSE_OP_CAPABILITIES && h->factor_protocol == 3)
            h->status = FSPARSE_ST_ERROR;
        else
            done = dispatch(h, c.region, &f);
        post_done(&c);
        if (done) break;
    }
    free_factors(&f);
    detach(&c);
    return 0;
}
