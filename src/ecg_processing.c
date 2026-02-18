#include "ecg_processing.h"
#include <stdlib.h>

struct ECG_Context {
    // TODO
    int dummy;
};

ECG_Context *ecg_create(const ECG_Params *params) {
    // minimal alloc for testing prupose
    ECG_Context *ctx = malloc(sizeof(ECG_Context));
    if (ctx) ctx->dummy = 0;
    return ctx;
}

void ecg_destroy(ECG_Context *ctx) {
    free(ctx);
}

ECG_Status ecg_analyze(ECG_Context *ctx,
                       const double *signal,
                       size_t n_samples,
                       int lead_idx,
                       ECG_Peaks *peaks,
                       ECG_Intervals *intervals) {
    // TODO treatment
    return ECG_OK;
}