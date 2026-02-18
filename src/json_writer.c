#include "json_writer.h"
#include <stdio.h>

static void write_int_array(FILE *f, const int *a, int n) {
    for (int i = 0; i < n; i++) {
        if (i) fputc(',', f), fputc(' ', f);
        fprintf(f, "%d", a[i]);
    }
}

static void write_double_array(FILE *f, const double *a, int n) {
    for (int i = 0; i < n; i++) {
        if (i) fputc(',', f), fputc(' ', f);
        fprintf(f, "%.2f", a[i]);
    }
}

int write_json(const char *filename, const ECG_Peaks *peaks, const ECG_Intervals *intervals) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("fopen");
        return -1;
    }

    static char outbuf[1 << 20];
    setvbuf(f, outbuf, _IOFBF, sizeof(outbuf));

    fputs("{\n", f);
    fputs("  \"peaks\": {\n", f);
    fputs("    \"R\": [", f);
    write_int_array(f, peaks->R, peaks->R_count);
    fputs("]\n  },\n", f);

    fputs("  \"intervals\": {\n", f);
    fputs("    \"RR\": [", f);
    write_double_array(f, intervals->RR, intervals->count);
    fputs("]\n  }\n", f);

    fputs("}\n", f);

    fclose(f);
    return 0;
}
