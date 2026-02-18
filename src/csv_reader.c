#include "csv_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>

double ecg_data[LEADS][MAX_SAMPLES];
int sample_count = 0;

static void skip_spaces(char **p) {
    while (**p && isspace((unsigned char)**p)) (*p)++;
}

int read_csv(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("fopen");
        return -1;
    }

    static char filebuf[1 << 20]; // 1 MiB
    setvbuf(f, filebuf, _IOFBF, sizeof(filebuf));

    int c;
    while ((c = fgetc(f)) != '\n' && c != EOF) { /* skip */ }

    char *line = NULL;
    size_t cap = 0;

    int lead = 0;
    int loaded_samples = -1;

    while (lead < LEADS) {
        ssize_t n = getline(&line, &cap, f);
        if (n == -1) break;

        char *p = line;

        while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) {
            line[--n] = '\0';
        }

        while (*p && *p != ',') p++;
        if (*p == ',') p++;

        int s = 0;
        while (*p && s < MAX_SAMPLES) {
            skip_spaces(&p);
            if (!*p) break;

            errno = 0;
            char *end = NULL;
            double v = strtod(p, &end);

            if (end == p) {
                while (*p && *p != ',') p++;
                if (*p == ',') p++;
                continue;
            }
            if (errno == ERANGE) {
                // overflow/underflow, on garde v mais tu peux décider de return error
            }

            ecg_data[lead][s++] = v;
            p = end;

            while (*p && *p != ',') p++;
            if (*p == ',') p++;
        }

        if (loaded_samples == -1) loaded_samples = s;

        lead++;
    }

    free(line);
    fclose(f);

    if (lead == 0 || loaded_samples <= 0) {
        fprintf(stderr, "Erreur: aucun lead/échantillon lu.\n");
        return -2;
    }

    sample_count = loaded_samples;
    printf("CSV chargé avec %d leads et %d échantillons.\n", lead, sample_count);
    return 0;
}
