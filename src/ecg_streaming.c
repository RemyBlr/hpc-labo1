// ecg_streaming.c — version corrigée
#include "ecg_streaming.h"
#include "ecg_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHUNK_SIZE    1000
#define OVERLAP       250
#define STEP          (CHUNK_SIZE - OVERLAP)   // 750

typedef struct {
    double hp_sum;
    size_t hp_w;
    double mwi_sum;
    size_t mwi_w;
    double signal_peak;
    double noise_peak;
    double threshold;
    int    last_r_global;
    long   global_offset;
    double hp[CHUNK_SIZE];
    double deriv[CHUNK_SIZE];
    double squared[CHUNK_SIZE];
    double mwi[CHUNK_SIZE];
} StreamState;

static int find_max_local(const double *sig, int n, int center, int half_win) {
    int start = center - half_win; if (start < 0) start = 0;
    int end   = center + half_win; if (end >= n)  end   = n - 1;
    int best  = start;
    for (int i = start + 1; i <= end; i++)
        if (sig[i] > sig[best]) best = i;
    return best;
}

/*
 * Calcule le vrai max_mwi sur le premier paquet en faisant tourner le pipeline
 * complet (sans état persistant, juste pour calibrer le seuil initial).
 * C'est exactement ce que fait ecg_analyze() avec son `max_mwi`.
 */
static double calibrate_threshold(const double *signal, int n, int fs) {
    const size_t hp_win  = (size_t)((130 * fs) / 1000);
    const size_t mwi_win = (size_t)((130 * fs) / 1000);

    static double hp[CHUNK_SIZE], deriv[CHUNK_SIZE],
                  squared[CHUNK_SIZE], mwi_buf[CHUNK_SIZE];

    ecg_highpass_ma(signal, hp, n, hp_win);
    ecg_derivative_1(hp, deriv, n);
    ecg_square(deriv, squared, n);
    ecg_mwi(squared, mwi_buf, n, mwi_win);

    double max_mwi = 0.0;
    for (int i = 0; i < n; i++)
        if (mwi_buf[i] > max_mwi) max_mwi = mwi_buf[i];

    return max_mwi;
}

static void stream_state_init(StreamState *s,
                               const double *signal, int n_samples, int fs) {
    memset(s, 0, sizeof(*s));

    int init_n = (n_samples < CHUNK_SIZE) ? n_samples : CHUNK_SIZE;
    double max_mwi = calibrate_threshold(signal, init_n, fs);

    printf("[STREAMING] max_mwi calibration = %f\n", max_mwi);

    s->signal_peak   = 0.25 * max_mwi;
    s->noise_peak    = 0.25 * max_mwi * 0.5;
    s->threshold     = s->noise_peak + 0.25 * (s->signal_peak - s->noise_peak);
    s->last_r_global = -9999;

    printf("[STREAMING] seuil initial = %f, signal_peak = %f\n",
           s->threshold, s->signal_peak);
}

static int process_chunk(StreamState *s,
                         const double *chunk, int n,
                         const ECG_Params *params,
                         int *all_peaks, int max_peaks, int peak_count)
{
    const int    fs       = params->sampling_rate_hz;
    const size_t hp_win   = (size_t)((130 * fs) / 1000);
    const size_t mwi_win  = (size_t)((130 * fs) / 1000);
    const int    refract  = (270 * fs) / 1000;
    const int    ref_win  = refract / 2;

    /* 1. Passe-haut avec état persistant */
    for (int i = 0; i < n; i++) {
        s->hp_sum += chunk[i];
        s->hp_w++;
        if (s->hp_w > hp_win) {
            s->hp_sum -= chunk[i - (int)hp_win];
            s->hp_w--;
        }
        double ma = (s->hp_w > 0) ? s->hp_sum / (double)s->hp_w : 0.0;
        s->hp[i] = chunk[i] - ma;
    }

    /* 2. Dérivée */
    s->deriv[0] = 0.0;
    for (int i = 1; i < n; i++)
        s->deriv[i] = s->hp[i] - s->hp[i-1];

    /* 3. Carré */
    for (int i = 0; i < n; i++)
        s->squared[i] = s->deriv[i] * s->deriv[i];

    /* 4. MWI avec état persistant */
    for (int i = 0; i < n; i++) {
        s->mwi_sum += s->squared[i];
        s->mwi_w++;
        if (s->mwi_w > mwi_win) {
            s->mwi_sum -= s->squared[i - (int)mwi_win];
            s->mwi_w--;
        }
        s->mwi[i] = (s->mwi_w > 0) ? s->mwi_sum / (double)s->mwi_w : 0.0;
    }

    /* 5. Détection — zone acceptée = après OVERLAP (sauf 1er paquet) */
    int accepted_from = (s->global_offset == 0) ? 0 : OVERLAP;

    for (int i = 1; i + 1 < n && peak_count < max_peaks; i++) {
        int is_local_max = (s->mwi[i] > s->mwi[i-1])
                        && (s->mwi[i] >= s->mwi[i+1]);
        if (!is_local_max) continue;

        int global_i = (int)s->global_offset + i;

        if (global_i - s->last_r_global < refract) {
            s->noise_peak = 0.875 * s->noise_peak + 0.125 * s->mwi[i];
            s->threshold  = s->noise_peak + 0.25*(s->signal_peak - s->noise_peak);
            continue;
        }

        if (s->mwi[i] < s->threshold) {
            s->noise_peak = 0.875 * s->noise_peak + 0.125 * s->mwi[i];
            s->threshold  = s->noise_peak + 0.25*(s->signal_peak - s->noise_peak);
            continue;
        }

        /* Pic accepté */
        s->signal_peak   = 0.875 * s->signal_peak + 0.125 * s->mwi[i];
        s->threshold     = s->noise_peak + 0.25*(s->signal_peak - s->noise_peak);
        s->last_r_global = global_i;

        if (i >= accepted_from) {
            int refined = find_max_local(chunk, n, i, ref_win);
            all_peaks[peak_count++] = (int)s->global_offset + refined;
        }
    }

    return peak_count;
}

int ecg_analyze_streaming(const ECG_Params *params,
                          const double *signal,
                          int n_samples,
                          ECG_Peaks *peaks,
                          ECG_Intervals *intervals)
{
    if (!params || !signal || !peaks) return -1;

    StreamState state;
    stream_state_init(&state, signal, n_samples, params->sampling_rate_hz);

    static int all_peaks[MAX_BEATS];
    int peak_count = 0;
    int pos = 0;

    while (pos < n_samples) {
        int n = n_samples - pos;
        if (n > CHUNK_SIZE) n = CHUNK_SIZE;

        peak_count = process_chunk(&state, signal + pos, n,
                                   params, all_peaks, MAX_BEATS, peak_count);

        if (n < CHUNK_SIZE) break;   /* dernier paquet partiel, on a fini */
        pos += STEP;
        state.global_offset += STEP;
    }

    peaks->R_count = peak_count;
    for (int i = 0; i < peak_count; i++)
        peaks->R[i] = all_peaks[i];

    printf("[STREAMING] %d pics R détectés\n", peak_count);

    if (intervals) {
        intervals->count = 0;
        const double inv_fs = 1.0 / (double)params->sampling_rate_hz;
        for (int i = 0; i+1 < peak_count && intervals->count < MAX_BEATS; i++) {
            double rr = (peaks->R[i+1] - peaks->R[i]) * inv_fs;
            if (rr >= 0.2 && rr <= 2.0)
                intervals->RR[intervals->count++] = rr;
        }
        if (intervals->count > 0) {
            double sum = 0;
            for (int i = 0; i < intervals->count; i++) sum += intervals->RR[i];
            printf("[STREAMING] BPM moyen: %.1f\n",
                   60.0 / (sum / intervals->count));
        }
    }

    return 0;
}