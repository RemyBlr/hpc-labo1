/**
 * @author: Bleuer Rémy
 * @file: ecg_processing.c
 * @brief: Implémentation de l'analyse ECG avvec une détéction de pics inpsiré par Pan-Tompkins.
 *
 * Voici les différentes optimisations mise en place afin d'améliorer les performances de l'analyse ECG :
 *
 * 1. Pas d'alloc dynmaique pendant l'analyse.
 * 2. Localité du cache.
 * 3. Streaming séquentiel.
 * 4. Capture des tendances avec seuils adaptatifs.
 * 5. Période réfractaire pour éviter les faux positifs.
 */
#include "ecg_processing.h"
#include "ecg_utils.h"

#include <stdlib.h>

/* ==================
 * Constantes
 * ================== */

/*
 * Période réfractaire minimale entre deux pics R.
 * Le coeur ne peut pas battre à plus de 220 bpm, soit environ 270 ms entre deux battements.
 * J'utilise donc 270 ms comme garde-fou pour éviter de détecter des pics R trop proches les uns des autres.
 */
#define REFRACTORY_PERIOD_MS 270
// Nombre d'échantillons qui correspond à la priode réfractaire
#define REFRACTORY_SAMPLES(heart_rate_hz) (heart_rate_hz * REFRACTORY_PERIOD_MS / 1000)

/*
 * Fenêtre du filtre passe-bas pour atténuer les hautes fréquences (bruit).
 * https://en.wikipedia.org/wiki/QRS_complex#:~:text=In%20adults,%20the%20QRS%20complex,wave%20follows%20the%20T%20wave.
 * La largeur du complexe QRS est généralement inférieure entre 70 et 110ms.
 * J'utilise donc 130 ms pour avoir tout le QRS et éviter de trop lisser les données.
 */
#define LOW_PASS_WINDOW_MS 130

/*
 * Seuil initial pour la détection des pics R.
 * Initilisation du seuil à 0.25 fois à l'amplitude max du signal (25% du max).
 */
#define THRESHOLD_INITIAL_FACTOR 0.25

/*
 * Mise à jour du seul adaptif après chaque pic R détecté.
 * J'utilise un facteur d'oubli expponentiel :
 * https://en.wikipedia.org/wiki/Moving_average
 * Valeurs définies par Pan-Tompkins :
 * signal_peak = 0.875 * signal_peak + 0.125 * new_peak
 * noise_peak = 0.875 * noise_peak + 0.125 * rejected_peak
 * Avantage : O(1) pour la màj du seuil, pas besoin de stocker les pics précédents.
 *
 */
#define SIGNAL_PEAK_DECAY_FACTOR 0.125
#define NOISE_PEAK_DECAY_FACTOR 0.125

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