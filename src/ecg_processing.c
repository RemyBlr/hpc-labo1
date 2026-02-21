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

/* ===============================================================================
 * Constantes
 * =============================================================================== */

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

/*
 * Fenêtre d'intégration en ms pour la MWI (Moving Window Integration).
 * Correspond à la durée typique du complexe QRS, soit environ 130 ms.
 * Pan-Thomkins précconise ~150ms.
 */
#define MWI_WINDOW_MS 130

/* ===============================================================================
 * Structures internes
 * =============================================================================== */

/*
 * Pré allocation des buffers pour éviter les alloc dynamiques pendant l'analyse.
 * Mémoire totale : 4 × MAX_SAMPLES × 8 bytes = 320 KB -> tient dans le cahce L2/L3.
 */
struct ECG_Context {
    // Copie locale des paramètres pour éviter les accès à la mémoire globale
    ECG_Params params;

    // Signal filtré (passe-bas)
    double *low_pass_buffer;
    // Signal dérivé (filtre passe-haut)
    double *high_pass_buffer;
    // Signal carré (non linéaire)
    double *squared_buffer;
    // Signal après fenêtre glissante (intégration)
    double *mwi_buffer;
};

/**
 * @brief Crée et init un contexte d'analyse ECG.
 *
 * Les appels à malloc sont couteux, on alloue donc les buffers une seule fois ici.
 * Ainsi, l'analyse se fait sans aucune alloc
 *
 * @param params Paramètres d'analyse (non NULL).
 * @return Pointeur vers le contexte en cas de succès, NULL sinon.
 */
ECG_Context *ecg_create(const ECG_Params *params) {
    if (!params) return NULL;

    ECG_Context *ctx = malloc(sizeof(ECG_Context));
    if (!ctx) return NULL;

    // Copie des paramètre pour éviter ptr externe qui pourrait être modifié.
    ctx->params = *params;

    // Alloc des buffers
    ctx->low_pass_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->high_pass_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->squared_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->mwi_buffer = malloc(sizeof(double) * MAX_SAMPLES);

    // Suffit d'un échec de malloc pour tout annuler
    if (!ctx->low_pass_buffer || !ctx->high_pass_buffer || !ctx->squared_buffer || !ctx->mwi_buffer) {
        free(ctx->low_pass_buffer);
        free(ctx->high_pass_buffer);
        free(ctx->squared_buffer);
        free(ctx->mwi_buffer);
        free(ctx);
        return NULL;
    }

    return ctx;
}

/**
 * @brief Libère un contexte d'analyse ECG et ses buffers internes.
 * @param ctx Le contexte à libérer. Les buffers internes sont également libérés.
 */
void ecg_destroy(ECG_Context *ctx) {
    if (!ctx) return;
    free(ctx->low_pass_buffer);
    free(ctx->high_pass_buffer);
    free(ctx->squared_buffer);
    free(ctx->mwi_buffer);
    free(ctx);
}

/* ===============================================================================
 * Fonctions utilitaires internes
 * =============================================================================== */

/**
 * Utilitaire pour affiner la position du pic R sur le signal NON filtré (car il y a un décalage).
 * La recherche se fait dans [center - half_window, center + half_window].
 *
 * @param signal Le signal ECG à analyser.
 * @param n_samples Le nombre d'échantillons dans le signal.
 * @param center L'indice central autour duquel chercher le maximum.
 * @param half_window La moitié de la fenêtre de recherche (en échantillons).
 * @return L'indice de l'échantillon avec la valeur maximale dans la fenêtre.
 */
static int find_max(const double *signal, size_t n_samples, int center, int half_window) {
    // On reste dans le tableau
    int start = (center - half_window > 0) ? center - half_window : 0;
    int end = (center + half_window < n_samples) ? center + half_window : (int)n_samples - 1;

    int best_id = start;
    double best_val = signal[start];

    // accès séquentiel
    for (int i = start + 1; i <= end; i++) {
        if (signal[i] > best_val) {
            best_val = signal[i];
            best_id = i;
        }
    }

    return best_id;
}

/* ===============================================================================
 * Fonction d'analyse principale
 * =============================================================================== */

ECG_Status ecg_analyze(ECG_Context *ctx,
                       const double *signal,
                       size_t n_samples,
                       int lead_idx,
                       ECG_Peaks *peaks,
                       ECG_Intervals *intervals) {
    // TODO treatment
    return ECG_OK;
}