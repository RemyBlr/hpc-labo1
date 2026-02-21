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

#include <stdio.h>

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
#define REFRACTORY_PERIOD_MS 200
// Nombre d'échantillons qui correspond à la priode réfractaire
#define REFRACTORY_SAMPLES(heart_rate_hz) ((REFRACTORY_PERIOD_MS * (heart_rate_hz)) / 1000)

/*
 * Fenêtre du filtre passe-bas pour atténuer les hautes fréquences (bruit).
 * https://en.wikipedia.org/wiki/QRS_complex#:~:text=In%20adults,%20the%20QRS%20complex,wave%20follows%20the%20T%20wave.
 * La largeur du complexe QRS est généralement inférieure entre 70 et 110ms.
 * J'utilise donc 130 ms pour avoir tout le QRS et éviter de trop lisser les données.
 */
#define LOW_PASS_WINDOW_MS 150

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
#define MWI_WINDOW_MS 150

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
    double *high_pass_buffer;
    // Signal dérivé (filtre passe-haut)
    double *derived_buffer;
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
    ctx->high_pass_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->derived_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->squared_buffer = malloc(sizeof(double) * MAX_SAMPLES);
    ctx->mwi_buffer = malloc(sizeof(double) * MAX_SAMPLES);

    // Suffit d'un échec de malloc pour tout annuler
    if (!ctx->high_pass_buffer || !ctx->derived_buffer || !ctx->squared_buffer || !ctx->mwi_buffer) {
        free(ctx->high_pass_buffer);
        free(ctx->derived_buffer);
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
    free(ctx->high_pass_buffer);
    free(ctx->derived_buffer);
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

/**
 * @brief Analyse un signal ECG et extrait les pics R et les intervalles RR.
 *
 * Pipeline d'analyse inspirée de Pan-Tompkins et optimisée pour la performance et la robustesse :
 *
 * 1. Filtre passe-haut pour atténuer les basses fréquences.
 * 2. Dérivation pour accentuer les transitions rapides (QRS).
 * 3. Mise au carré pour réctifier le signal et accentuer les pics.
 * 4. Intégration sur une fenêtre glissante pour lisser l'énergie du signal.
 * 5. Détectopm des pics R avec un seuil adaptif + période réfractaire.
 * 6. Affinage de la position des pics R sur le signal original.
 * 7. Calcul des intervalles RR à partir des indices des pics R détectés.
 *
 * @param ctx Contexte d'analyse (non NULL).
 * @param signal Pointeur vers les échantillons ECG (non NULL).
 * @param n_samples Nombre d'échantillons dans signal.
 * @param lead_idx Index de la dérivation à analyser (0..LEADS-1).
 * @param peaks Résultats de détection des pics (non NULL).
 * @param intervals Intervalles calculés (peut être NULL si non demandé).
 * @return ECG_OK en cas de succès, sinon un code d'erreur négatif.
 */
ECG_Status ecg_analyze(ECG_Context *ctx,
                       const double *signal,
                       size_t n_samples,
                       int lead_idx,
                       ECG_Peaks *peaks,
                       ECG_Intervals *intervals) {

    // Vérifications de base
    if (!ctx || !signal || !peaks) return ECG_ERR_NULL;
    if (n_samples == 0 || n_samples > MAX_SAMPLES) return ECG_ERR_PARAM;
    if (lead_idx < 0 || lead_idx >= ctx->params.leads) return ECG_ERR_PARAM;

    // fréquence en Hz
    const int fs = ctx->params.sampling_rate_hz;
    // Passe-haut
    double *hp = ctx->high_pass_buffer;
    // Dérivation
    double *deriv = ctx->derived_buffer;
    // Carré
    double *squared = ctx->squared_buffer;
    // Intégration
    double *mwi = ctx->mwi_buffer;

    //Calcul des fenêtres
    const size_t low_pass_window = (size_t)((LOW_PASS_WINDOW_MS * fs) / 1000);
    const size_t mwi_window = (size_t)((MWI_WINDOW_MS * fs) / 1000);
    const int refractory_samples = REFRACTORY_SAMPLES(fs);

    printf("[ECG] fs=%d Hz, low_pass_window=%zu samples, mwi_window=%zu samples, refractory_samples=%d samples\n",
           fs, low_pass_window, mwi_window, refractory_samples);

    // 1. Filtre passe-haut
    // Objectif : supprimer dérive lente de la ligne basse (< 1Hz)
    // Méthode : soustraction moyenne glissante (passe-haute =  x - MA(x))
    // HPC : O(n), zéro alloc
    ecg_highpass_ma(signal, hp, n_samples, low_pass_window);

    // 2. Dérivation discrète
    // Objectif : accentuer les transitions rapides (QRS a des pentes très raides par rapport aux ondes P et T)
    // Méthode : y[i] = x[i] - x[i-1], différence premier ordre pour simplicité et rapidité.
    // HPC : O(n), accès séquentiel, zéro alloc
    ecg_derivative_1(hp, deriv, n_samples);

    // 3. Mise au carré
    // Objectif : réctification, tout devient positif, accentuation non-linéaire des pics.
    // Méthode : y[i] = x[i]^2
    // HPC : O(n), multiplication simple par élément
    ecg_square(deriv, squared, n_samples);

    // 4. Intégration sur une fenêtre glissante (Moving Window Integration)
    // Objectif : lisser l'énergie du signal, faire ressortir les régions où le QRS est présent.
    // Méthode : moyenne glissante sur une fenêtre de taille mwi_window
    // HPC : O(n), somme glissante
    ecg_mwi(squared, mwi, n_samples, mwi_window);

    // 5. Détection des pics R avec seuil adaptatif + période réfractaire
    // Objectif : chercher les pic  R locaux qui dépassent un seuil dynamique.
    // Seuil adaptif:
    // - signal_peak : moyenne des vrais pics détectés.
    // - noise_peak : moyenne des candidats rejetés.
    // - threshold = noise_peak + 0.25 * (signal_peak - noise_peak)
    // HPC : O(1), pas de tableau d'historique à maintenir, pas de recalcul
    double max_mwi = 0.0;
    for (size_t i = 0; i < n_samples; i++)
        if (mwi[i] > max_mwi) max_mwi = mwi[i];

    double signal_peak = THRESHOLD_INITIAL_FACTOR * max_mwi;
    double noise_peak = THRESHOLD_INITIAL_FACTOR * max_mwi * 0.5;
    double threshold = noise_peak + 0.25 * (signal_peak - noise_peak);

    // Compteur de pics R détectés
    peaks->R_count = 0;

    // Indice dernier pic R détecté, init à -inf pour ne pas bloquer les premiers pics
    int last_r_index = -refractory_samples;

    // 6. Demi-fenêtre pour affinage local
    // On cherche le vrai max dans une fenetre de +-refinement_window autour du pic détecté sur le signal intégré
    const int refinement_window = refractory_samples / 2;

    // Boucle de détection
    // Pic local
    // HPC : O(n), accès séquentiel, on ne cible que le sommet
    for (size_t i = 1; i + 1 < n_samples && peaks->R_count < MAX_BEATS; i++)
    {
        // Pic local si (mwi[i] > mwi[i-1] et mwi[i] >= mwi[i+1])
        // Cette façon de faire permet de retourner que le 1er point d'un "plateau"
        int is_local_max = (mwi[i] > mwi[i-1]) && (mwi[i] >= mwi[i+1]);
        if (!is_local_max) continue;

        // Période refractaore
        if ((int)i - last_r_index < refractory_samples) {
            // Trop proche du dernier pic R détecté, on ignore pour éviter les faux positifs
            // O(1), pas de buffer
            noise_peak = (1.0 - NOISE_PEAK_DECAY_FACTOR) * noise_peak + NOISE_PEAK_DECAY_FACTOR * mwi[i];
            threshold = noise_peak + 0.25 * (signal_peak - noise_peak);
            continue;
        }

        // Test du seuil
        if (mwi[i] < threshold) {
            // Pas assez grand pour être un pic R, on considère comme du bruit
            noise_peak = (1.0 - NOISE_PEAK_DECAY_FACTOR) * noise_peak + NOISE_PEAK_DECAY_FACTOR * mwi[i];
            threshold = noise_peak + 0.25 * (signal_peak - noise_peak);
            continue;
        }

        // Pic R détecté
        // Pic dans mwi[] est décalé temporellement à cause du Moving Window Integration,
        // On affine en cherchant le vrai max local
        int r_index = find_max(signal, n_samples, (int)i, refinement_window);
        // TODO probablement changer ça
        //int r_index = (int)i - 13;
        //if (r_index < 0) r_index = 0;

        // Màj speak_noise avec le nouveau pic R détecté
        signal_peak = (1.0 - SIGNAL_PEAK_DECAY_FACTOR) * signal_peak + SIGNAL_PEAK_DECAY_FACTOR * mwi[i];
        threshold = noise_peak + 0.25 * (signal_peak - noise_peak);

        // Save du pic local
        peaks->R[peaks->R_count] = r_index;
        peaks->R_count++;

        // Indice dans mwi[] pour la période réfractaire
        last_r_index = (int)i;
    }

    printf("[ECG] Pics R détectés: %d\n", peaks->R_count);

    // 7. Calcul des intervalles RR
    // Objectif : calculer la durée entre chaque pic R détecté, en secondes
    // Méthode : RR[i] = (R[i+1] - R[i]) / fs
    // HPC : O(nbr pics), accès séquentiel, une seule passe
    if (intervals) {
        intervals->count = 0;

        // Pré calcul
        const double interval_fs = 1.0 / (double)fs;

        for (int i = 0; i + 1 < peaks->R_count && intervals->count < MAX_BEATS; i++) {
            int delta = peaks->R[i+1] - peaks->R[i];

            // Filtre des val aberrantes : on ignore les intervalles trop courts (< 200 ms) ou trop longs (> 2s)
            double rr = delta * interval_fs;
            if (rr >= 0.2 && rr <= 2.0) {
                intervals->RR[intervals->count] = rr;
                intervals->count++;
            }
        }

        printf("[ECG] Intervalles RR valides: %d\n", intervals->count);

        // Fréquence cardiaque moyenne en BPM
        if (intervals->count > 0) {
            double rr_sum = 0.0;
            for (int i = 0; i < intervals->count; i++)
                rr_sum += intervals->RR[i];

            double rr_mean = rr_sum / intervals->count;
            double bpm = 60.0 / rr_mean;
            printf("[ECG] Fréquence cardiaque moyenne: %.1f BPM\n", bpm);
            printf("[ECG] RR moyen: %.3f s\n", rr_mean);
        }
    }

    return ECG_OK;
}