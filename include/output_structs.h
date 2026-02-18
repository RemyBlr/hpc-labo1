/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    output_structs.h
 * @brief   Définitions des structures et constantes de sortie pour l’analyse ECG.
 * @details Ce fichier regroupe les constantes globales ainsi que les structures
 *          utilisées pour stocker les résultats de détection des pics ECG
 *          et le calcul des intervalles temporels.
 *
 */

#ifndef OUTPUT_STRUCTS_H
#define OUTPUT_STRUCTS_H

/* ================================
 * Constantes globales
 * ================================ */

/** @brief Nombre maximal d’échantillons ECG traités. */
#define MAX_SAMPLES     10000

/** @brief Nombre de dérivations ECG. */
#define LEADS           12

/** @brief Taille maximale de l’historique RR. */
#define MAX_RR_HISTORY  8

/** @brief Fréquence d’échantillonnage ECG (en Hz). */
#define SAMPLING_RATE   500

/** @brief Fréquence cardiaque maximale théorique (en BPM). */
#define HR_MAX_BPM      240

/** @brief Durée maximale des données ECG (en secondes). */
#define DURATION_S      (MAX_SAMPLES / SAMPLING_RATE)

/** @brief Nombre maximal de battements cardiaques attendus. */
#define MAX_BEATS       ((HR_MAX_BPM * DURATION_S) / 60 + 16)

/* ================================
 * Structures de données
 * ================================ */
/**
 * @brief Structure contenant les indices des pics ECG détectés.
 *
 * Chaque tableau stocke les indices (échantillons) correspondant aux
 * différents types de pics (P, Q, R, S, T).
 */
typedef struct {
    int R[MAX_BEATS]; /**< Indices des pics R détectés. */
    int P[MAX_BEATS]; /**< Indices des pics P détectés. */
    int Q[MAX_BEATS]; /**< Indices des pics Q détectés. */
    int S[MAX_BEATS]; /**< Indices des pics S détectés. */
    int T[MAX_BEATS]; /**< Indices des pics T détectés. */

    int R_count; /**< Nombre de pics R détectés. */
    int P_count; /**< Nombre de pics P détectés. */
    int Q_count; /**< Nombre de pics Q détectés. */
    int S_count; /**< Nombre de pics S détectés. */
    int T_count; /**< Nombre de pics T détectés. */
} ECG_Peaks;

/**
 * @brief Structure contenant les intervalles RR.
 *
 * Les intervalles RR représentent le temps entre deux pics R successifs,
 * généralement exprimé en secondes ou en échantillons selon l’implémentation.
 */
typedef struct {
    double RR[MAX_BEATS];   /**< Valeurs des intervalles RR. */
    int    count;           /**< Nombre d’intervalles RR valides. */
} ECG_Intervals;

#endif /* OUTPUT_STRUCTS_H */
