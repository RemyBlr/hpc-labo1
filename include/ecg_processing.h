/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    ecg_processing.h
 * @brief   API d'analyse ECG : pré-traitement et détection de caractéristiques (P, Q, R, S, T).
 * @details Ce module définit une interface générique pour analyser un signal ECG.
 *          L'implémentation est libre (traitement offline, streaming, filtres au choix, etc.).

 */

#ifndef ECG_PROCESSING_H
#define ECG_PROCESSING_H

#include <stddef.h>
#include "output_structs.h"

/* ================================
 * Types
 * ================================ */

/**
 * @brief Contexte interne de traitement ECG (opaque).
 *
 * L'étudiant est libre de choisir les structures internes (buffers, filtres, états, etc.).
 */
typedef struct ECG_Context ECG_Context;


/**
 * @brief Paramètres d'analyse ECG.
 *
 * Cette structure permet de configurer l'analyse ECG avant son exécution.
 */
typedef struct
{
    int sampling_rate_hz;   /**< Fréquence d'échantillonnage (Hz). */
    int leads;              /**< Nombre de dérivations disponibles. */

    double gain; /**< Gain de l'amplification (optionnel). */

    /*ADD NECESSARY PARAMETERS*/

    double r_threshold_hint; /**< Seuil initial pour la détection des pics R (optionnel). */

} ECG_Params;


/**
 * @brief Code de retour standard pour l'analyse ECG.
 */
typedef enum
{
    ECG_OK = 0,          /**< Succès. */
    ECG_ERR_NULL = -1,   /**< Pointeur NULL. */
    ECG_ERR_PARAM = -2,  /**< Paramètres invalides. */
    ECG_ERR_ALLOC = -3,  /**< Échec d'allocation mémoire. */
    ECG_ERR_FAIL  = -4   /**< Erreur générique. */
} ECG_Status;

/* ================================
 * API Publique
 * ================================ */

 /**
 * @brief Crée un contexte d'analyse ECG.
 *
 * @param[in] params Paramètres d'analyse (non NULL).
 * @return Pointeur vers le contexte en cas de succès, NULL sinon.
 */
ECG_Context *ecg_create(const ECG_Params *params);

/**
 * @brief Libère un contexte d'analyse ECG.
 *
 * @param[in,out] ctx Contexte à libérer (peut être NULL).
 */
void ecg_destroy(ECG_Context *ctx);


/**
 * @brief Analyse un signal ECG et extrait les pics/caractéristiques.
 *
 * @details Fonction principale du module. Elle doit détecter les pics R les ondes P/Q/S/T et 
 *          les intervalles, peut-être fait ultérieurement.
 *
 * @param[in]  ctx        Contexte d'analyse (non NULL).
 * @param[in]  signal     Pointeur vers les échantillons ECG (non NULL).
 * @param[in]  n_samples  Nombre d'échantillons dans @p signal.
 * @param[in]  lead_idx   Index de la dérivation à analyser (0..LEADS-1).
 * @param[out] peaks      Résultats de détection des pics (non NULL).
 * @param[out] intervals  Intervalles calculés (peut être NULL si non demandé).
 *
 * @return ECG_OK en cas de succès, sinon un code d'erreur négatif.
 *
 * @note L'implémentation est libre (analyse globale ou incrémentale).
 * @note Les indices stockés dans @p peaks sont des indices d'échantillons (0..n_samples-1).
 */
ECG_Status ecg_analyze(
    ECG_Context *ctx,
    const double *signal,
    size_t n_samples,
    int lead_idx,
    ECG_Peaks *peaks,
    ECG_Intervals *intervals
);


/* AJOUTER N'IMPORTE QU'ELLE FONCTION UTILE */

#endif ECG_PROCESSING_H