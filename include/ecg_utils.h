/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    ecg_utils.h
 * @brief   Outils optionnels pour le pré-traitement ECG (filtres, dérivée, intégration).
 * @details Ces fonctions sont fournies comme aides. Leur usage est optionnel :
 *          les étudiants peuvent les utiliser, les modifier, ou implémenter leurs propres méthodes.
 *
*/

#ifndef ECG_UTILS_H
#define ECG_UTILS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Applique un gain au signal (in-place).
 *
 * @param[in,out] x  Signal.
 * @param[in]     n  Nombre d'échantillons.
 * @param[in]  gain  Gain multiplicatif.
 */
void ecg_apply_gain(double *x, size_t n, double gain);

/**
 * @brief Supprime une composante DC (offset) simple (in-place).
 *
 * @details Soustrait la moyenne du signal.
 *
 * @param[in,out] x Signal.
 * @param[in]     n Nombre d'échantillons.
 */
void ecg_remove_dc(double *x, size_t n);

/**
 * @brief Filtre moyenne glissante (low-pass simple).
 *
 * @param[in]  x     Signal d'entrée.
 * @param[out] y     Signal filtré.
 * @param[in]  n     Nombre d'échantillons.
 * @param[in]  win   Taille de la fenêtre (>= 1).
 */
void ecg_moving_average(const double *x, double *y, size_t n, size_t win);

/**
 * @brief Filtre passe-haut simple par soustraction de moyenne glissante.
 *
 * @details y = x - moving_average(x, win)
 *
 * @param[in]  x     Signal d'entrée.
 * @param[out] y     Signal filtré.
 * @param[in]  n     Nombre d'échantillons.
 * @param[in]  win   Taille de la fenêtre (>= 1).
 */
void ecg_highpass_ma(const double *x, double *y, size_t n, size_t win);

/**
 * @brief Dérivée discrète simple.
 *
 * @details y[i] = x[i] - x[i-1] (avec y[0]=0)
 *
 * @param[in]  x Signal d'entrée.
 * @param[out] y Signal dérivé.
 * @param[in]  n Nombre d'échantillons.
 */
void ecg_derivative_1(const double *x, double *y, size_t n);

/**
 * @brief Met au carré un signal (rectification énergie).
 *
 * @param[in]  x Signal d'entrée.
 * @param[out] y Signal de sortie.
 * @param[in]  n Nombre d'échantillons.
 */
void ecg_square(const double *x, double *y, size_t n);

/**
 * @brief Intégration sur fenêtre glissante (moving window integration).
 *
 * @param[in]  x     Signal d'entrée.
 * @param[out] y     Signal intégré.
 * @param[in]  n     Nombre d'échantillons.
 * @param[in]  win   Fenêtre d'intégration (>= 1).
 */
void ecg_mwi(const double *x, double *y, size_t n, size_t win);

#ifdef __cplusplus
}
#endif

#endif /* ECG_UTILS_H */
