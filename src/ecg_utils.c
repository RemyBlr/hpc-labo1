/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    ecg_utils.c
 * @brief   Implémentation d'outils optionnels pour le pré-traitement ECG.
 * @details Ces fonctions sont fournies comme aides (optionnelles) : filtres simples,
 *          dérivée, squaring, intégration fenêtre glissante. Aucune allocation interne.
 *
 */

#include "ecg_utils.h"

#include <stddef.h>

/* ================================
 * Helpers (internes)
 * ================================ */

static double ecg_mean(const double *x, size_t n)
{
    if (!x || n == 0) return 0.0;

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum += x[i];
    }
    return sum / (double)n;
}

/* ================================
 * API
 * ================================ */

void ecg_apply_gain(double *x, size_t n, double gain)
{
    if (!x || n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        x[i] *= gain;
    }
}

void ecg_remove_dc(double *x, size_t n)
{
    if (!x || n == 0) return;

    const double m = ecg_mean(x, n);
    for (size_t i = 0; i < n; ++i) {
        x[i] -= m;
    }
}

void ecg_moving_average(const double *x, double *y, size_t n, size_t win)
{
    if (!x || !y || n == 0) return;
    if (win == 0) win = 1;

    /* Somme glissante: O(n) */
    double sum = 0.0;
    size_t w = 0;

    for (size_t i = 0; i < n; ++i) {
        sum += x[i];
        ++w;

        if (w > win) {
            sum -= x[i - win];
            --w;
        }

        y[i] = (w > 0) ? (sum / (double)w) : 0.0;
    }
}

void ecg_highpass_ma(const double *x, double *y, size_t n, size_t win)
{
    if (!x || !y || n == 0) return;
    if (win == 0) win = 1;

    /* y = x - MA(x) */
    double sum = 0.0;
    size_t w = 0;

    for (size_t i = 0; i < n; ++i) {
        sum += x[i];
        ++w;

        if (w > win) {
            sum -= x[i - win];
            --w;
        }

        const double ma = (w > 0) ? (sum / (double)w) : 0.0;
        y[i] = x[i] - ma;
    }
}

void ecg_derivative_1(const double *x, double *y, size_t n)
{
    if (!x || !y || n == 0) return;

    y[0] = 0.0;
    for (size_t i = 1; i < n; ++i) {
        y[i] = x[i] - x[i - 1];
    }
}

void ecg_square(const double *x, double *y, size_t n)
{
    if (!x || !y || n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        y[i] = x[i] * x[i];
    }
}

void ecg_mwi(const double *x, double *y, size_t n, size_t win)
{
    if (!x || !y || n == 0) return;
    if (win == 0) win = 1;

    /* Moyenne glissante (intégration fenêtre) */
    double sum = 0.0;
    size_t w = 0;

    for (size_t i = 0; i < n; ++i) {
        sum += x[i];
        ++w;

        if (w > win) {
            sum -= x[i - win];
            --w;
        }

        y[i] = (w > 0) ? (sum / (double)w) : 0.0;
    }
}
