/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    csv_reader.h
 * @brief   Interface du module de lecture de fichiers CSV pour données ECG.
 * @details Ce module fournit les fonctions et variables nécessaires pour
 *          charger des données ECG depuis un fichier CSV et les stocker
 *          dans des structures globales exploitables par le reste du programme.
 *
 */


#ifndef CSV_READER_H
#define CSV_READER_H

#include "output_structs.h"

/* ================================
 * Variables globales
 * ================================ */
/**
 * @brief Données ECG lues depuis le fichier CSV.
 *
 * Tableau contenant les échantillons ECG pour chaque dérivation.
 * - Dimension 1 : nombre de dérivations (LEADS)
 * - Dimension 2 : nombre maximal d'échantillons (MAX_SAMPLES)
 */
extern double ecg_data[LEADS][MAX_SAMPLES];

/**
 * @brief Nombre d'échantillons effectivement lus.
 *
 * Cette variable indique combien d'échantillons valides ont été chargés
 * dans le tableau ecg_data.
 */
extern int sample_count;

/* ================================
 * API Publique
 * ================================ */
/**
 * @brief Lit un fichier CSV contenant des données ECG.
 *
 * Cette fonction ouvre un fichier CSV, parse son contenu et remplit
 * le tableau ecg_data ainsi que la variable sample_count.
 *
 * @param[in] filename Chemin vers le fichier CSV à lire.
 *
 * @return 0 en cas de succès,
 * @return valeur négative en cas d'erreur (fichier introuvable, format invalide, etc.).
 */
int read_csv(const char *filename);

#endif /* CSV_READER_H */
