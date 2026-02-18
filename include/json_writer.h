/**
 * @authors
 * - Bruno Carvalho
 *
 * @file    json_writer.h
 * @brief   Interface du module d’écriture JSON des résultats d’analyse ECG.
 * @details Ce module fournit la fonction permettant d’exporter les résultats
 *          (pics détectés, intervalles calculés, etc.) dans un fichier JSON
 *          à partir des structures définies dans output_structs.h.
 *
 */


#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "output_structs.h"

/* ================================
 * API Publique
 * ================================ */
/**
 * @brief Écrit les résultats d’analyse ECG dans un fichier JSON.
 *
 * Cette fonction sérialise les structures de résultats (pics ECG et intervalles)
 * et les écrit dans un fichier JSON au chemin fourni.
 *
 * @param[in] filename   Chemin vers le fichier JSON de sortie.
 * @param[in] peaks      Pointeur vers la structure contenant les pics détectés (non NULL).
 * @param[in] intervals  Pointeur vers la structure contenant les intervalles calculés (non NULL).
 *
 * @return 0 en cas de succès,
 * @return valeur négative en cas d’erreur
 *
 * @note Le format exact du JSON dépend de l’implémentation (.c) et des champs
 *       présents dans ECG_Peaks et ECG_Intervals.
 */
int write_json(const char *filename, const ECG_Peaks *peaks, const ECG_Intervals *intervals);

#endif /* JSON_WRITER_H */
