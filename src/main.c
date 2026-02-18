#include <stdio.h>
#include <string.h>

#include "csv_reader.h"
#include "ecg_processing.h"
#include "json_writer.h"
#include "output_structs.h"

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_csv> <output_json>\n", argv[0]);
        return 1;
    }

    if (read_csv(argv[1]) != 0) {
        fprintf(stderr, "Erreur lecture CSV.\n");
        return 2;
    }

    ECG_Peaks peaks;
    ECG_Intervals intervals;
    memset(&peaks, 0, sizeof(peaks));
    memset(&intervals, 0, sizeof(intervals));

    ECG_Params params;
    memset(&params, 0, sizeof(params));
    params.sampling_rate_hz = SAMPLING_RATE;
    params.leads            = LEADS;
    params.gain             = 100.0; // <-- Ajuster le gain si nécessaire
    params.r_threshold_hint = 0.0; // <-- Optionnel, peut être 0.0, et peut-être adaptatif au long du code.

    ECG_Context *ctx = ecg_create(&params);
    if (!ctx) {
        fprintf(stderr, "Erreur: ecg_create() a échoué.\n");
        return 4;
    }

    int lead_index = 1; // Analyser la LEAD II (index 1)
    if (lead_index < 0 || lead_index >= LEADS) {
        fprintf(stderr, "Erreur: lead_index invalide.\n");
        ecg_destroy(ctx);
        return 5;
    }

    /* Ici vous êtes libre de déconstruire en chunk ou d'analyser le signal dans son entiéreté
       Dans la réalité vous serez plus ammené a avoir un flux continus plutôt qu'un gros chunk de données */
    ECG_Status st = ecg_analyze(
        ctx,
        ecg_data[lead_index],
        (size_t)sample_count,
        lead_index,
        &peaks,
        &intervals
    );

    if (st != ECG_OK) {
        fprintf(stderr, "Erreur: ecg_analyze() a retourné %d.\n", (int)st);
        ecg_destroy(ctx);
        return 6;
    }

    printf("%d pics R détectés.\n", peaks.R_count);

    ecg_destroy(ctx);

    if (write_json(argv[2], &peaks, &intervals) != 0) {
        fprintf(stderr, "Erreur écriture JSON.\n");
        return 3;
    }

    printf("Analyse terminée. Résultats sauvegardés dans %s\n", argv[2]);
    return 0;
}
