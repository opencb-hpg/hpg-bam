
#include <stdlib.h>
#include <string.h>

#include "bam.h"
#include "bam_qc_batch.h"

/* ******************************************************
 *       	Function implementations      		*
 * ******************************************************/

bam_qc_batch_t* bam_qc_batch_new(int id) {
    bam_qc_batch_t* bam_qc_batch_p = (bam_qc_batch_t*) calloc(1, sizeof(bam_qc_batch_t));
    bam_qc_batch_p->id = id;

    return bam_qc_batch_p;
}

void bam_qc_batch_free(bam_qc_batch_t* batch_p, int all) {
    if (batch_p == NULL) {
        return;
    }

    if (batch_p->strand_counter_p != NULL) {
        free(batch_p->strand_counter_p);
    }

    if (batch_p->map_quality_p != NULL) {
        free(batch_p->map_quality_p);
    }

    if (batch_p->alignment_length_p != NULL) {
        free(batch_p->alignment_length_p);
    }

    if (batch_p->qc_alignment_p != NULL) {
        free(batch_p->qc_alignment_p);
    }

    if (all) {
        free_bam1(batch_p->alignments_p, batch_p->num_alignments);
    }

    free(batch_p);
}
