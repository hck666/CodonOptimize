#ifndef _SIMAN_H_
#define _SIMAN_H_

#include "rna.h"




/* Simulated Annealing Metric */
#define SM_C        0x002
#define SM_CC       0x004
#define SM_CA       0x008
#define SM_N        0x010
#define SM_NN       0x020
#define SM_NNC      0x040
#define SM_NNN      0x080
#define SM_NNNN     0x100
#define SM_NNNNN    0x200
#define SM_NDIST    0x400



extern rna_metrics_t sSimAn;
extern rna_metrics_t sSimAnTmp;
extern rna_metrics_t sSimAnUC;

extern uint16_t sSimAnMetrics;

extern double sCFrequency[NUM_CODON_TYPE];
extern double sCCFrequency[NUM_CODON_TYPE][NUM_CODON_TYPE];
extern double sCAFrequency[NUM_CODON_TYPE][NUM_TRANSLATIONS];
extern double sNucleoFrequency[NUM_NUCLEO_TYPE];
extern double sNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
extern double sNNCFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_CODON_TYPE];
extern double sNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
extern double sNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
extern double sNNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
extern double sNucleoDistanceFrequency[NUM_NUCLEO_TYPE][DIST_REPORT_DEPTH];

void optimizeSimAn(bool verbose, unsigned int seed);
void checkSimAn();

#endif /* defined _SIMAN_H_ */