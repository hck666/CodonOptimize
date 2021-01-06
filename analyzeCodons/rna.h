#ifndef _RNA_H_
#define _RNA_H_

#include <stdint.h>

#define DIST_REPORT_DEPTH       20
#define MAX_CODONS              5000
#define MAX_NUCLEO              (MAX_CODONS * 3)
#define NUM_NUCLEO_TYPE         4
#define NUM_CODON_TYPE          64

typedef enum
{
    AminoAcid_A,
    AminoAcid_C,
    AminoAcid_D,
    AminoAcid_E,
    AminoAcid_F,
    AminoAcid_G,
    AminoAcid_H,
    AminoAcid_I,
    AminoAcid_K,
    AminoAcid_L,
    AminoAcid_M, /* start! */
    AminoAcid_N,
    AminoAcid_P,
    AminoAcid_Q,
    AminoAcid_R,
    AminoAcid_S,
    AminoAcid_T,
    AminoAcid_V,
    AminoAcid_W,
    AminoAcid_Y,
    StopCode,
    NUM_TRANSLATIONS
} codon_translation_t;

typedef struct
{
    uint8_t                 lastN[4];
    uint8_t                 lastC;

    codon_translation_t     histA[MAX_CODONS];
    uint8_t                 histC[MAX_CODONS];
    uint8_t                 histN[MAX_NUCLEO];

    uint32_t                freqN[NUM_NUCLEO_TYPE];
    uint32_t                freqC[NUM_CODON_TYPE];
    uint32_t                freqCC[NUM_CODON_TYPE][NUM_CODON_TYPE];
    uint32_t                freqCA[NUM_CODON_TYPE][NUM_TRANSLATIONS];
    uint32_t                freqNN[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    uint32_t                freqNNC[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_CODON_TYPE];
    uint32_t                freqNNN[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    uint32_t                freqNNNN[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    uint32_t                freqNNNNN[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];

    uint32_t                maxDistN[NUM_NUCLEO_TYPE];
    uint32_t                lastPosN[NUM_NUCLEO_TYPE];
    uint32_t                repDistN[NUM_NUCLEO_TYPE][DIST_REPORT_DEPTH];

    double                  freqNf[NUM_NUCLEO_TYPE];
    double                  freqCf[NUM_CODON_TYPE];
    double                  freqCCf[NUM_CODON_TYPE][NUM_CODON_TYPE];
    double                  freqCAf[NUM_CODON_TYPE][NUM_TRANSLATIONS];
    double                  freqNNf[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    double                  freqNNCf[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_CODON_TYPE];
    double                  freqNNNf[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    double                  freqNNNNf[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
    double                  freqNNNNNf[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];

    double                  ratioRepDistNf[NUM_NUCLEO_TYPE][DIST_REPORT_DEPTH];

    uint32_t                countN;
    uint32_t                countC;
} rna_metrics_t;

extern char sNucleotideC[4];
extern char sCodon2String[NUM_CODON_TYPE][NUM_NUCLEO_TYPE];
extern codon_translation_t sCodonTranslation[NUM_CODON_TYPE];
extern const char* sCodonTranslationNames[NUM_TRANSLATIONS];

extern rna_metrics_t sOrig;
extern rna_metrics_t sVaccin;
extern rna_metrics_t sOptimized;

bool isNucleotide(char c);
uint8_t nucleotideValue(char c);
uint8_t codonValue(char c1, char c2, char c3);
char codonChar(uint8_t codon, uint32_t pos);
void initRna(rna_metrics_t* r);
void processNucleotide(rna_metrics_t* r, char n);
void processCodon(rna_metrics_t* r, uint8_t c);
double sequenceScoreNucleo(rna_metrics_t* r);
double sequenceScoreCodon(rna_metrics_t* r);
void dumpStatistics();
void dumpSolution();
void reportCodonSet(uint64_t list);
void reportCodonSetBin(uint64_t list);
void traceAminoAcid(codon_translation_t a);
uint64_t findCodonRelation(uint8_t codonOrig, uint8_t codonVaccin, uint8_t codonPrev);
void showAllInstances(char* t);
uint8_t parseCodon(char* str);
void loadCodons();

#endif /* defined _RNA_H_ */
