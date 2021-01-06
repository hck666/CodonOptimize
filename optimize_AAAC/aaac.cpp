#include "rna.h"
#include "aaac.h"

#ifdef GENERATE_TABLES
uint8_t sPreferedCodonWithPrev[NUM_TRANSLATIONS][NUM_TRANSLATIONS][NUM_TRANSLATIONS];
uint8_t sPreferedCodon[NUM_TRANSLATIONS];
#else
#include "tables_AAAC.h"
#endif

void optimizeAAAC()
{
    int i;
    char o1, o2, o3;
    uint8_t prevAminoAcid = 0;
    uint8_t prevAminoAcid2 = 0;

    for (i = 0; i < (int)sOrig.countC; i++)
    {
        uint8_t c;
        uint8_t a = sOrig.histA[i];

        // do the magic here.

        // do not translate with prev for the first two
        if (i > 1)
        {
            // start with doing the translation to amino acid
            c = sPreferedCodonWithPrev[a][prevAminoAcid][prevAminoAcid2];
        }
        else
        {
            c = sPreferedCodon[a];
        }

        o1 = (c >> 4) & 3;
        o2 = (c >> 2) & 3;
        o3 = (c) & 3;
        processNucleotide(&sOptimized, sNucleotideC[o1]);
        processNucleotide(&sOptimized, sNucleotideC[o2]);
        processNucleotide(&sOptimized, sNucleotideC[o3]);
        processCodon(&sOptimized, c);

        prevAminoAcid2 = prevAminoAcid;
        prevAminoAcid = a;
    }
}
