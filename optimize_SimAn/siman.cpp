#include "rna.h"
#include "siman.h"

#include <stdio.h>
#include <math.h>

static rna_metrics_t sSimAn;
static rna_metrics_t sSimAnTmp;
static rna_metrics_t sSimAnUC;

uint16_t sSimAnMetrics = SM_N | SM_C | SM_CA | SM_NNNNN | SM_NDIST;

#ifdef GENERATE_TABLES
double sCFrequency[NUM_CODON_TYPE];
double sCCFrequency[NUM_CODON_TYPE][NUM_CODON_TYPE];
double sCAFrequency[NUM_CODON_TYPE][NUM_TRANSLATIONS];
double sNucleoFrequency[NUM_NUCLEO_TYPE];
double sNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
double sNNCFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_CODON_TYPE];
double sNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
double sNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
double sNNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE];
double sNucleoDistanceFrequency[NUM_NUCLEO_TYPE][DIST_REPORT_DEPTH];
#else
#include "tables_SimAn.h"
#endif



uint32_t nextRand = 1;

static int myrand(void)
{
    nextRand = nextRand * 1103515243 + 12345;
    return (unsigned int)(nextRand >> 16) & 0x7fff;
}

static void mysrand(unsigned int seed)
{
    nextRand = seed;
}

static int rnd(uint32_t max)
{
    return (myrand() * max) >> 15;
}


static double scoreSimAn(rna_metrics_t* r)
{
#if 1
    double rv = 0.0f, tmp;
    int32_t i, j, k, l, m;

    if (sSimAnMetrics & SM_C)
    {
        for (i = 0; i < NUM_CODON_TYPE; i++)
        {
            tmp = (double)r->freqC[i] / (double)r->countC;
            rv += fabs(sCFrequency[i] - tmp) * 100;
        }
    }

    if (sSimAnMetrics & SM_CA)
    {
        for (i = 0; i < NUM_CODON_TYPE; i++)
        {
            for (j = 0; j < NUM_TRANSLATIONS; j++)
            {
                tmp = (double)r->freqCA[i][j] / ((double)r->countC - 1);
                rv += fabs(sCAFrequency[i][j] - tmp) * 100;
            }
        }
    }

    if (sSimAnMetrics & SM_CC)
    {
        for (i = 0; i < NUM_CODON_TYPE; i++)
        {
            for (j = 0; j < NUM_CODON_TYPE; j++)
            {
                tmp = (double)r->freqCC[i][j] / ((double)r->countC - 1);
                rv += fabs(sCCFrequency[i][j] - tmp) * 100;
            }
        }
    }

    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        if (sSimAnMetrics & SM_N)
        {
            tmp = (double)r->freqN[i] / (double)r->countN;
            rv += fabs(sNucleoFrequency[i] - tmp) * 64;
        }
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            if (sSimAnMetrics & SM_NN)
            {
                tmp = (double)r->freqNN[i][j] / (double)(r->countN - 1);
                rv += fabs(sNNFrequency[i][j] - tmp) * 16;
            }
            for (k = 0; k < NUM_NUCLEO_TYPE; k++)
            {
                if (sSimAnMetrics & SM_NNN)
                {
                    tmp = (double)r->freqNNN[i][j][k] / (double)(r->countN - 2);
                    rv += fabs(sNNNFrequency[i][j][k] - tmp) * 4;
                }
                for (l = 0; l < NUM_NUCLEO_TYPE; l++)
                {
                    if (sSimAnMetrics & SM_NNNN)
                    {
                        tmp = (double)r->freqNNNN[i][j][k][l] / (double)(r->countN - 3);
                        rv += fabs(sNNNNFrequency[i][j][k][l] - tmp);
                    }
                    for (m = 0; m < NUM_NUCLEO_TYPE; m++)
                    {
                        if (sSimAnMetrics & SM_NNNNN)
                        {
                            tmp = (double)r->freqNNNNN[i][j][k][l][m] / (double)(r->countN - 3);
                            rv += fabs(sNNNNNFrequency[i][j][k][l][m] - tmp);
                        }
                    }
                }
            }
            if (sSimAnMetrics & SM_NNC)
            {
                for (k = 0; k < NUM_CODON_TYPE; k++)
                {
                    tmp = (double)r->freqNNC[i][j][k] / (double)(r->countC - 1);
                    rv += fabs(sNNCFrequency[i][j][k] - tmp) * 1;
                }
            }
        }
    }

    if (sSimAnMetrics & SM_NDIST)
    {
        for (i = 0; i < NUM_NUCLEO_TYPE; i++)
        {
            for (j = 0; j < DIST_REPORT_DEPTH; j++)
            {
                tmp = (double)r->repDistN[i][j] / (double)(r->freqN[i] - 1);
                rv += fabs(sNucleoDistanceFrequency[i][j] - tmp) * ((double)j + DIST_REPORT_DEPTH / 2) / ((double)20 * DIST_REPORT_DEPTH);
            }
        }
    }

    return rv;
#else
    return 100 - sequenceScore(r);
#endif
}


void checkSimAn()
{
    if (scoreSimAn(&sVaccin) > 0.000001)
    {
        fprintf(stderr, "Final check fails: sVaccin returns %9.7f on SimAn score function.\n", scoreSimAn(&sVaccin));
    }
}


void optimizeSimAn(bool verbose, unsigned int seed)
{
    int i, j;
    bool lastReportedTemp = false;
    double temp;
    double currentScore;
    double lastTemp = 1000;
    uint32_t amino2codon[NUM_TRANSLATIONS][10] = { 0 };

    mysrand(seed);

    for (i = 0; i < NUM_CODON_TYPE; i++)
    {
        uint8_t a = sCodonTranslation[i];
        amino2codon[a][0]++;
        amino2codon[a][amino2codon[a][0]] = i;
    }
    sSimAn = sOrig;
    currentScore = scoreSimAn(&sSimAn);
    for (temp = 100; temp > 0.01f; temp -= 0.25f)
    {
        uint32_t iter;
        if (!lastReportedTemp)
        {
            lastReportedTemp = true;
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "Temp %6.2f\r", temp);
        for (iter = 0; iter < 5000; iter++)
        {
            /* copy sSimAn to tmp RNA */
            sSimAnTmp = sSimAn;

            /* modify it, it is OK if we only change codons but not adapt nucleotides and other metrics */
            int numChanges = 1 + rnd(2) + rnd((int)temp / 10);
            for (j = 0; j < numChanges; j++)
            {
                /* select a position to change */
                uint32_t pos = rnd(sSimAnTmp.countC);
                uint8_t  a = sSimAnTmp.histA[pos];
                uint8_t  r = 1 + rnd(amino2codon[a][0]);
                uint8_t  c = amino2codon[a][r];
                sSimAnTmp.histC[pos] = c;
            }

            /* Construct the RNA metrics correctly by adding all codons and nucleotides separately */
            initRna(&sSimAnUC);
            for (i = 0; i < (int)sSimAnTmp.countC; i++)
            {
                uint8_t c = sSimAnTmp.histC[i];
                uint8_t a = sSimAnTmp.histA[i];
                char o1 = (c >> 4) & 3;
                char o2 = (c >> 2) & 3;
                char o3 = (c) & 3;
                processNucleotide(&sSimAnUC, sNucleotideC[o1]);
                processNucleotide(&sSimAnUC, sNucleotideC[o2]);
                processNucleotide(&sSimAnUC, sNucleotideC[o3]);
                processCodon(&sSimAnUC, c);
            }

            /* calculate the score */
            double score = scoreSimAn(&sSimAnUC);
            uint32_t weight = rnd(1000);
            if ((score < currentScore) ||
                ((score * 100000) < (weight * (temp + 1))))
            {
                currentScore = score;
                sSimAn = sSimAnUC;
                if (verbose)
                {
                    printf("New solution at temp %6.1f, SimAn-score = %10.6f, sequence conformity = %6.2f(N)/%6.2f(C)\r",
                        temp, currentScore, sequenceScoreNucleo(&sSimAn), sequenceScoreCodon(&sSimAn));
                    lastTemp = temp;
                    lastReportedTemp = false;
                }
            }
        }
    }
    sOptimized = sSimAn;
}
