#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"

char sNucleotideC[4] = { 'A', 'C', 'G', 'T' };
char sCodon2String[NUM_CODON_TYPE][NUM_NUCLEO_TYPE];

rna_metrics_t sOrig;
rna_metrics_t sVaccin;
rna_metrics_t sOptimized;


const char* sCodonTranslationNames[NUM_TRANSLATIONS] =
{
    "A(Ala)",
    "C(Cys)",
    "D(Asp)",
    "E(Glu)",
    "F(Phe)",
    "G(Gly)",
    "H(His)",
    "I(Ile)",
    "K(Lys)",
    "L(Leu)",
    "M(Met)",
    "N(Asn)",
    "P(Pro)",
    "Q(Gln)",
    "R(Arg)",
    "S(Ser)",
    "T(Thr)",
    "V(Val)",
    "W(Trp)",
    "Y(Tyr)",
    " -Stop"
};

codon_translation_t sCodonTranslation[NUM_CODON_TYPE] =
{
    AminoAcid_K,        // AAA
    AminoAcid_N,        // AAC
    AminoAcid_K,        // AAG
    AminoAcid_N,        // AAU
    AminoAcid_T,        // ACA
    AminoAcid_T,        // ACC
    AminoAcid_T,        // ACG
    AminoAcid_T,        // ACU
    AminoAcid_R,        // AGA
    AminoAcid_S,        // AGC
    AminoAcid_R,        // AGG
    AminoAcid_S,        // AGU
    AminoAcid_I,        // AUA
    AminoAcid_I,        // AUC
    AminoAcid_M,        // AUG
    AminoAcid_I,        // AUU

    AminoAcid_Q,        // CAA
    AminoAcid_H,        // CAC
    AminoAcid_Q,        // CAG
    AminoAcid_H,        // CAU
    AminoAcid_P,        // CCA
    AminoAcid_P,        // CCC
    AminoAcid_P,        // CCG
    AminoAcid_P,        // CCU
    AminoAcid_R,        // CGA
    AminoAcid_R,        // CGC
    AminoAcid_R,        // CGG
    AminoAcid_R,        // CGU
    AminoAcid_L,        // CUA
    AminoAcid_L,        // CUC
    AminoAcid_L,        // CUG
    AminoAcid_L,        // CUU

    AminoAcid_E,        // GAA
    AminoAcid_D,        // GAC
    AminoAcid_E,        // GAG
    AminoAcid_D,        // GAU
    AminoAcid_A,        // GCA
    AminoAcid_A,        // GCC
    AminoAcid_A,        // GCG
    AminoAcid_A,        // GCU
    AminoAcid_G,        // GGA
    AminoAcid_G,        // GGC
    AminoAcid_G,        // GGG
    AminoAcid_G,        // GGU
    AminoAcid_V,        // GUA
    AminoAcid_V,        // GUC
    AminoAcid_V,        // GUG
    AminoAcid_V,        // GUU

    StopCode   ,        // UAA
    AminoAcid_Y,        // UAC
    StopCode   ,        // UAG
    AminoAcid_Y,        // UAU
    AminoAcid_S,        // UCA
    AminoAcid_S,        // UCC
    AminoAcid_S,        // UCG
    AminoAcid_S,        // UCU
    StopCode   ,        // UGA
    AminoAcid_C,        // UGC
    AminoAcid_W,        // UGG
    AminoAcid_C,        // UGU
    AminoAcid_L,        // UUA
    AminoAcid_F,        // UUC
    AminoAcid_L,        // UUG
    AminoAcid_F,        // UUU
};


bool isNucleotide(char c)
{
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    {
        return true;
    }
    else
    {
        return true;
    }
}

uint8_t nucleotideValue(char c)
{
    uint8_t rv = 0xff;
    if (c == 'A')
    {
        rv = 0;
    }
    else if (c == 'C')
    {
        rv = 1;
    }
    else if (c == 'G')
    {
        rv = 2;
    }
    else if (c == 'T')
    {
        rv = 3;
    }
    else
    {
        fprintf(stderr, "Unknown nucleotide '%c'.\n", c);
    }
    return rv;
}

uint8_t codonValue(char c1, char c2, char c3)
{
    return (nucleotideValue(c1) << 4) + (nucleotideValue(c2) << 2) + nucleotideValue(c3);
}

char codonChar(uint8_t codon, uint32_t pos)
{
    return sNucleotideC[(codon >> (2 * pos)) & 3];
}

void initRna(rna_metrics_t* r)
{
    uint32_t i;

    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        r->lastPosN[i] = 0xFFFFFFFF;
        r->maxDistN[i] = 0;
    }

    memset(r->repDistN, 0, sizeof(r->repDistN));
    memset(r->freqN, 0, sizeof(r->freqN));
    memset(r->freqNN, 0, sizeof(r->freqNN));
    memset(r->freqNNC, 0, sizeof(r->freqNNC));
    memset(r->freqNNN, 0, sizeof(r->freqNNN));
    memset(r->freqNNNN, 0, sizeof(r->freqNNNN));
    memset(r->freqNNNNN, 0, sizeof(r->freqNNNNN));
    memset(r->freqC, 0, sizeof(r->freqC));
    memset(r->freqCC, 0, sizeof(r->freqCC));
    memset(r->freqCA, 0, sizeof(r->freqCA));
    memset(r->lastN, 0, sizeof(r->lastN));

    r->lastC = 0;
    r->countN = 0;
    r->countC = 0;
}

void processNucleotide(rna_metrics_t* r, char n)
{
    uint8_t in = nucleotideValue(n);

    if (r->lastPosN[in] != 0xFFFFFFFF)
    {
        uint32_t diff = r->countN - r->lastPosN[in];
        if (diff > r->maxDistN[in])
        {
            r->maxDistN[in] = diff;
        }
        if (diff >= DIST_REPORT_DEPTH)
        {
            diff = DIST_REPORT_DEPTH - 1;
        }
        r->repDistN[in][diff]++;
    }

    r->freqN[in]++;
    if (r->countN > 0)
    {
        r->freqNN[in][r->lastN[0]]++;
    }
    if (r->countN > 1)
    {
        r->freqNNN[in][r->lastN[0]][r->lastN[1]]++;
    }
    if (r->countN > 2)
    {
        r->freqNNNN[in][r->lastN[0]][r->lastN[1]][r->lastN[2]]++;
    }
    if (r->countN > 3)
    {
        r->freqNNNNN[in][r->lastN[0]][r->lastN[1]][r->lastN[2]][r->lastN[3]]++;
    }
    r->lastN[3] = r->lastN[2];
    r->lastN[2] = r->lastN[1];
    r->lastN[1] = r->lastN[0];
    r->lastN[0] = in;
    r->lastPosN[in] = r->countN;
    r->histN[r->countN] = in;
    r->countN++;
}

void processCodon(rna_metrics_t* r, uint8_t c)
{
    r->freqC[c]++;
    if (r->countC > 0)
    {
        uint8_t oldC = r->lastC;
        uint8_t n0 = (oldC >> 2) & 3;
        uint8_t n1 = (oldC) & 3;
        r->freqNNC[n1][n0][c]++;
        r->freqCC[c][oldC]++;
        r->freqCA[c][sCodonTranslation[oldC]]++;
    }
    r->lastC = c;
    r->histA[r->countC] = sCodonTranslation[c];
    r->histC[r->countC] = c;
    r->countC++;
}

double sequenceScoreNucleo(rna_metrics_t* r)
{
    int i;
    int count = 0;
    for (i = 0; i < (int)r->countN; i++)
    {
        if (r->histN[i] == sVaccin.histN[i])
        {
            count++;
        }
    }
    return 100.0f * (double)count / (double)r->countN;
}

double sequenceScoreCodon(rna_metrics_t* r)
{
    int i;
    int count = 0;
    for (i = 0; i < (int)r->countC; i++)
    {
        if (r->histC[i] == sVaccin.histC[i])
        {
            count++;
        }
    }
    return 100.0f * (double)count / (double)r->countC;
}

void dumpStatistics()
{
    int i, j;
    fprintf(stderr, "Freq: orig  vaccin opt'd:\n");
    for (i = 0; i < 64; i++)
    {
        fprintf(stderr, "%s : %5d %5d %5d\n", sCodon2String[i], sOrig.freqC[i], sVaccin.freqC[i], sOptimized.freqC[i]);
    }

    fprintf(stderr, "XX: --Vaccin-------      --Optimized----\n");
    fprintf(stderr, "XX:   A   C   G   T        A   C   G   T\n");
    for (i = 0; i < DIST_REPORT_DEPTH; i++)
    {
        if (sVaccin.repDistN[0][i] ||
            sVaccin.repDistN[1][i] ||
            sVaccin.repDistN[2][i] ||
            sVaccin.repDistN[3][i] ||
            sOptimized.repDistN[0][i] ||
            sOptimized.repDistN[1][i] ||
            sOptimized.repDistN[2][i] ||
            sOptimized.repDistN[3][i])
        {
            fprintf(stderr, "%02d: %3d %3d %3d %3d      %3d %3d %3d %3d\n", i,
                sVaccin.repDistN[0][i],
                sVaccin.repDistN[1][i],
                sVaccin.repDistN[2][i],
                sVaccin.repDistN[3][i],
                sOptimized.repDistN[0][i],
                sOptimized.repDistN[1][i],
                sOptimized.repDistN[2][i],
                sOptimized.repDistN[3][i]);
        }
    }

    uint32_t scoreOptN = 0;
    uint32_t scoreOptC = 0;
    uint32_t scoreOrigN = 0;
    uint32_t scoreOrigC = 0;
    uint32_t total[NUM_TRANSLATIONS] = { 0 };
    uint32_t offenders[NUM_TRANSLATIONS] = { 0 };
    for (i = 0; i < (int)sVaccin.countC; i++)
    {
        total[sCodonTranslation[sOptimized.histC[i]]]++;
        if (sVaccin.histC[i] == sOptimized.histC[i])
        {
            scoreOptC++;
        }
        else
        {
            offenders[sOptimized.histA[i]]++;
            fprintf(stderr, "[%04d]:", i);
            for (j = i - 2; j <= i + 2; j++)
            {
                if (j < 0)
                {
                    fprintf(stderr, " XXX ");
                }
                else if (j == i)
                {
                    fprintf(stderr, "|%s%s|", sCodon2String[sVaccin.histC[j]], &sCodonTranslationNames[sCodonTranslation[sVaccin.histC[j]]][1]);
                }
                else
                {
                    fprintf(stderr, " %s%s ", sCodon2String[sVaccin.histC[j]], &sCodonTranslationNames[sCodonTranslation[sVaccin.histC[j]]][1]);
                }
            }

            fprintf(stderr, "      ");
            for (j = i - 2; j <= i + 2; j++)
            {
                if (j < 0)
                {
                    fprintf(stderr, " XXX ");
                }
                else if (j == i)
                {
                    fprintf(stderr, "|%s%s|", sCodon2String[sOptimized.histC[j]], &sCodonTranslationNames[sCodonTranslation[sOptimized.histC[j]]][1]);
                }
                else
                {
                    fprintf(stderr, " %s%s ", sCodon2String[sOptimized.histC[j]], &sCodonTranslationNames[sCodonTranslation[sOptimized.histC[j]]][1]);
                }
            }

            fprintf(stderr, "\n");
        }
        if (sVaccin.histC[i] == sOrig.histC[i])
        {
            scoreOrigC++;
        }
    }

    for (i = 0; i < NUM_TRANSLATIONS; i++)
    {
        fprintf(stderr, "[%d- %d/%d] %s:", i, offenders[i], total[i], sCodonTranslationNames[i]);
        for (j = 0; j < 64; j++)
        {
            if (sCodonTranslation[j] == i)
            {
                uint32_t c;
                uint32_t tot = 0;
                fprintf(stderr, " %s:", sCodon2String[j]);
                for (c = 0; c < sOptimized.countC; c++)
                {
                    if (sOptimized.histC[c] == j)
                    {
                        tot++;
                    }
                }
                fprintf(stderr, "O%d", tot);
                tot = 0;
                for (c = 0; c < sVaccin.countC; c++)
                {
                    if (sVaccin.histC[c] == j)
                    {
                        tot++;
                    }
                }
                fprintf(stderr, "/V%d", tot);
            }
        }
        fprintf(stderr, "\n");
    }

    for (i = 0; i < (int)sVaccin.countN; i++)
    {
        if (sVaccin.histN[i] == sOptimized.histN[i])
        {
            scoreOptN++;
        }
        if (sVaccin.histN[i] == sOrig.histN[i])
        {
            scoreOrigN++;
        }
    }

    fprintf(stderr, "<<check>>\n");
    fprintf(stderr, "Similarity per nucleotide(Orig)      : %7.3f (%d / %d)\n", 100.0f * ((double)scoreOrigN) / (double)(sVaccin.countN), scoreOrigN, sVaccin.countN);
    fprintf(stderr, "Similarity per codon(Orig)           : %7.3f (%d / %d)\n", 100.0f * ((double)scoreOrigC) / (double)(sVaccin.countC), scoreOrigC, sVaccin.countC);
    fprintf(stderr, "<<results>>\n");
    fprintf(stderr, "Similarity per nucleotide(Optimized) : %7.3f (%d / %d) -> %d\n", 100.0f * ((double)scoreOptN) / (double)(sOptimized.countN), scoreOptN, sOptimized.countN, sOptimized.countN - scoreOptN);
    fprintf(stderr, "Similarity per codon(Optimized)      : %7.3f (%d / %d) -> %d\n", 100.0f * ((double)scoreOptC) / (double)(sOptimized.countC), scoreOptC, sOptimized.countC, sOptimized.countC - scoreOptC);
}


void dumpSolution()
{
    uint32_t i;
    printf("abspos,codonOrig,codonVaccine,codonOptimized\n");
    for (i = 0; i < sOrig.countC; i++)
    {
        printf("%d,%s,%s,%s\n", i*3, sCodon2String[sOrig.histC[i]], sCodon2String[sVaccin.histC[i]], sCodon2String[sOptimized.histC[i]]);
    }
}


void traceAminoAcid(codon_translation_t a)
{
    int32_t i, j;
    fprintf(stderr, "Tracing amino acid[%d] %s:", a, sCodonTranslationNames[a]);
    for (j = 0; j < 64; j++)
    {
        if (sCodonTranslation[j] == a)
        {
            uint32_t c;
            uint32_t tot = 0;
            for (c = 0; c < sVaccin.countC; c++)
            {
                if (sVaccin.histC[c] == j)
                {
                    tot++;
                }
            }
            fprintf(stderr, " %s:%d", sCodon2String[j], tot);
        }
    }
    fprintf(stderr, "\n");
    for (i = 0; i < (int32_t)sVaccin.countC; i++)
    {
        if (sCodonTranslation[sVaccin.histC[i]] == a)
        {
            //fprintf(stderr, "[%04d]:", i);
            for (j = i - 2; j <= i + 2; j++)
            {
                if (j < 0)
                {
                    fprintf(stderr, " XXX ");
                }
                else if (j == i)
                {
                    fprintf(stderr, "|%s%s|", sCodon2String[sVaccin.histC[j]], &sCodonTranslationNames[sCodonTranslation[sVaccin.histC[j]]][1]);
                }
                else
                {
                    fprintf(stderr, " %s%s ", sCodon2String[sVaccin.histC[j]], &sCodonTranslationNames[sCodonTranslation[sVaccin.histC[j]]][1]);
                }
            }
            fprintf(stderr, "\n");
        }
    }
}

uint64_t findCodonRelation(uint8_t codonOrig, uint8_t codonVaccin, uint8_t codonPrev)
{
    int i, j;
    uint64_t prevs = 0;
    const char* strOrig = (codonOrig == 0xff) ? "XXX" : sCodon2String[codonOrig];
    const char* strVaccin = (codonVaccin == 0xff) ? "XXX" : sCodon2String[codonVaccin];
    const char* strPrev = (codonPrev == 0xff) ? "XXX" : sCodon2String[codonPrev];

    fprintf(stderr, "Looking for %s->%s with prev %s\n", strOrig, strVaccin, strPrev);
    for (i = 0; i < (int)sOrig.countC; i++)
    {
        if (((codonOrig == 0xFF) || (sOrig.histC[i] == codonOrig)) &&
            ((codonVaccin == 0xFF) || (sVaccin.histC[i] == codonVaccin)) &&
            ((codonPrev == 0xFF) || ((i == 0) || (sOrig.histC[i - 1] == codonPrev)))
            )
        {
            if (i > 0)
            {
                prevs |= ((uint64_t)1) << sVaccin.histC[i - 1];
            }
            fprintf(stderr, "%4d:", i);
            for (j = i - 2; j <= i + 2; j++)
            {
                if (j < 0)
                {
                    fprintf(stderr, " XXX ");
                }
                else if (j == i)
                {
                    fprintf(stderr, "|%s|", sCodon2String[sOrig.histC[j]]);
                }
                else
                {
                    fprintf(stderr, " %s ", sCodon2String[sOrig.histC[j]]);
                }
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "%4d:", i);
            for (j = i - 2; j <= i + 2; j++)
            {
                if (j < 0)
                {
                    fprintf(stderr, " XXX ");
                }
                else if (j == i)
                {
                    fprintf(stderr, "|%s|", sCodon2String[sVaccin.histC[j]]);
                }
                else
                {
                    fprintf(stderr, " %s ", sCodon2String[sVaccin.histC[j]]);
                }
            }
            fprintf(stderr, "\n\n");
        }
    }
    return prevs;
}

void showAllInstances(char* t)
{
    uint8_t a = 0xFF;
    int i, i2, j, k;
    for (i = 0; i < NUM_TRANSLATIONS; i++)
    {
        if (strncmp(&sCodonTranslationNames[i][2], t, 3) == 0)
        {
            a = i;
        }
    }
    if (a != 0xFF)
    {
        /* Go through all codon types */
        for (i2 = 0; i2 < 64; i2++)
        {
            for (i = 0; i < 64; i++)
            {
                if (sCodonTranslation[i] != a)
                {
                    continue;
                }
                for (j = 1; j < (int)sVaccin.countC; j++)
                {
                    if ((sVaccin.histC[j - 1] == i2) && (sVaccin.histC[j] == i))
                    {
                        fprintf(stderr, "%4d:", j);
                        for (k = j - 2; k <= j + 2; k++)
                        {
                            if (k < 0)
                            {
                                fprintf(stderr, " XXX ");
                            }
                            else if (k == j)
                            {
                                fprintf(stderr, "|%s|", sCodon2String[sVaccin.histC[k]]);
                            }
                            else
                            {
                                fprintf(stderr, " %s ", sCodon2String[sVaccin.histC[k]]);
                            }
                        }
                        fprintf(stderr, "       ");
                        for (k = j - 2; k <= j + 2; k++)
                        {
                            if (k < 0)
                            {
                                fprintf(stderr, " XXX ");
                            }
                            else if (k == j)
                            {
                                fprintf(stderr, "|%s|", sCodon2String[sOptimized.histC[k]]);
                            }
                            else
                            {
                                fprintf(stderr, " %s ", sCodon2String[sOptimized.histC[k]]);
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                }
            }
        }
    }
}

void reportCodonSet(uint64_t list)
{
    uint8_t i;
    for (i = 0; i < 64; i++)
    {
        if (list & (((uint64_t)1) << i))
        {
            fprintf(stderr, "%s ", sCodon2String[i]);
        }
    }
    fprintf(stderr, "\n");
}

void reportCodonSetBin(uint64_t list)
{
    uint8_t i;
    for (i = 0; i < 64; i++)
    {
        if (list & (((uint64_t)1) << i))
        {
            fprintf(stderr, "1");
        }
        else
        {
            fprintf(stderr, "0");
        }
    }
    fprintf(stderr, "\n");
}

uint8_t parseCodon(char* str)
{
    uint8_t c = 0xFF;
    if (strlen(str) != 3)
    {
        fprintf(stderr, "Codon should have length 3\n");
    }
    else
    {
        uint8_t o1 = str[0];
        uint8_t o2 = str[1];
        uint8_t o3 = str[2];

        if (!isNucleotide(o1) || !isNucleotide(o2) || !isNucleotide(o3))
        {
            fprintf(stderr, "Invalid codon '%s' as argument\r\n", str);
        }
        else
        {
            c = codonValue(o1, o2, o3);
        }
    }
    return c;
}

void loadCodons()
{
    char line[100];
    uint8_t codonOrig;
    uint8_t codonVaccin;
    int i;
    char o1, o2, o3, v1, v2, v3;

    for (i = 0; i < 64; i++)
    {
        o1 = (i >> 4) & 3;
        o2 = (i >> 2) & 3;
        o3 = (i) & 3;
        sCodon2String[i][0] = sNucleotideC[o1];
        sCodon2String[i][1] = sNucleotideC[o2];
        sCodon2String[i][2] = sNucleotideC[o3];
        sCodon2String[i][3] = 0;
    }

    initRna(&sOrig);
    initRna(&sVaccin);
    initRna(&sOptimized);

    FILE* f;
    fopen_s(&f, "codons.txt", "r");
    if (!f)
    {
        fopen_s(&f, "..\\codons.txt", "r");
    }
    if (f)
    {
        fgets(line, 100, f); /* skip first line */
        fgets(line, 100, f); /* get first relevant line */
        do
        {
            /* find the first comma */
            i = 0;
            while (line[i] && (line[i] != ','))
            {
                i++;
            }
            o1 = line[i + 1];
            o2 = line[i + 2];
            o3 = line[i + 3];
            v1 = line[i + 5];
            v2 = line[i + 6];
            v3 = line[i + 7];
            codonOrig = codonValue(o1, o2, o3);
            codonVaccin = codonValue(v1, v2, v3);

            processNucleotide(&sVaccin, v1);
            processNucleotide(&sVaccin, v2);
            processNucleotide(&sVaccin, v3);
            processCodon(&sVaccin, codonVaccin);
            processNucleotide(&sOrig, o1);
            processNucleotide(&sOrig, o2);
            processNucleotide(&sOrig, o3);
            processCodon(&sOrig, codonOrig);
            fgets(line, 100, f); /* get next line */
        } while (!feof(f));
        fclose(f);
    }
    else
    {
        fprintf(stderr, "Could not open file.\n");
    }
}

