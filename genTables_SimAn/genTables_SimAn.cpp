// genTables_SimAn.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "rna.h"
#include "siman.h"
#include "genTables_SimAn.h"
#include <stdio.h>

void generateTablesSimAn()
{
    int32_t i, j, k, l, m;
    /* We want tables of:
     * frequency codons
     * frequency codon followed by codon
     * frequency nucleotides
     * frequency nucleotide followed by nucleotide
     */
    printf("static double sCFrequency[NUM_CODON_TYPE] =\n{\n");
    for (i = 0; i < NUM_CODON_TYPE; i++)
    {
        if (i == NUM_CODON_TYPE - 1)    printf(" %9.7f", sVaccin.freqCf[i]);
        else if ((i & 7) == 0)          printf("    %9.7f,", sVaccin.freqCf[i]);
        else                            printf(" %9.7f,", sVaccin.freqCf[i]);
        if ((i & 7) == 7)               printf("\n");
    }
    printf("};\n\n");

    printf("static double sCCFrequency[NUM_CODON_TYPE][NUM_CODON_TYPE] =\n{\n");
    for (i = 0; i < NUM_CODON_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_CODON_TYPE; j++)
        {
            if (j == NUM_CODON_TYPE - 1)    printf(" %9.7f", sVaccin.freqCCf[i][j]);
            else if ((j & 7) == 0)          printf("        %9.7f,", sVaccin.freqCCf[i][j]);
            else                            printf(" %9.7f,", sVaccin.freqCCf[i][j]);
            if ((j & 7) == 7)               printf("\n");
        }
        if (i == NUM_CODON_TYPE - 1)    printf(" }\n");
        else                            printf(" },\n");
    }
    printf("};\n\n");

    printf("static double sCAFrequency[NUM_CODON_TYPE][NUM_TRANSLATIONS] =\n{\n");
    for (i = 0; i < NUM_CODON_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_TRANSLATIONS; j++)
        {
            if (j == NUM_TRANSLATIONS - 1)  printf(" %9.7f", sVaccin.freqCAf[i][j]);
            else if ((j & 7) == 0)          printf("        %9.7f,", sVaccin.freqCAf[i][j]);
            else                            printf(" %9.7f,", sVaccin.freqCAf[i][j]);
            if ((j & 7) == 7)               printf("\n");
        }
        if (i == NUM_CODON_TYPE - 1)    printf(" }\n");
        else                            printf(" },\n");
    }
    printf("};\n\n");

    printf("static double sNucleoFrequency[NUM_NUCLEO_TYPE] = {");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        if (i == NUM_NUCLEO_TYPE - 1)   printf(" %9.7f", sVaccin.freqNf[i]);
        else                            printf(" %9.7f,", sVaccin.freqNf[i]);
    }
    printf(" };\n\n");

    printf("static double sNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE] =\n {\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    { ");
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            if (j == NUM_NUCLEO_TYPE - 1)   printf(" %9.7f", sVaccin.freqNNf[i][j]);
            else                            printf(" %9.7f,", sVaccin.freqNNf[i][j]);
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf(" }\n");
        else                            printf(" },\n");
    }
    printf("};\n\n");

    printf("static double sNNCFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_CODON_TYPE] =\n {\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            printf("        {\n");
            for (k = 0; k < NUM_CODON_TYPE; k++)
            {
                if (!(k & 7))                   printf("            ");
                if (k == NUM_CODON_TYPE - 1)    printf(" %9.7f", sVaccin.freqNNCf[i][j][k]);
                else if ((k & 7) == 0)          printf("    %9.7f,", sVaccin.freqNNCf[i][j][k]);
                else                            printf(" %9.7f,", sVaccin.freqNNCf[i][j][k]);
                if ((k & 7) == 7)               printf("\n");
            }
            if (j == NUM_NUCLEO_TYPE - 1)   printf("        }\n");
            else                            printf("        },\n");
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf("    }\n");
        else                            printf("    },\n");
    }
    printf("};\n\n");

    printf("static double sNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE] =\n {\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            printf("        { ");
            for (k = 0; k < NUM_NUCLEO_TYPE; k++)
            {
                if (k == NUM_NUCLEO_TYPE - 1)   printf(" %9.7f", sVaccin.freqNNNf[i][j][k]);
                else                            printf(" %9.7f,", sVaccin.freqNNNf[i][j][k]);
            }
            if (j == NUM_NUCLEO_TYPE - 1)   printf(" }\n");
            else                            printf(" },\n");
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf("    }\n");
        else                            printf("    },\n");
    }
    printf("};\n\n");

    printf("static double sNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE] =\n {\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            printf("        { ");
            for (k = 0; k < NUM_NUCLEO_TYPE; k++)
            {
                printf("            { ");
                for (l = 0; l < NUM_NUCLEO_TYPE; l++)
                {
                    if (l == NUM_NUCLEO_TYPE - 1)   printf(" %9.7f", sVaccin.freqNNNNf[i][j][k][l]);
                    else                            printf(" %9.7f,", sVaccin.freqNNNNf[i][j][k][l]);
                }
                if (k == NUM_NUCLEO_TYPE - 1)   printf(" }\n");
                else                            printf(" },\n");
            }
            if (j == NUM_NUCLEO_TYPE - 1)   printf("        }\n");
            else                            printf("        },\n");
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf("    }\n");
        else                            printf("    },\n");
    }
    printf("};\n\n");

    printf("static double sNNNNNFrequency[NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE][NUM_NUCLEO_TYPE] =\n{\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            printf("        {\n");
            for (k = 0; k < NUM_NUCLEO_TYPE; k++)
            {
                printf("            {\n");
                for (l = 0; l < NUM_NUCLEO_TYPE; l++)
                {
                    printf("                { ");
                    for (m = 0; m < NUM_NUCLEO_TYPE; m++)
                    {
                        if (m == NUM_NUCLEO_TYPE - 1)   printf(" %9.7f", sVaccin.freqNNNNNf[i][j][k][l][m]);
                        else                            printf(" %9.7f,", sVaccin.freqNNNNNf[i][j][k][l][m]);
                    }
                    if (l == NUM_NUCLEO_TYPE - 1)   printf(" }\n");
                    else                            printf(" },\n");
                }
                if (k == NUM_NUCLEO_TYPE - 1)   printf("            }\n");
                else                            printf("            },\n");
            }
            if (j == NUM_NUCLEO_TYPE - 1)   printf("        }\n");
            else                            printf("        },\n");
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf("    }\n");
        else                            printf("    },\n");
    }
    printf("};\n\n");

    printf("static double sNucleoDistanceFrequency[NUM_NUCLEO_TYPE][DIST_REPORT_DEPTH] =\n {\n");
    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        printf("    { ");
        for (j = 0; j < DIST_REPORT_DEPTH; j++)
        {
            if (j == DIST_REPORT_DEPTH - 1) printf(" %9.7f", sVaccin.ratioRepDistNf[i][j]);
            else if ((j & 7) == 0)          printf("\n        %9.7f,", sVaccin.ratioRepDistNf[i][j]);
            else                            printf(" %9.7f,", sVaccin.ratioRepDistNf[i][j]);
        }
        if (i == NUM_NUCLEO_TYPE - 1)   printf("\n    }\n");
        else                            printf("\n    },\n");
    }
    printf("};\n\n");
}


void initializeSimAn()
{
    int32_t i, j, k, l, m;

    for (i = 0; i < NUM_CODON_TYPE; i++)
    {
        sVaccin.freqCf[i] = (double)sVaccin.freqC[i] / (double)sVaccin.countC;
        sCFrequency[i] = sVaccin.freqCf[i];
        for (j = 0; j < NUM_TRANSLATIONS; j++)
        {
            sVaccin.freqCAf[i][j] = (double)sVaccin.freqCA[i][j] / (double)(sVaccin.countC - 1);
            sCAFrequency[i][j] = sVaccin.freqCAf[i][j];
        }
        for (j = 0; j < NUM_CODON_TYPE; j++)
        {
            sVaccin.freqCCf[i][j] = (double)sVaccin.freqCC[i][j] / (double)(sVaccin.countC - 1);
            sCCFrequency[i][j] = sVaccin.freqCCf[i][j];
        }
    }

    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        sVaccin.freqNf[i] = (double)sVaccin.freqN[i] / (double)sVaccin.countN;
        sNucleoFrequency[i] = sVaccin.freqNf[i];
        for (j = 0; j < NUM_NUCLEO_TYPE; j++)
        {
            sVaccin.freqNNf[i][j] = (double)sVaccin.freqNN[i][j] / (double)(sVaccin.countN - 1);
            sNNFrequency[i][j] = sVaccin.freqNNf[i][j];
            for (k = 0; k < NUM_NUCLEO_TYPE; k++)
            {
                sVaccin.freqNNNf[i][j][k] = (double)sVaccin.freqNNN[i][j][k] / (double)(sVaccin.countN - 2);
                sNNNFrequency[i][j][k] = sVaccin.freqNNNf[i][j][k];
                for (l = 0; l < NUM_NUCLEO_TYPE; l++)
                {
                    sVaccin.freqNNNNf[i][j][k][l] = (double)sVaccin.freqNNNN[i][j][k][l] / (double)(sVaccin.countN - 3);
                    sNNNNFrequency[i][j][k][l] = sVaccin.freqNNNNf[i][j][k][l];
                    for (m = 0; m < NUM_NUCLEO_TYPE; m++)
                    {
                        sVaccin.freqNNNNNf[i][j][k][l][m] = (double)sVaccin.freqNNNNN[i][j][k][l][m] / (double)(sVaccin.countN - 3);
                        sNNNNNFrequency[i][j][k][l][m] = sVaccin.freqNNNNNf[i][j][k][l][m];
                    }
                }
            }
            for (k = 0; k < NUM_CODON_TYPE; k++)
            {
                sVaccin.freqNNCf[i][j][k] = (double)sVaccin.freqNNC[i][j][k] / (double)(sVaccin.countN - 3);
                sNNCFrequency[i][j][k] = sVaccin.freqNNCf[i][j][k];
            }
        }
    }

    for (i = 0; i < NUM_NUCLEO_TYPE; i++)
    {
        for (j = 0; j < DIST_REPORT_DEPTH; j++)
        {
            sVaccin.ratioRepDistNf[i][j] = (double)sVaccin.repDistN[i][j] / (double)(sVaccin.freqN[i] - 1);
            sNucleoDistanceFrequency[i][j] = sVaccin.ratioRepDistNf[i][j];
        }
    }
}

