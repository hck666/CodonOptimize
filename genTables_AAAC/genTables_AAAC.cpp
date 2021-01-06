// genTables_AAAC.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "rna.h"
#include "aaac.h"
#include "genTables_AAAC.h"
#include <stdio.h>

void initializePreferedCodon()
{
    uint32_t i, j, j2, k, l;
    for (i = 0; i < NUM_TRANSLATIONS; i++)
    {
        int32_t max = -1;
        for (j = 0; j < NUM_CODON_TYPE; j++)
        {
            if (sCodonTranslation[j] == i)
            {
                /* this one translates to it. Now count it in the vaccine */
                int32_t cnt = 0;
                for (k = 0; k < sVaccin.countC; k++)
                {
                    if (sVaccin.histC[k] == j)
                    {
                        cnt++;
                    }
                }
                if (cnt > max)
                {
                    sPreferedCodon[i] = j;
                    max = cnt;
                }
            }
        }
    }

    /* We now have the max in general. But now we are going to correlate it to the previous as well */
    uint32_t totalFail = 0;
    for (i = 0; i < NUM_TRANSLATIONS; i++) // cur
    {
        for (j = 0; j < NUM_TRANSLATIONS; j++) // prev
        {
            for (j2 = 0; j2 < NUM_TRANSLATIONS; j2++) // prev2
            {
                int32_t max = -1;
                for (k = 0; k < NUM_CODON_TYPE; k++) // codon
                {
                    if (sCodonTranslation[k] == i)
                    {
                        int32_t cnt = 0;
                        for (l = 2; l < sVaccin.countC; l++)
                        {
                            if ((sVaccin.histC[l] == k) && (sVaccin.histA[l - 1] == j) && (sVaccin.histA[l - 2] == j2))
                            {
                                cnt++;
                            }
                        }
                        if (cnt > max)
                        {
                            if (max > 0)
                            {
                                totalFail += max;
                            }
                            sPreferedCodonWithPrev[i][j][j2] = k;
                            max = cnt;
                        }
                        else
                        {
                            if (cnt > 0)
                            {
                                totalFail += cnt;
                            }
                        }
                    }
                }

                if (max < 0)
                {
                    sPreferedCodonWithPrev[i][j][j2] = sPreferedCodon[i];
                }
            }
        }
    }
    //printf("Expected total fails: %d\n", totalFail);
}


void generateTables()
{
    //
    //static uint8_t sPreferedCodon[NUM_TRANSLATIONS];
    int i, j, k;

    printf("static uint8_t sPreferedCodonWithPrev[NUM_TRANSLATIONS][NUM_TRANSLATIONS][NUM_TRANSLATIONS] =\n");
    printf("{\n");
    for (i = 0; i < NUM_TRANSLATIONS; i++)
    {
        printf("    {\n");
        for (j = 0; j < NUM_TRANSLATIONS; j++)
        {
            printf("        {");
            for (k = 0; k < NUM_TRANSLATIONS; k++)
            {
                if (k == 0)
                {
                    printf(" %2d", sPreferedCodonWithPrev[i][j][k]);
                }
                else
                {
                    printf(", %2d", sPreferedCodonWithPrev[i][j][k]);
                }
            }
            if (j == (NUM_TRANSLATIONS - 1))
            {
                printf("}\n");
            }
            else
            {
                printf("},\n");
            }
        }
        if (i == (NUM_TRANSLATIONS - 1))
        {
            printf("    }\n");
        }
        else
        {
            printf("    },\n");
        }
    }
    printf("};\n\n");

    printf("static uint8_t sPreferedCodon[NUM_TRANSLATIONS] =\n");
    printf("{\n    ");
    for (k = 0; k < NUM_TRANSLATIONS; k++)
    {
        if (k == 0)
        {
            printf(" %2d", sPreferedCodon[k]);
        }
        else
        {
            printf(", %2d", sPreferedCodon[k]);
        }
    }
    printf("\n};\n\n");
}