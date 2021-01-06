// optimize_SimAn.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "rna.h"
#include "siman.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[])
{
    int i;
    unsigned int seed = 1;
    bool verbose = false;

    for (i=1; i<argc; i++)
    {
        if (strcmp(argv[i], "-h") == 0)
        {
            printf("Usage: %s [-v] [-h] [-s <seed>] [-m <metriclist>].\n", argv[0]);
            printf("Valid metrics are N,2N,3N,4N,5N,NDIST,C,CC,CA,NNC\n");
            printf("Default is N,C,CA,4N,NDIST\n");
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            verbose = true;
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            i++;
            if (i < argc)
            {
                seed = atoi(argv[i]);
            }
            else
            {
                printf("%s: -s needs an argument (seed).\r\n", argv[0]);
            }
        }
        else if (strcmp(argv[i], "-m") == 0)
        {
            i++;
            if (i < argc)
            {
                uint16_t mask = 0;
                char s[1000] = ",";
                strcat_s(s, argv[i]);
                strcat_s(s, ",");
                if (strstr(s, ",N,"))       { mask |= SM_N; }
                if (strstr(s, ",2N,"))      { mask |= SM_NN; }
                if (strstr(s, ",3N,"))      { mask |= SM_NNN; }
                if (strstr(s, ",4N,"))      { mask |= SM_NNNN; }
                if (strstr(s, ",5N,"))      { mask |= SM_NNNNN; }
                if (strstr(s, ",NDIST,"))   { mask |= SM_NDIST; }
                if (strstr(s, ",C,"))       { mask |= SM_C; }
                if (strstr(s, ",CC,"))      { mask |= SM_CC; }
                if (strstr(s, ",CA,"))      { mask |= SM_CA; }
                if (strstr(s, ",NNC,"))     { mask |= SM_NNC; }
                sSimAnMetrics = mask;
            }
            else
            {
                printf("%s: -m needs an argument (metrics list).\r\n", argv[0]);
            }
        }
    }

    loadCodons();
    optimizeSimAn(verbose, seed);
    if (verbose)
    {
        dumpStatistics();
    }
    else
    {
        dumpSolution();
    }
    return 0;
}
