// optimize_AAAC.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "rna.h"
#include "aaac.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    int i;
    bool verbose = false;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-h") == 0)
        {
            printf("Usage: %s [-v] [-h].\n", argv[0]);
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            verbose = true;
        }
    }

    loadCodons();
    optimizeAAAC();
    if (verbose)
    {
        dumpStatistics();
    }
    else
    {
        dumpSolution();
    }
}