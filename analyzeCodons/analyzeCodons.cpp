// analyzeCodons.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "rna.h"
#include "aaac.h"
#include "siman.h"
#include "genTables_AAAC.h"
#include "genTables_SiMan.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[])
{
    int i;
    uint8_t codonOrig;

    loadCodons();

    if ((argc == 2) && (strcmp(argv[1], "generate") == 0))
    {
        initializePreferedCodon();
        optimizeAAAC();
        generateTables();
        return 0;
    }

    if ((argc == 2) && (strcmp(argv[1], "generateSimAn") == 0))
    {
        initializeSimAn();
        generateTablesSimAn();
        return 0;
    }

    if ((argc == 2) && (strcmp(argv[1], "SimAn") == 0))
    {
        initializeSimAn();
        generateTablesSimAn();
        optimizeSimAn(true, 1);
        dumpStatistics();
        return 0;
    }

    initializePreferedCodon();
    optimizeAAAC();

    if (argc == 1)
    {
        dumpStatistics();
    }
    else if (argc == 2)
    {
        if (argv[1][0] == 'a')
        {
            codon_translation_t a = NUM_TRANSLATIONS;
            for (i = 0; i < NUM_TRANSLATIONS; i++)
            {
                if (strncmp(&sCodonTranslationNames[i][2], &argv[1][1], 3) == 0)
                {
                    a = (codon_translation_t)i;
                }
            }
            if (a != NUM_TRANSLATIONS)
            {
                traceAminoAcid(a);
            }
            else
            {
                printf("Unknown amino acid %s\n", &argv[1][1]);
            }
        }
        else if (argv[1][0] == 'c')
        {
            codonOrig = parseCodon(&argv[1][1]);
            if (codonOrig != 0xFF)
            {
                uint64_t prevs = findCodonRelation(codonOrig, 0xFF, 0xFF);
                reportCodonSet(prevs);
            }
        }
        else if (argv[1][0] == 's')
        {
            showAllInstances(&argv[1][1]);
        }
    }
    else
    {
        uint8_t match[4][2] = { { 0xff, 0xff },{ 0xff, 0xff },{ 0xff, 0xff },{ 0xff, 0xff } };
        uint8_t cnt[2] = { 0,0 };
        bool ok = true;
        for (i = 1; i < argc; i++)
        {
            if (strlen(argv[i]) == 5)
            {
                int idx = argv[i][1] - '0';
                int slot;
                uint8_t codon = parseCodon(&argv[i][2]);
                if (argv[i][0] == 'o')
                {
                    slot = 0;
                }
                else if (argv[i][0] == 'v')
                {
                    slot = 1;
                }
                else if (argv[i][0] == 'p')
                {
                    slot = 2;
                }
                else
                {
                    printf("Slot must be o(riginal), v(accin) or p(revious).\n");
                    ok = false;
                }

                if ((idx != 0) && (idx != 1))
                {
                    printf("Index must be 0 or 1.\n");
                    ok = false;
                }

                if (ok)
                {
                    match[slot][idx] = codon;
                }
            }
        }

        if (ok)
        {
            if ((match[0][1] != 0xFF) || (match[1][1] != 0xff) || (match[2][1] != 0xff))
            {
                uint64_t prevs[2];
                prevs[0] = findCodonRelation(match[0][0], match[1][0], match[2][0]);
                prevs[1] = findCodonRelation(match[0][1], match[1][1], match[2][1]);
                reportCodonSetBin(prevs[0]);
                reportCodonSetBin(prevs[1]);
                printf("----------------------------------------------------------------\n");
                reportCodonSetBin(prevs[0] & prevs[1]);
                printf("Conflicts: ");
                reportCodonSet(prevs[0] & prevs[1]);
            }
            else
            {
                uint64_t p = findCodonRelation(match[0][0], match[1][0], match[2][0]);
                reportCodonSetBin(p);
            }
        }
    }

    return 0;
}

