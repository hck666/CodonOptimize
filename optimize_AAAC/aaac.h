#ifndef _AAAC_H_
#define _AAAC_H_

#include "rna.h"

extern uint8_t sPreferedCodonWithPrev[NUM_TRANSLATIONS][NUM_TRANSLATIONS][NUM_TRANSLATIONS];
extern uint8_t sPreferedCodon[NUM_TRANSLATIONS];

void initializePreferedCodon();
void optimizeAAAC();

#endif /* defined _AAAC_H_ */
