#include "rna.h"
#include "aaac.h"
#include "genTables_AAAC.h"


int main(int argc, char* argv[])
{
    loadCodons();
    initializePreferedCodon();
    generateTables();
}
