#include "rna.h"
#include "siman.h"
#include "genTables_SimAn.h"


int main(int argc, char* argv[])
{
    loadCodons();
    initializeSimAn();
    generateTablesSimAn();
}
