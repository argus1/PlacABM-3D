#pragma once

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./placabm_setup.h"
#include "./placabm_substrates.h"

using namespace BioFVM;
using namespace PhysiCell;

void create_cell_types(void);
void setup_tissue(void);
void setup_microenvironment(void);

std::vector<std::string> my_coloring_function(Cell* pCell);

void phenotype_function(Cell* pCell, Phenotype& phenotype, double dt);
void custom_function(Cell* pCell, Phenotype& phenotype, double dt);
void contact_function(Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt);
