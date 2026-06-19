#pragma once

#include "../core/PhysiCell.h"

#include <string>
#include <vector>

namespace placabm
{

enum TissueLayerId
{
    LAYER_DECIDUA = 0,
    LAYER_SYNCYTIOTROPHOBLAST = 1,
    LAYER_VILLOUS_STROMA = 2,
    LAYER_FETAL_CAPILLARY = 3
};

void initialize_geometry_and_masks();
bool load_tissue_masks(const std::string& mask_directory);
void generate_procedural_geometry(int depth, double branch_angle_degrees, double branch_length);

const std::vector<int>& decidua_mask_voxels();
const std::vector<int>& syncytiotrophoblast_mask_voxels();
const std::vector<int>& villous_stroma_mask_voxels();
const std::vector<int>& fetal_capillary_mask_voxels();
const std::vector<int>& maternal_inlet_voxels();

} // namespace placabm
