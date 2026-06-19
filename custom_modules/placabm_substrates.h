#pragma once

#include "../core/PhysiCell.h"

namespace placabm
{

enum SubstrateId
{
    SUBSTRATE_OXYGEN = 0,
    SUBSTRATE_GLUCOSE = 1,
    SUBSTRATE_PATHOGEN_TOXIN = 2,
    SUBSTRATE_IL6 = 3,
    SUBSTRATE_TNF_ALPHA = 4,
    SUBSTRATE_CXCL8 = 5,
    SUBSTRATE_CXCL10 = 6,
    SUBSTRATE_IL10 = 7,
    SUBSTRATE_IFN_GAMMA = 8,
    SUBSTRATE_VEGF = 9
};

void register_placabm_substrates();
void set_uniform_diffusion_coefficients(double diffusion_coefficient, double decay_rate);
void apply_maternal_boundary_conditions(double oxygen_value, double glucose_value);
void apply_maternal_boundary_conditions_to_voxels(const std::vector<int>& inlet_voxels,
                                                  double oxygen_value,
                                                  double glucose_value);
void set_tissue_diffusion_coefficients(const std::vector<int>& syncytiotrophoblast_mask,
                                       const std::vector<int>& fetal_capillary_mask);

} // namespace placabm
