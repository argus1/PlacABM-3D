#include "./placabm_substrates.h"

#include <algorithm>
#include <stdexcept>

namespace
{
std::vector<std::string> required_substrates = {
    "oxygen",
    "glucose",
    "pathogen_toxin",
    "il6",
    "tnf_alpha",
    "cxcl8",
    "cxcl10",
    "il10",
    "ifn_gamma",
    "vegf"};

int require_substrate_index(const std::string& substrate_name)
{
    int index = BioFVM::microenvironment.find_density_index(substrate_name);
    if (index < 0)
    {
        throw std::runtime_error("Missing required substrate in XML: " + substrate_name);
    }
    return index;
}

std::vector<std::vector<double> > voxel_bulk_uptake_overrides;

void placabm_bulk_uptake_function(BioFVM::Microenvironment* pMicroenvironment,
                                  int voxel_index,
                                  std::vector<double>* write_destination)
{
    if (!pMicroenvironment || !write_destination)
    {
        return;
    }

    write_destination->assign(pMicroenvironment->number_of_densities(), 0.0);

    if (voxel_index < 0 || voxel_index >= (int)voxel_bulk_uptake_overrides.size())
    {
        return;
    }

    *write_destination = voxel_bulk_uptake_overrides[voxel_index];
}

} // namespace

namespace placabm
{

void register_placabm_substrates()
{
    std::cout << "[PlacABM] Verifying substrate registry..." << std::endl;
    for (const auto& substrate_name : required_substrates)
    {
        int index = require_substrate_index(substrate_name);
        std::cout << "  - " << substrate_name << " @ index " << index << std::endl;
    }
}

void set_uniform_diffusion_coefficients(double diffusion_coefficient, double decay_rate)
{
    for (size_t i = 0; i < BioFVM::microenvironment.diffusion_coefficients.size(); ++i)
    {
        BioFVM::microenvironment.diffusion_coefficients[i] = diffusion_coefficient;
        BioFVM::microenvironment.decay_rates[i] = decay_rate;
    }

    std::cout << "[PlacABM] Applied uniform diffusion coefficients: D=" << diffusion_coefficient
              << " decay=" << decay_rate << std::endl;
}

void apply_maternal_boundary_conditions(double oxygen_value, double glucose_value)
{
    int oxygen_index = require_substrate_index("oxygen");
    int glucose_index = require_substrate_index("glucose");

    for (size_t voxel = 0; voxel < BioFVM::microenvironment.mesh.voxels.size(); ++voxel)
    {
        BioFVM::microenvironment.update_dirichlet_node((int)voxel, oxygen_index, oxygen_value);
        BioFVM::microenvironment.update_dirichlet_node((int)voxel, glucose_index, glucose_value);
    }

    std::cout << "[PlacABM] Applied uniform maternal boundary conditions for oxygen and glucose." << std::endl;
}

void apply_maternal_boundary_conditions_to_voxels(const std::vector<int>& inlet_voxels,
                                                  double oxygen_value,
                                                  double glucose_value)
{
    if (inlet_voxels.empty())
    {
        apply_maternal_boundary_conditions(oxygen_value, glucose_value);
        return;
    }

    int oxygen_index = require_substrate_index("oxygen");
    int glucose_index = require_substrate_index("glucose");

    int applied = 0;
    for (int voxel : inlet_voxels)
    {
        if (voxel < 0 || voxel >= (int)BioFVM::microenvironment.number_of_voxels())
        {
            continue;
        }

        BioFVM::microenvironment.update_dirichlet_node(voxel, oxygen_index, oxygen_value);
        BioFVM::microenvironment.update_dirichlet_node(voxel, glucose_index, glucose_value);
        ++applied;
    }

    std::cout << "[PlacABM] Applied maternal inlet boundary conditions on " << applied << " voxels." << std::endl;
}

void set_tissue_diffusion_coefficients(const std::vector<int>& syncytiotrophoblast_mask,
                                       const std::vector<int>& fetal_capillary_mask)
{
    // BioFVM's default API uses global diffusion coefficients per substrate.
    // Here we approximate tissue-level heterogeneity by assigning voxel-level
    // bulk uptake rates that act as localized sinks/sources in selected masks.

    const int oxygen_index = require_substrate_index("oxygen");
    const int glucose_index = require_substrate_index("glucose");
    const int il6_index = require_substrate_index("il6");
    const int cxcl8_index = require_substrate_index("cxcl8");

    voxel_bulk_uptake_overrides.assign(BioFVM::microenvironment.number_of_voxels(),
                                       std::vector<double>(BioFVM::microenvironment.number_of_densities(), 0.0));

    // STB barrier: stronger local sink for key inflammatory signals and oxygen.
    for (int voxel : syncytiotrophoblast_mask)
    {
        if (voxel < 0 || voxel >= (int)voxel_bulk_uptake_overrides.size())
        {
            continue;
        }

        voxel_bulk_uptake_overrides[voxel][oxygen_index] = std::max(voxel_bulk_uptake_overrides[voxel][oxygen_index], 0.010);
        voxel_bulk_uptake_overrides[voxel][il6_index] = std::max(voxel_bulk_uptake_overrides[voxel][il6_index], 0.020);
        voxel_bulk_uptake_overrides[voxel][cxcl8_index] = std::max(voxel_bulk_uptake_overrides[voxel][cxcl8_index], 0.020);
    }

    // Fetal capillary: high uptake sink for oxygen and glucose exchange.
    for (int voxel : fetal_capillary_mask)
    {
        if (voxel < 0 || voxel >= (int)voxel_bulk_uptake_overrides.size())
        {
            continue;
        }

        voxel_bulk_uptake_overrides[voxel][oxygen_index] = std::max(voxel_bulk_uptake_overrides[voxel][oxygen_index], 0.050);
        voxel_bulk_uptake_overrides[voxel][glucose_index] = std::max(voxel_bulk_uptake_overrides[voxel][glucose_index], 0.050);
    }

    BioFVM::microenvironment.bulk_uptake_rate_function = placabm_bulk_uptake_function;

    std::cout << "[PlacABM] Configured tissue transport overrides via voxel bulk uptake"
              << " (stb=" << syncytiotrophoblast_mask.size()
              << ", capillary=" << fetal_capillary_mask.size() << ")." << std::endl;
}

} // namespace placabm
