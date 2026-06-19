#include "./placabm_setup.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace
{

std::vector<int> voxel_layer_id;

std::vector<int> mask_decidua;
std::vector<int> mask_syncytiotrophoblast;
std::vector<int> mask_villous_stroma;
std::vector<int> mask_fetal_capillary;
std::vector<int> mask_maternal_inlet;

void clear_masks()
{
    mask_decidua.clear();
    mask_syncytiotrophoblast.clear();
    mask_villous_stroma.clear();
    mask_fetal_capillary.clear();
    mask_maternal_inlet.clear();

    voxel_layer_id.assign(BioFVM::microenvironment.number_of_voxels(), -1);
}

double clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

void assign_layer(int voxel_index, int layer)
{
    if (voxel_index < 0 || voxel_index >= (int)BioFVM::microenvironment.number_of_voxels())
    {
        return;
    }

    if (voxel_layer_id[voxel_index] != -1)
    {
        return;
    }

    voxel_layer_id[voxel_index] = layer;

    switch (layer)
    {
    case placabm::LAYER_DECIDUA:
        mask_decidua.push_back(voxel_index);
        break;
    case placabm::LAYER_SYNCYTIOTROPHOBLAST:
        mask_syncytiotrophoblast.push_back(voxel_index);
        break;
    case placabm::LAYER_VILLOUS_STROMA:
        mask_villous_stroma.push_back(voxel_index);
        break;
    case placabm::LAYER_FETAL_CAPILLARY:
        mask_fetal_capillary.push_back(voxel_index);
        break;
    default:
        break;
    }
}

void sort_and_unique(std::vector<int>& values)
{
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

bool load_single_mask_file(const std::string& filename, int layer)
{
    std::ifstream stream(filename.c_str());
    if (!stream)
    {
        return false;
    }

    std::string line;
    int loaded = 0;
    while (std::getline(stream, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        for (char& c : line)
        {
            if (c == ',' || c == ';' || c == '\t')
            {
                c = ' ';
            }
        }

        std::stringstream parser(line);
        int voxel = -1;
        if (!(parser >> voxel))
        {
            continue;
        }

        if (voxel >= 0 && voxel < (int)BioFVM::microenvironment.number_of_voxels())
        {
            assign_layer(voxel, layer);
            ++loaded;
        }
    }

    return loaded > 0;
}

void build_maternal_inlet_mask()
{
    mask_maternal_inlet.clear();

    if (mask_decidua.empty())
    {
        return;
    }

    double y_max_decidua = -9e99;
    for (int voxel : mask_decidua)
    {
        y_max_decidua = std::max(y_max_decidua, BioFVM::microenvironment.mesh.voxels[voxel].center[1]);
    }

    double y_spacing = 0.0;
    if (BioFVM::microenvironment.mesh.y_coordinates.size() > 1)
    {
        y_spacing = std::fabs(BioFVM::microenvironment.mesh.y_coordinates[1]
                              - BioFVM::microenvironment.mesh.y_coordinates[0]);
    }

    const double inlet_threshold = y_max_decidua - std::max(1.0, 0.75 * y_spacing);

    for (int voxel : mask_decidua)
    {
        const double y = BioFVM::microenvironment.mesh.voxels[voxel].center[1];
        if (y >= inlet_threshold)
        {
            mask_maternal_inlet.push_back(voxel);
        }
    }

    sort_and_unique(mask_maternal_inlet);
}

} // namespace

namespace placabm
{

void initialize_geometry_and_masks()
{
    clear_masks();

    const bool loaded = load_tissue_masks("./config/tissue_masks");
    if (!loaded)
    {
        std::cout << "[PlacABM] Tissue masks not found. Using procedural fallback geometry." << std::endl;
        generate_procedural_geometry(4, 25.0, 120.0);
    }

    build_maternal_inlet_mask();

    std::cout << "[PlacABM] Tissue masks ready: decidua=" << mask_decidua.size()
              << " stb=" << mask_syncytiotrophoblast.size()
              << " stroma=" << mask_villous_stroma.size()
              << " capillary=" << mask_fetal_capillary.size()
              << " inlet=" << mask_maternal_inlet.size() << std::endl;
}

bool load_tissue_masks(const std::string& mask_directory)
{
    const bool has_decidua = load_single_mask_file(mask_directory + "/decidua.csv", LAYER_DECIDUA);
    const bool has_stb = load_single_mask_file(mask_directory + "/syncytiotrophoblast.csv", LAYER_SYNCYTIOTROPHOBLAST);
    const bool has_stroma = load_single_mask_file(mask_directory + "/villous_stroma.csv", LAYER_VILLOUS_STROMA);
    const bool has_capillary = load_single_mask_file(mask_directory + "/fetal_capillary.csv", LAYER_FETAL_CAPILLARY);

    sort_and_unique(mask_decidua);
    sort_and_unique(mask_syncytiotrophoblast);
    sort_and_unique(mask_villous_stroma);
    sort_and_unique(mask_fetal_capillary);

    return has_decidua && has_stb && has_stroma && has_capillary;
}

void generate_procedural_geometry(int depth, double branch_angle_degrees, double branch_length)
{
    (void)branch_length;

    const auto& bb = BioFVM::microenvironment.mesh.bounding_box;
    const double x_min = bb[0];
    const double y_min = bb[1];
    const double x_max = bb[3];
    const double y_max = bb[4];

    const double x_range = std::max(1.0, x_max - x_min);
    const double y_range = std::max(1.0, y_max - y_min);

    const double angle_phase = branch_angle_degrees * 3.14159265358979323846 / 180.0;

    for (unsigned int voxel = 0; voxel < BioFVM::microenvironment.number_of_voxels(); ++voxel)
    {
        const auto& c = BioFVM::microenvironment.mesh.voxels[voxel].center;
        const double x_norm = clamp01((c[0] - x_min) / x_range);
        const double y_norm = clamp01((c[1] - y_min) / y_range);

        if (y_norm >= 0.78)
        {
            assign_layer((int)voxel, LAYER_DECIDUA);
            continue;
        }

        if (y_norm >= 0.64)
        {
            assign_layer((int)voxel, LAYER_SYNCYTIOTROPHOBLAST);
            continue;
        }

        // Lightweight procedural villous branching fallback:
        // create meandering branch centerlines in lower tissue and classify nearby voxels as capillary.
        const double frequency = 2.5 + 0.75 * std::max(1, depth);
        const double centerline = 0.5
            + 0.16 * std::sin(2.0 * 3.14159265358979323846 * frequency * y_norm + angle_phase)
            + 0.08 * std::sin(2.0 * 3.14159265358979323846 * (frequency * 0.5) * y_norm + 0.5 * angle_phase);

        const double radial_threshold = 0.035 + 0.02 * (1.0 - y_norm);

        if (std::fabs(x_norm - centerline) <= radial_threshold)
        {
            assign_layer((int)voxel, LAYER_FETAL_CAPILLARY);
        }
        else
        {
            assign_layer((int)voxel, LAYER_VILLOUS_STROMA);
        }
    }

    // Any unassigned voxels default to stroma.
    for (unsigned int voxel = 0; voxel < BioFVM::microenvironment.number_of_voxels(); ++voxel)
    {
        if (voxel_layer_id[voxel] < 0)
        {
            assign_layer((int)voxel, LAYER_VILLOUS_STROMA);
        }
    }

    sort_and_unique(mask_decidua);
    sort_and_unique(mask_syncytiotrophoblast);
    sort_and_unique(mask_villous_stroma);
    sort_and_unique(mask_fetal_capillary);
}

const std::vector<int>& decidua_mask_voxels()
{
    return mask_decidua;
}

const std::vector<int>& syncytiotrophoblast_mask_voxels()
{
    return mask_syncytiotrophoblast;
}

const std::vector<int>& villous_stroma_mask_voxels()
{
    return mask_villous_stroma;
}

const std::vector<int>& fetal_capillary_mask_voxels()
{
    return mask_fetal_capillary;
}

const std::vector<int>& maternal_inlet_voxels()
{
    return mask_maternal_inlet;
}

} // namespace placabm
