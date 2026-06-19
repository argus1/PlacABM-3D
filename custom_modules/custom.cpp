#include "./custom.h"

void create_cell_types(void)
{
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;
    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype = phenotype_function;
    cell_defaults.functions.custom_cell_rule = custom_function;
    cell_defaults.functions.contact_function = contact_function;

    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;

    initialize_cell_definitions_from_pugixml();
    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    setup_cell_rules();

    display_cell_definitions(std::cout);
}

void setup_microenvironment(void)
{
    initialize_microenvironment();

    placabm::register_placabm_substrates();
    placabm::set_uniform_diffusion_coefficients(100000.0, 0.1);

    placabm::initialize_geometry_and_masks();
    placabm::set_tissue_diffusion_coefficients(placabm::syncytiotrophoblast_mask_voxels(),
                                               placabm::fetal_capillary_mask_voxels());

    double maternal_oxygen = 38.0;
    double maternal_glucose = 5.0;
    if (default_microenvironment_options.Dirichlet_condition_vector.size() > placabm::SUBSTRATE_GLUCOSE)
    {
        maternal_oxygen = default_microenvironment_options.Dirichlet_condition_vector[placabm::SUBSTRATE_OXYGEN];
        maternal_glucose = default_microenvironment_options.Dirichlet_condition_vector[placabm::SUBSTRATE_GLUCOSE];
    }

    placabm::apply_maternal_boundary_conditions_to_voxels(placabm::maternal_inlet_voxels(),
                                                          maternal_oxygen,
                                                          maternal_glucose);
}

void setup_tissue(void)
{
    load_cells_from_pugixml();
    set_parameters_from_distributions();

    int number_of_cells = 0;
    if (parameters.ints.find_index("number_of_cells") != -1)
    {
        number_of_cells = parameters.ints("number_of_cells");
    }

    if (all_cells->empty() && number_of_cells > 0)
    {
        const double x_min = microenvironment.mesh.bounding_box[0];
        const double y_min = microenvironment.mesh.bounding_box[1];
        double z_min = microenvironment.mesh.bounding_box[2];

        const double x_max = microenvironment.mesh.bounding_box[3];
        const double y_max = microenvironment.mesh.bounding_box[4];
        double z_max = microenvironment.mesh.bounding_box[5];

        if (default_microenvironment_options.simulate_2D)
        {
            z_min = 0.0;
            z_max = 0.0;
        }

        for (int n = 0; n < number_of_cells; n++)
        {
            Cell* pC = create_cell(cell_defaults);
            std::vector<double> position = {
                x_min + UniformRandom() * (x_max - x_min),
                y_min + UniformRandom() * (y_max - y_min),
                z_min + UniformRandom() * (z_max - z_min)};
            pC->assign_position(position);
        }
    }
}

std::vector<std::string> my_coloring_function(Cell* pCell)
{
    return paint_by_number_cell_coloring(pCell);
}

void phenotype_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    return;
}

void custom_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    return;
}

void contact_function(Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt)
{
    return;
}
