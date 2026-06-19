#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"
#include "./custom_modules/custom.h"

using namespace BioFVM;
using namespace PhysiCell;

int main(int argc, char* argv[])
{
    bool xml_status = false;
    char copy_command[1024];

    if (argc > 1)
    {
        xml_status = load_PhysiCell_config_file(argv[1]);
        std::snprintf(copy_command, sizeof(copy_command), "cp %s %s", argv[1], PhysiCell_settings.folder.c_str());
    }
    else
    {
        xml_status = load_PhysiCell_config_file("./config/placabm_default.xml");
        std::snprintf(copy_command, sizeof(copy_command), "cp ./config/placabm_default.xml %s", PhysiCell_settings.folder.c_str());
    }

    if (!xml_status)
    {
        std::cerr << "Failed to load XML settings." << std::endl;
        return -1;
    }

    system(copy_command);

    omp_set_num_threads(PhysiCell_settings.omp_num_threads);

    setup_microenvironment();

    const double mechanics_voxel_size = 30.0;
    create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);

    create_cell_types();
    setup_tissue();

    set_save_biofvm_mesh_as_matlab(true);
    set_save_biofvm_data_as_matlab(true);
    set_save_biofvm_cell_data(true);
    set_save_biofvm_cell_data_as_custom_matlab(true);

    char filename[1024];
    if (!all_cells->empty())
    {
        std::snprintf(filename, sizeof(filename), "%s/initial", PhysiCell_settings.folder.c_str());
        save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
    }
    else
    {
        std::cout << "[PlacABM] No initial cells present; skipping initial MultiCellDS cell snapshot." << std::endl;
    }

    PhysiCell_SVG_options.length_bar = 200;
    std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
    std::string (*substrate_coloring_function)(double, double, double) = paint_by_density_percentage;

    std::snprintf(filename, sizeof(filename), "%s/initial.svg", PhysiCell_settings.folder.c_str());
    SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function);

    std::snprintf(filename, sizeof(filename), "%s/legend.svg", PhysiCell_settings.folder.c_str());
    create_plot_legend(filename, cell_coloring_function);

    display_citations();

    BioFVM::RUNTIME_TIC();
    BioFVM::TIC();

    std::ofstream report_file;
    if (PhysiCell_settings.enable_legacy_saves)
    {
        std::snprintf(filename, sizeof(filename), "%s/simulation_report.txt", PhysiCell_settings.folder.c_str());
        report_file.open(filename);
        report_file << "simulated time\tnum cells\tnum division\tnum death\twall time" << std::endl;
    }

    while (PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1 * diffusion_dt)
    {
        if (PhysiCell_globals.current_time > PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt)
        {
            display_simulation_status(std::cout);

            if (PhysiCell_settings.enable_legacy_saves)
            {
                log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
            }

            if (PhysiCell_settings.enable_full_saves)
            {
                if (!all_cells->empty())
                {
                    std::snprintf(filename, sizeof(filename), "%s/output%08u", PhysiCell_settings.folder.c_str(), PhysiCell_globals.full_output_index);
                    save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
                }
            }

            PhysiCell_globals.full_output_index++;
            PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
        }

        if (PhysiCell_globals.current_time > PhysiCell_globals.next_SVG_save_time - 0.5 * diffusion_dt)
        {
            if (PhysiCell_settings.enable_SVG_saves)
            {
                std::snprintf(filename, sizeof(filename), "%s/snapshot%08u.svg", PhysiCell_settings.folder.c_str(), PhysiCell_globals.SVG_output_index);
                SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function);
                PhysiCell_globals.SVG_output_index++;
                PhysiCell_globals.next_SVG_save_time += PhysiCell_settings.SVG_save_interval;
            }
        }

        microenvironment.simulate_bulk_sources_and_sinks(diffusion_dt);
        microenvironment.simulate_diffusion_decay(diffusion_dt);
        ((Cell_Container*)microenvironment.agent_container)->update_all_cells(PhysiCell_globals.current_time);
        PhysiCell_globals.current_time += diffusion_dt;
    }

    if (PhysiCell_settings.enable_legacy_saves)
    {
        log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
        report_file.close();
    }

    if (!all_cells->empty())
    {
        std::snprintf(filename, sizeof(filename), "%s/final", PhysiCell_settings.folder.c_str());
        save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, PhysiCell_globals.current_time);
    }

    std::snprintf(filename, sizeof(filename), "%s/final.svg", PhysiCell_settings.folder.c_str());
    SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function);

    std::cout << std::endl << "Total simulation runtime:" << std::endl;
    BioFVM::display_stopwatch_value(std::cout, BioFVM::runtime_stopwatch_value());

    return 0;
}
