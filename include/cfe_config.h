/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1 
 *                and earlier
 *
 * Originally conceived and developed by: 
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS 
 *         Office of Water Prediction, Tuscaloosa, AL
 *
 */
 

#ifndef CFE_CONFIG_H
#define CFE_CONFIG_H

#include <stddef.h>

#include "nash_cascade.h"
#include "giuh.h"

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

// NDISC is the canonical constant (defined in cfe_soil_discrete.h).
// NDISCS is an alias used in the config struct (plural form).
#include "cfe_soil_discrete.h"
#ifndef NDISCS
#define NDISCS NDISC
#endif

/* String buffer sizes for config struct fields.
 * UNIT/CAT_NAME/OPTION are 64 bytes — sufficient for current usage.
 * Override via compiler flag (e.g. -DUNIT_STRING_LENGTH=256) if needed. */
#ifndef UNIT_STRING_LENGTH
#define UNIT_STRING_LENGTH 64
#endif

#ifndef CAT_NAME_STRING_LENGTH
#define CAT_NAME_STRING_LENGTH 64
#endif

#ifndef PATH_FILENAME_STRING_LENGTH
#define PATH_FILENAME_STRING_LENGTH 1024
#endif

#ifndef OPTION_STRING_LENGTH
#define OPTION_STRING_LENGTH 64
#endif
//===============================================================
// CFE v3 Configuration Structure — temporary exchange between parser and model init.
// Temporary exchange structure.   These elements are assigned to
// one of the following structures as defined in cfe_types.h   
//          const cfe_options_struct*    options
//          const cfe_parameters_struct* parms
//          const cfe_forcing_struct*    forcings
//          const cfe_state_struct*      states
//          const cfe_output_struct*     outputs
//===============================================================
typedef struct CFE_CONFIG
{
    double version;
    double timestep_h;
    char   timestep_units[UNIT_STRING_LENGTH];
    char   cat_id[CAT_NAME_STRING_LENGTH];
    double cat_latitude;
    char   cat_latitude_units[UNIT_STRING_LENGTH];
    double cat_longitude;
    char   cat_longitude_units[UNIT_STRING_LENGTH];
    double cat_elev;
    char   cat_elev_units[UNIT_STRING_LENGTH];
    double cat_area_km2;
    char   cat_area_units[UNIT_STRING_LENGTH];
    double cat_impervious_fraction;
    char   cat_impervious_units[UNIT_STRING_LENGTH];
    int    total_timesteps;
    int    verbosity;
    double et_alpha_pt;
    double soil_depth_m;
    char   soil_depth_units[UNIT_STRING_LENGTH];
    double soil_Clapp_Hornberger_exponent_b;
    double soil_sat_hydraulic_conductivity_cm_per_h;
    char   soil_sat_hydraulic_conductivity_units[UNIT_STRING_LENGTH];
    double soil_sat_capillary_head_cm;
    char   soil_sat_capillary_head_units[UNIT_STRING_LENGTH];
    double soil_to_gw_percolation_rate_limiter_0_to_1;
    double soil_effective_porosity;
    double soil_wilting_point_moisture_content;
    double soil_field_capacity_Pcap_over_Patm_0_1;
    double soil_reservoir_init_storage_m;
    char   soil_reservoir_init_storage_units[UNIT_STRING_LENGTH];
    double soil_reservoir_init_discrete_storage_theta[NDISCS];  // New array field
    char   soil_reservoir_init_discrete_storage_theta_units[UNIT_STRING_LENGTH];
    double soil_reservoir_rate_const_to_subsurface_lateral_flow;
    double soil_ice_content_impervious_threshold;
    double gw_reservoir_max_storage_m;
    char   gw_reservoir_max_storage_units[UNIT_STRING_LENGTH];
    double gw_reservoir_init_storage_m;
    char   gw_reservoir_init_storage_units[UNIT_STRING_LENGTH];
    double gw_discharge_coeff_m_per_timestep;
    char   gw_discharge_coeff_m_per_timestep_units[UNIT_STRING_LENGTH];
    double gw_discharge_exponent;
    double subsurface_routing_nash_K;
    double subsurface_routing_nash_cascade_init_storage_m[2];  // Changed from 3 to 2
    char   subsurface_routing_nash_cascade_init_storage_units[UNIT_STRING_LENGTH];
    char   control_input_forcing_filename[PATH_FILENAME_STRING_LENGTH];
    int    control_soil_simulate_freeze_thaw_true_false;
    int    control_soil_simulate_discrete_soil_moisture_true_false;  // New field ffor version 3
    int    control_ET_deepest_root_zone_discretization;
    int    control_soil_use_lookup_table_num_points;
    char   partitioning_scheme_name[OPTION_STRING_LENGTH];
    char   surface_routing_scheme_name[OPTION_STRING_LENGTH];
    int    surface_routing_num_giuh_ordinates;  
    double surface_routing_giuh_ordinates[MAX_NUM_GIUH_ORDINATES];   // much longer than needed, typ. < 10 ordinates
    char   surface_routing_giuh_ordinates_units[UNIT_STRING_LENGTH];
    double surface_routing_init_giuh_convolution_queue_m[MAX_NUM_GIUH_ORDINATES];
    char   surface_routing_init_giuh_convolution_queue_units[UNIT_STRING_LENGTH];
    double soil_Xinanjiang_tension_water_inflection_point;
    double soil_Xinanjiang_tension_water_soil_moist_distrib_exponent;
    double soil_Xinanjiang_free_water_soil_moist_distrib_exponent;


    // Output configuration fields
    char   output_path_name[PATH_FILENAME_STRING_LENGTH];
    char   output_status_warnings_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_internal_fluxes_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_internal_storages_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_volume_balance_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_soil_moisture_theta_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_discharge_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_time_standard_format[UNIT_STRING_LENGTH];
    char   output_file_delimiter[OPTION_STRING_LENGTH];
    char   output_value_format[OPTION_STRING_LENGTH];
    char   output_total_discharge_m3_per_sec_filename[PATH_FILENAME_STRING_LENGTH];
    char   output_new_config_filename[PATH_FILENAME_STRING_LENGTH];  // this enables checkpointing and hotstarting
} CFE_CONFIG;

#endif // CFE_CONFIG_STRUCTS_H
