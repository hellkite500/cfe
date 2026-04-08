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
 

#ifndef CFE_HELPERS_H
#define CFE_HELPERS_H

#include "cfe_types.h"
#include "cfe_config.h"
#include "cfe.h"
#include "cfe_soil_discrete.h"

/* Defaults */

int is_fabs_less_than_epsilon(double a,double epsilon);

void set_parameters_defaults_2_1(cfe_parameters_struct* params);  // sets all defaultt parameter values

void set_options_defaults_2_1(cfe_options_struct* opts); // sets defaults that are only meaningful cfe version>=2.1

void set_state_defaults(cfe_state_struct* state);        // sets all state variable defaults

int normalize_config_units(CFE_CONFIG* config, int is_legacy);

int map_config_to_cfe_structs(CFE_CONFIG* config, cfe_parameters_struct* params, cfe_options_struct* options);

/* Validate cfe required parameters (errors -> -1, warnings printed to stderr) */
int validate_required_parameters(const CFE_CONFIG* cfg, const int verbosity);   // all versions

/* Map config struct elements to parameters, options, and states (intelligently copy array elements -hopefully) */
int map_config_to_parameters_and_options(const CFE_CONFIG* cfg,              // all versions
                                             cfe_parameters_struct* params,
                                             cfe_options_struct* opts);

/* Parse, validate, config file for all cfe versions */
int parse_config_driver(const char* config_file, double config_file_version,  // calls different parser funcs
                             CFE_CONFIG* config,
                             cfe_parameters_struct* params,                   // ffor cfe versions <2.1 and
                             cfe_options_struct* opts);                       // versions > 2.0

/* Initialize state from params/options, allocating queues as needed */
int cfe_initialize(const cfe_parameters_struct* params,
                   const cfe_options_struct* opts,
                   cfe_state_struct* state);

/* Run one step. TODO: wire to cfe() kernel inside this function. */
int cfe_step(const cfe_parameters_struct* params,
             const cfe_options_struct* opts,
             cfe_state_struct* state,
             const cfe_forcing_struct* forcing,
             double dt_seconds,
             cfe_outputs_struct* outputs,
             cfe_volbal_struct* volbal);

                                   
double cfe_get_last_qout_m(const cfe_outputs_struct* outputs);

// ADDED FOR DEBUGGING PURPOSES:
void print_cfe_input_debug(const cfe_options_struct*    o,
                           const cfe_parameters_struct* p,
                           const cfe_state_struct*      s,
                           const char*                  cfg_path,
                           const char*                  forcing_path,
                           const char*                  qout_path,
                           const char*                  volbal_path);

// ADDED FOR DEBUGGING PURPOSES:
void print_exchange_values(int timestep,
                          const cfe_parameters_struct* p,
                          const cfe_options_struct* o,
                          const cfe_state_struct* s,
                          const cfe_forcing_struct* forcing,
                          double dt_seconds,
                          const cfe_outputs_struct* outputs);

// ADDED FOR DEBUGGING PURPOSES:
void check_dsbm_local_volume_balance(
    TimestepSoilVolbal *soil_volbal,
    SoilFluxes *soil_fluxes,
    SoilStateIn *soil_state_in,
    SoilStateOut *soil_state_out,
    SoilGeometry *soil_geometry,
    double infiltration_depth_m,
    double actual_et_from_soil_m,
    double balance_tolerance
);

/* Output writing helper functions */
const char* get_delimiter_string(const char* delimiter_name);

/* Cleanup any allocations made in state */
int cfe_finalize(cfe_state_struct* s);


#endif
