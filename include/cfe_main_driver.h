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
 
 
#ifndef CFE_MAIN_DRIVER_H
#define CFE_MAIN_DRIVER_H

#include <stdio.h>
#include "cfe_types.h"


// Function prototypes

static int parse_command_line(int argc, char* argv[], cfe_cmdline_args_struct* args);

static int check_command_line_overrides(cfe_options_struct* options, const cfe_cmdline_args_struct* args);

static int open_output_files(const cfe_options_struct* options, const cfe_cmdline_args_struct* args,
                            FILE** q_fptr, FILE** Q_fptr, FILE** fluxes_fptr, FILE** storages_fptr, 
                            FILE** volbal_fptr, FILE** theta_fptr);
                            
static void close_output_files(FILE* q_fptr, FILE* Q_fptr, FILE* fluxes_fptr, FILE* storages_fptr, 
                              FILE* volbal_fptr, FILE* theta_fptr, FILE* warnings_fptr, 
                              FILE* forcing_fptr);


double calculate_pet_priestley_taylor(const cfe_forcing_struct* forcing, int dt_seconds, double alpha_pt);

static void print_usage(const char* prog);

static int open_forcing_file(const char* path, FILE** fptr, aorc_cols_t* cols);

static int read_next_forcing_aorc(FILE* f, const aorc_cols_t* cols, int dt_seconds,
                                  cfe_forcing_struct* forcing, 
                                  char* time_str);

static int calculate_time_delta(const char* time1_str, const char* time2_str);

static int parse_time_string(const char* time_str, aorc_forcing_time_struct* forcing_time);


static int analyze_forcing_file(const char* forcing_filename, aorc_forcing_time_struct* forcing_time, int verbosity);

static void format_timestamp(char* timestamp_str, int timestep, 
                             const aorc_forcing_time_struct* forcing_time,
                             const char* time_format);

int write_hotstart_config(const CFE_CONFIG* cfg,
                          const cfe_state_struct* state,
                          const aorc_forcing_time_struct* current_time);
                          
static void write_output_headers(const cfe_options_struct* options,
                                const char* time_format,
                                const char* delimiter,
                                FILE* q_fptr,
                                FILE* Q_fptr,
                                FILE* fluxes_fptr,
                                FILE* storages_fptr,
                                FILE* theta_fptr);

static void write_all_outputs(int timestep,
                              const cfe_options_struct* options,
                              const cfe_parameters_struct* params,
                              const cfe_state_struct* state,
                              const cfe_outputs_struct* outputs,
                              const cfe_forcing_struct* forcing,
                              const aorc_forcing_time_struct* forcing_time,
                              const cfe_volbal_struct* volbal,
                              FILE* q_fptr,
                              FILE* Q_fptr,
                              FILE* fluxes_fptr,
                              FILE* storages_fptr,
                              FILE* volbal_fptr,
                              FILE* theta_fptr);

// This function sets initial values of the volbal structt
void cfe_initialize_volume_balance(const cfe_parameters_struct* params,
                                   const cfe_options_struct* options, 
                                   const cfe_state_struct* state,
                                   cfe_volbal_struct* volbal);

static void write_volume_balance_summary(FILE* output_fptr,
                                        const cfe_options_struct* options,
                                        const cfe_volbal_struct* volbal,
                                        const cfe_state_struct* final_state,
                                        const cfe_parameters_struct* params,
                                        double *dz);

static int calculate_julian_day(int year, int month, int day);

static double calculate_julian_date(int year, int month, int day, int hour, int minute, int second);
#endif // CFE_MAIN_DRIVER_H
