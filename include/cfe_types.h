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
 

#ifndef CFE_TYPES_H
#define CFE_TYPES_H

#include <stddef.h>

#include "nash_cascade.h"
#include "giuh.h"
#include "cfe_soil_discrete.h"
#include "cfe_config.h"


#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

// used in parsing command line args
typedef struct {
    char* cfg_path;
    char* forcing_path;
    char* qout_path;
    char* volbal_path;
    char* fluxes_path;
    char* stores_path;
    char* thetas_path;
    int command_line_verbosity;
    int run_without_forcing;
} cfe_cmdline_args_struct;

/* struct to store the temporal properties of the forcing file */
typedef struct {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    int delta_t_seconds;
    int num_valid_lines;
    char time_format[64];  // Store the original time format from forcing file
} aorc_forcing_time_struct;

//###########
typedef struct {
    int    ndisc;                 // must equal NDISC
    int    deepest_root_disc;     // 1..ndisc
    int    use_ch_lookup_table;   // 1 => use LUT; 0 => analytic CH
    int    is_sft_coupled;        // TRUE iff coupled to soil freeze-thaw 
    double dt_hours;              // usually 1.0
} SoilControl;

//###########
typedef struct {
    double dz_m[NDISC];
    double zc_m[NDISC];
    double depth_m;
} SoilGeometry;

//###########
typedef struct {
    double theta_r;               // residual saturation (m3/m3)
    double theta_sat;             // saturation (m3/m3)
    double theta_fc;              // field capacity (m3/m3)
    double theta_wp;              // wilting point (m3/m3)
    double theta_aet_eq_pet;      // AET = PET at or above this (m3/m3)

    double K_sat_cm_per_h;        // cm/h
    double phi_sat_cm;            // cm
    double b_exp;                 // Clapp-Hornberger exponent

    double perc_limiter_0_to_1;   // 0..1 bottom drainage limiter
    double klf_per_h;           // lateral removal rate constant (m/h)
} SoilParameters;

// CH lookup tables over Theta=(theta-theta_r)/(theta_sat-theta_r)
//###########
typedef struct {
    int    n;
    double lnTheta_min;
    double dlnTheta;
    double inv_dlnTheta;
    double *lnpsi;                // ln(psi[m]) length n
    double *lnK;                  // ln(K[m/h]) length n
} SoilLookupTables;

//###########
typedef struct {
     double theta_in[NDISC];
     double psi_in[NDISC];   // m
     double K_in[NDISC];     // m/h
     int    ch_lut_hint_in[NDISC];
     double ice_fraction;
     double total_storage_m;
     double storage_deficit_m;
} SoilStateIn;

//###########
typedef struct {
    double theta_out[NDISC];
    int    ch_lut_hint_out[NDISC];
    double total_storage_m;
    double storage_deficit_m;
} SoilStateOut;

//###########
typedef struct {
    double rain_mm_per_h;         // mm/h
    double pet_mm_per_h;          // mm/h
} SoilForcing;

//###########
typedef struct {
    // Step-integrated exchanges (m)
    double AET_by_disc_m[NDISC];
    double lateral_by_disc_m[NDISC];
    double percolation_to_gw_m;
    double rain_into_soil_m;
    double rain_excess_m;

    // Internal vertical exchanges: [0..NDISC-2] interfaces i to i+1; [NDISC-1] bottom perc
    double interface_vol_m[NDISC];
    double interface_rate_m_per_h[NDISC];

    int    n_substeps_used;
} SoilFluxes;

typedef struct {
    double in_rain_m;         // infiltrated
    double excess_m;          // rejected
    double perc_m;
    double AET_m;
    double lateral_m;
    double delta_storage_m;   // Sum(theta_out-theta_in)*dz
    double residual_m;        // in - (outs) - Delta_S
} TimestepSoilVolbal;

/* Partitioning and routing options */
//###########
typedef enum {
    PARTITION_SCHAAKE = 0,
    PARTITION_XINANJIANG = 1
} cfe_partition_scheme_t;

//############
typedef enum {
    SURF_ROUTE_GIUH = 0,
    SURF_ROUTE_NASH_CASCADE = 1
} cfe_surface_route_t;

/* Options: switches and run controls (from config) */
//############
typedef struct {
    double cfe_version;           // 1.0, 1.1, 2.0, 2.1, 3.0.  If 0.0, then assume 2.0
    cfe_partition_scheme_t liquid_partitioning_scheme;
    cfe_surface_route_t surface_routing_scheme;
    int time_step_seconds;
    int num_timesteps;
    int verbosity;
    int enable_ET_Priestley_Taylor;
    int enable_freeze_thaw;
    int simulate_discrete_soil_moisture;
    int use_soil_lookup_table;
    char input_forcing_filename[PATH_FILENAME_STRING_LENGTH];

    // Output configuration
    char output_discharge_filename[PATH_FILENAME_STRING_LENGTH];
    char output_total_discharge_m3_per_sec_filename[PATH_FILENAME_STRING_LENGTH];
    char output_status_warnings_filename[PATH_FILENAME_STRING_LENGTH];
    char output_internal_fluxes_filename[PATH_FILENAME_STRING_LENGTH];
    char output_internal_storages_filename[PATH_FILENAME_STRING_LENGTH];
    char output_volume_balance_filename[PATH_FILENAME_STRING_LENGTH];
    char output_soil_moisture_theta_filename[PATH_FILENAME_STRING_LENGTH];
    char output_time_standard_format[64];
    char output_file_delimiter[64];
    char output_value_format[64];
    char output_path_name[PATH_FILENAME_STRING_LENGTH];
    char output_new_config_filename[PATH_FILENAME_STRING_LENGTH];
} cfe_options_struct;

/* Parameters: physical/empirical constants and arrays (from config) and initial states */
//############
typedef struct {
    /* Catchment */
    double catchment_area_km2;
    
    /* Soil */
    double soil_depth_m;
    double soil_b;
    double ksat_m_per_s;
    double sat_capillary_head_m;
    double effective_porosity;
    double wilting_point;
    double field_capacity_Pcap_over_Patm;
    double field_capacity_moisture_content; 
    double field_capacity_storage_m;        
    double soil_init_storage_m;
    double soil_k_lateral_per_h;
    double soil_ice_imperv_threshold;
    double refkdt;  // strictly not a parameter, it is a constant = 3.0
    double schaake_magic_constant;
    double soil_to_gw_percolation_rate_limiter_0_1;
    
    /* for discrete soil simulation */
    double soil_discrete_init_theta[NDISC];
    int    lut_n_points;
    double lut_theta_min;
    
    /* for using Priestley-Taylor method to calculate PET  */
    double alpha_pt;   // orginally 1.26, but often smaller for deserts and can be larger

    /* Xinanjiang, read even if not used */
    double xj_tension_inflection_0_1;
    double xj_tension_b;
    double xj_free_b;

    /* Groundwater */
    double gw_max_storage_m;
    double gw_init_storage_m;
    double gw_discharge_coeff_m_per_s;
    double gw_discharge_exponent;

    /* Surface routing */
    int giuh_num_ordinates;
    double giuh_ordinates[MAX_NUM_GIUH_ORDINATES];
    double giuh_init_queue_m[MAX_NUM_GIUH_ORDINATES];

    int nash_surface_N;
    double nash_surface_K_per_h;
    double nash_surface_init_storage_m[MAX_NUM_SURFACE_NASH_CASCADE];

    /* Subsurface Nash */
    int nash_subsurface_N;
    double nash_subsurface_K_per_h;
    double nash_subsurface_init_storage_m[2];

    /* Optional extras for surface Nash */
    double surface_Kinf_per_h;
    double surface_retention_depth_m;
} cfe_parameters_struct;

/* State: storages and internal queues that evolve over time */
//############
typedef struct {
    double soil_storage_m;
    double soil_storage_deficit_m;
    double gw_storage_m;
    double gw_storage_deficit_m;

    double nash_surface_storage_m[MAX_NUM_SURFACE_NASH_CASCADE];
    double nash_subsurface_storage_m[MAX_NUM_SUBSURFACE_NASH_CASCADE];
    double giuh_queue_m[MAX_NUM_GIUH_ORDINATES];
    double soil_discrete_storage_theta[NDISC];

    // Lookup tables ffor discrete soil moisture (calculate once, use many times)
    SoilLookupTables ch_lookup_tables;  

   // Needed ffor discrete soil simulation
    SoilControl      soil_control;
    SoilGeometry     soil_geometry; 
    SoilParameters   soil_parameters;
    SoilStateIn      soil_state_in;
    SoilStateOut     soil_state_out;
    SoilFluxes       soil_fluxes;
    
    int current_time_step;
} cfe_state_struct;    


/* Forcing: inputs per step */
//############
typedef struct {
    double rainfall_depth_m;      // PUT THIS FIRST - CFE expects it here!
    double et_potential_m;        // PUT THIS SECOND - CFE expects it here!
    double APCP_surface;          // mm (converted to m for rain_m)
    double DLWRF_surface;         // W/m^2
    double DSWRF_surface;         // W/m^2  
    double PRES_surface;          // Pa
    double SPFH_2maboveground;    // kg/kg
    double TMP_2maboveground;     // K
    double UGRD_10maboveground;   // m/s
    double VGRD_10maboveground;   // m/s
    double precip_rate;           // m/s
} cfe_forcing_struct;

typedef struct {  //<----- needed in forcing reader...
    int time_idx;
    int apcp_idx;
    int precip_rate_idx;
    int dlwrf_idx;
    int dswrf_idx;
    int pres_idx;
    int spfh_idx;
    int tmp_idx;
    int ugrd_idx;
    int vgrd_idx;
} aorc_cols_t;

//typedef struct { <- old pre 3.0
//    double rainfall_depth_m;
//    double et_potential_m;
//} cfe_forcing_struct;



// In cfe_types.h, add this BEFORE the CFE_Model_Context definition:

typedef struct {
    double volstart            ;
    double volstart_soil       ;
    double volstart_gw         ;
    double volstart_surface    ;
    double volstart_subsurface ;
    double vol_direct_runoff   ;  // this is water that won't fit into the soil suring a timestep because it is full
                                  // particularly in the discretized soil situation, where the upper disc fills
    double vol_runoff          ;
    double vol_infilt          ;
    double vol_runon_infilt    ;
    double vol_out_surface     ;
    double vol_end_surface     ;
    double vol_to_gw           ;
    double vol_in_gw_start     ;
    double vol_in_gw_end       ;
    double vol_from_gw         ;
    double vol_in_subsurf_nash ;   
    double vol_out_subsurf_nash;   
    double vol_soil_start      ;
    double vol_to_soil         ;
    double vol_soil_to_lat_flow;
    double vol_soil_to_gw      ;
    double vol_soil_end        ;
    double vol_et_from_soil    ;
    double vol_et_from_rain    ;
    double vol_et_to_atm       ;
    double volin               ;
    double volin_PET           ;
    double volout              ;
    double volend_soil         ;
    double volend_gw           ;
    double volend_surface      ;
    double volend_subsurface   ;
    double volend              ;
    double vol_et_from_retention_depth;
} cfe_volbal_struct;

/* output configuration helper structure */
typedef struct {
    char time_format[64];      // "timestep", "datetime", "juliandate"
    char delimiter[8];         // ",", " ", "\t", "|"
    int current_timestep;
    aorc_forcing_time_struct* forcing_time;
} cfe_output_config_struct;

/* Outputs: results per step */
//############
typedef struct {
    double surface_runoff_generated_m;
    double surface_routed_to_outlet_m;
    double lateral_flow_m;
    double baseflow_m;
    double qout_m;
    double actual_et_m;
} cfe_outputs_struct;

/* CFE model context: the whole model instance bundled together */
//############
typedef struct {
    CFE_CONFIG              config;        /* everything read in from a config file in one package */
    cfe_options_struct      options;       /* switches and run control */
    cfe_parameters_struct   parameters;    /* physical parameters */
    cfe_state_struct        state;         /* model storages and internal state */
    cfe_forcing_struct      forcing;       /* current step forcing */
    cfe_outputs_struct      last_outputs;  /* results from last step */
    cfe_volbal_struct       volbal;        /* the volume balance struct */
    double timestep_storage_start_m;       /* four per-timestep volume balance measures*/
    double timestep_input_m;
    double timestep_output_m; 
    double timestep_storage_end_m;
} CFE_Model_Context;

#endif
