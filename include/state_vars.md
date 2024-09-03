double timestep_rainfall_input_m;
double soil_reservoir_storage_deficit_m;
double infiltration_depth_m;
double gw_reservoir_storage_deficit_m;

double timestep_h;

// ***********************************************************
// ***************** Non-dynamic allocations *****************
// ***********************************************************
struct conceptual_reservoir soil_reservoir;
    int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
    double gw_storage;   // Initial Storage - LKC: added since I need to keep track of it when changing parameters
    double storage_max_m;   // maximum storage in this reservoir
    double storage_m;       // state variable.
    double storage_change_m; // storage change in the current step
    double coeff_primary;    // the primary outlet
    double exponent_primary;
    double storage_threshold_primary_m;
    double storage_threshold_secondary_m;
    double coeff_secondary;
    double exponent_secondary;
    double ice_fraction_schaake, ice_fraction_xinanjiang;
    int   is_sft_coupled; // boolean - true if SFT is ON otherwise OFF (default is OFF)
    
    //---Root zone adjusted AET development -rlm -ahmad -------------
    double *smc_profile; //soil moisture content profile
    int n_soil_layers; // number of soil layers
    double *soil_layer_depths_m; // soil layer depths defined in the config file in units of [m]
    int aet_root_zone; // boolean - true if aet_root_zone is ON otherwise OFF (default is OFF)
    int max_root_zone_layer;  // maximum root zone layer is used to identify the maximum layer to remove water from for AET
    double *delta_soil_layer_depth_m; // used to calculate the total soil moisture in each layer
    double soil_water_content_field_capacity;  // water content [m/m] at field capacity.  Used in AET routine
struct conceptual_reservoir gw_reservoir;
    int    is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
    double gw_storage;   // Initial Storage - LKC: added since I need to keep track of it when changing parameters
    double storage_max_m;   // maximum storage in this reservoir
    double storage_m;       // state variable.
    double storage_change_m; // storage change in the current step
    double coeff_primary;    // the primary outlet
    double exponent_primary;
    double storage_threshold_primary_m;
    double storage_threshold_secondary_m;
    double coeff_secondary;
    double exponent_secondary;
    double ice_fraction_schaake, ice_fraction_xinanjiang;
    int   is_sft_coupled; // boolean - true if SFT is ON otherwise OFF (default is OFF)
    
    //---Root zone adjusted AET development -rlm -ahmad -------------
    double *smc_profile; //soil moisture content profile
    int n_soil_layers; // number of soil layers
    double *soil_layer_depths_m; // soil layer depths defined in the config file in units of [m]
    int aet_root_zone; // boolean - true if aet_root_zone is ON otherwise OFF (default is OFF)
    int max_root_zone_layer;  // maximum root zone layer is used to identify the maximum layer to remove water from for AET
    double *delta_soil_layer_depth_m; // used to calculate the total soil moisture in each layer
    double soil_water_content_field_capacity;  // water content [m/m] at field capacity.  Used in AET routine
struct NWM_soil_parameters NWM_soil_params;
    // using same variable names as used in NWM.  <sorry>
    double smcmax;  // effective porosity [V/V]
    double wltsmc;  // wilting point soil moisture content [V/V]
    double satdk;   // saturated hydraulic conductivity [m s-1]
    double satpsi;	// saturated capillary head [m]
    double bb;      // beta exponent on Clapp-Hornberger (1978) soil water relations [-]
    double mult;    // the multiplier applied to satdk to route water rapidly downslope
    double slop;   // this factor (0-1) modifies the gradient of the hydraulic head at the soil bottom.  0=no-flow.
    double D;       // soil depth [m]
    double wilting_point_m;
    // LKC: Add this two parameters since they belong to soils. Makes the parameter specification consistent
    double alpha_fc;
    double refkdt;
    double soil_storage;
struct evapotranspiration_structure et_struct;
    double potential_et_m_per_s;
    double potential_et_m_per_timestep;
    double reduced_potential_et_m_per_timestep;
    double actual_et_from_rain_m_per_timestep;
    double actual_et_from_soil_m_per_timestep;
    double actual_et_m_per_timestep;
struct massbal vol_struct;
    double volstart            ;
    double vol_runoff          ;   
    double vol_infilt          ;   
    double vol_out_giuh        ;
    double vol_end_giuh        ;
    double vol_to_gw           ;
    double vol_in_gw_start     ;
    double vol_in_gw_end       ;
    double vol_from_gw         ;
    double vol_in_nash         ;
    double vol_in_nash_end     ;  // note the nash cascade is empty at start of simulation.
    double vol_out_nash        ;
    double vol_soil_start      ;
    double vol_to_soil         ;
    double vol_soil_to_lat_flow;
    double vol_soil_to_gw      ;  // this should equal vol_to_gw
    double vol_soil_end        ;
    double vol_et_from_soil    ;
    double vol_et_from_rain    ; 
    double vol_et_to_atm       ;   
    double volin               ;
    double volout              ;
    double volend              ;

/* xinanjiang_dev */
struct direct_runoff_parameters_structure direct_runoff_params_struct;
    surface_water_partition_type surface_partitioning_scheme;
        typedef enum {Schaake=1, Xinanjiang=2} surface_water_partition_type;
    double Schaake_adjusted_magic_constant_by_soil_type;
    double a_Xinanjiang_inflection_point_parameter;
    double b_Xinanjiang_shape_parameter;
    double x_Xinanjiang_shape_parameter;
    double urban_decimal_fraction;
    double ice_content_threshold; // ice content above which soil is impermeable
// Epoch-based start time (BMI start time is considered 0.0)
long epoch_start_time;
int num_timesteps;
int current_time_step;
int time_step_size;
double time_step_fraction;
// an integer to flag the forcings coming as "set_value" from BMI
int is_forcing_from_bmi;

char* forcing_file;

/* xinanjiang_dev
double Schaake_adjusted_magic_constant_by_soil_type;    */

//LKC Changed this to N_nash for consistency
//int num_lateral_flow_nash_reservoirs;

//LKC: added N_nash the same way as the other Nash parameters - making this consistent
double K_lf;
double K_nash;   
int N_nash;

int num_giuh_ordinates;

// ***********************************************************
// ******************* Dynamic allocations *******************
// ***********************************************************
struct aorc_forcing_data_cfe aorc;
    // struct NAME                          DESCRIPTION                                            ORIGINAL AORC NAME
    //______________________________________________________________________________________________________
    double precip_kg_per_m2;                // Surface precipitation "kg/m^2"                         | APCP_surface
    double incoming_longwave_W_per_m2 ;     // Downward Long-Wave Rad. Flux at 0m height, W/m^2       | DLWRF_surface
    double incoming_shortwave_W_per_m2;     // Downward Short-Wave Radiation Flux at 0m height, W/m^2 | DSWRF_surface
    double surface_pressure_Pa;             // Surface atmospheric pressure, Pa                       | PRES_surface
    double specific_humidity_2m_kg_per_kg;  // Specific Humidity at 2m height, kg/kg                  | SPFH_2maboveground
    double air_temperature_2m_K;            // Air temparture at 2m height, K                         | TMP_2maboveground
    double u_wind_speed_10m_m_per_s;        // U-component of Wind at 10m height, m/s                 | UGRD_10maboveground
    double v_wind_speed_10m_m_per_s;        // V-component of Wind at 10m height, m/s                 | VGRD_10maboveground
    double latitude;                        // degrees north of the equator.  Negative south          | latitude
    double longitude;                       // degrees east of prime meridian. Negative west          | longitude
    long int time; //TODO: type?           // seconds since 1970-01-01 00:00:00.0 0:00               | time

double* forcing_data_precip_kg_per_m2;
long* forcing_data_time;

//result_fluxes* fluxes;

double* giuh_ordinates;
double* nash_storage;
double* runoff_queue_m_per_timestep;

/* xinanjiang_dev
    changing the name to the more general "direct runoff"
double* flux_Schaake_output_runoff_m;*/
double* flux_output_direct_runoff_m ;

double* flux_giuh_runoff_m;
double* flux_nash_lateral_runoff_m;
double* flux_from_deep_gw_to_chan_m;
double* flux_perc_m;
double* flux_lat_m;
double* flux_Qout_m;

int verbosity;