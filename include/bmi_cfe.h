#ifndef CFE_BMI_CFE_H
#define CFE_BMI_CFE_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "./cfe.h"
#include "bmi.h"

//--------------------------------------------------
// Experiment to simplify BMI implementation (SDP)
//--------------------------------------------------
typedef struct Variable{
    unsigned int index;
    char name[80];
    char type[80];
    unsigned int size;
    // char units[80];
    // char role[80];
    // bool is_pointer;
} Variable;

/** Read number of lines in file and max line length, returning -1 if it does not exist or could not be read. */
int read_file_line_counts_cfe(const char* file_name, int* line_count, int* max_line_length);

/*  DEFINE A STRUCTURE TO HOLD THE ENTIRE MODEL STATE  */
//DATA STRUCTURE TO HOLD AORC FORCING DATA
struct aorc_forcing_data_cfe
{
// struct NAME                          DESCRIPTION                                            ORIGINAL AORC NAME
//____________________________________________________________________________________________________________________
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
} ;
static inline uint16_t sizeof_aorc_data(){
    uint16_t bytes = 0;
    bytes += 10*sizeof(double);
    bytes += 1*sizeof(long int);
    return bytes;
}

static void serialize_aorc_data(struct aorc_forcing_data_cfe* p, char* buffer){
    if(p!= NULL && buffer != NULL){
        char* pos = buffer;
        pos = copy_to_double(&(p->precip_kg_per_m2), pos);
        pos = copy_to_double(&(p->incoming_longwave_W_per_m2), pos);
        pos = copy_to_double(&(p->incoming_shortwave_W_per_m2), pos);
        pos = copy_to_double(&(p->surface_pressure_Pa), pos);
        pos = copy_to_double(&(p->specific_humidity_2m_kg_per_kg), pos);
        pos = copy_to_double(&(p->air_temperature_2m_K), pos);
        pos = copy_to_double(&(p->u_wind_speed_10m_m_per_s), pos);
        pos = copy_to_double(&(p->v_wind_speed_10m_m_per_s), pos);
        pos = copy_to_double(&(p->latitude), pos);
        pos = copy_to_double(&(p->longitude), pos);
        memcpy(pos, &(p->time), sizeof(long int));
        pos += sizeof(long int);
    }
}

static void deserialize_aorc_data(struct aorc_forcing_data_cfe* p, char* buffer){
    if(p!= NULL && buffer != NULL){
        char* pos = buffer;
        pos = copy_from_double(pos, &(p->precip_kg_per_m2));
        pos = copy_from_double(pos, &(p->incoming_longwave_W_per_m2));
        pos = copy_from_double(pos, &(p->incoming_shortwave_W_per_m2));
        pos = copy_from_double(pos, &(p->surface_pressure_Pa));
        pos = copy_from_double(pos, &(p->specific_humidity_2m_kg_per_kg));
        pos = copy_from_double(pos, &(p->air_temperature_2m_K));
        pos = copy_from_double(pos, &(p->u_wind_speed_10m_m_per_s));
        pos = copy_from_double(pos, &(p->v_wind_speed_10m_m_per_s));
        pos = copy_from_double(pos, &(p->latitude));
        pos = copy_from_double(pos, &(p->longitude));
        memcpy(&(p->time), pos, sizeof(long int));
        pos += sizeof(long int);
    }
}
typedef struct aorc_forcing_data_cfe aorc_forcing_data_cfe;

struct cfe_state_struct {

    // *************************************
    double timestep_rainfall_input_m;
    double soil_reservoir_storage_deficit_m;
    double infiltration_depth_m;
    double gw_reservoir_storage_deficit_m;

    double timestep_h;

    // ***********************************************************
    // ***************** Non-dynamic allocations *****************
    // ***********************************************************
    struct conceptual_reservoir soil_reservoir;
    struct conceptual_reservoir gw_reservoir;
    struct NWM_soil_parameters NWM_soil_params;
    struct evapotranspiration_structure et_struct;
    struct massbal vol_struct;

    /* xinanjiang_dev */
    struct direct_runoff_parameters_structure direct_runoff_params_struct;

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

};

static inline uint16_t sizeof_cfe_state(struct cfe_state_struct* p){
    uint16_t bytes = 0;
    bytes += 8*sizeof(double);
    bytes += 7*sizeof(int);
    bytes += 1*sizeof(long);
    bytes += strlen(p->forcing_file); //TODO might need to write this length to flatt buffer
    bytes += sizeof_reservoir(&(p->soil_reservoir));
    bytes += sizeof_reservoir(&(p->gw_reservoir));
    bytes += sizeof_soil_params();
    bytes += sizeof_et();
    bytes += sizeof_massbal();
    bytes += sizeof_runoff_params();
    bytes += sizeof_aorc_data();

    if(p->is_forcing_from_bmi){
        bytes += sizeof(double); //forcing_data pointer is dynamically allocated to single double
        bytes += sizeof(long); //focing_data_time pointer is dyanmically allocated to single long
    }
    else{
        //Error??? Or try to serialize the num_time_steps+1 data???
        assert(0); //not using bmi!!!
    }
    bytes += p->num_giuh_ordinates*sizeof(double);
    bytes += p->N_nash*sizeof(double);
    bytes += (p->num_giuh_ordinates+1)*sizeof(double);
    bytes += 7*sizeof(double); //dynamically allocated "flux" vars are single double, 7 flux vars
    //printf("Computed %ld bytes internally\n", bytes);
    return bytes;
}

static void serialize_cfe_state(struct cfe_state_struct* p, char* buffer){
    if(p!= NULL && buffer != NULL){
        char* pos = buffer;
        //printf("Putting %lf into buff\n", p->timestep_rainfall_input_m);
        pos = copy_to_double(&(p->timestep_rainfall_input_m), pos);
        //printf("Reading %lf from buff\n", *(double*) buffer);
        double tmp = -1;
        copy_to_double(buffer, &tmp);
        //printf("Reading %lf from var\n", tmp);
        //printf("Setting next %lf \n", p->soil_reservoir_storage_deficit_m);
        char* tmp2 = pos;
        pos = copy_to_double(&(p->soil_reservoir_storage_deficit_m), pos);
        //printf("Reading next %lf \n", *(double*) tmp2);
        pos = copy_to_double(&(p->infiltration_depth_m), pos);
        pos = copy_to_double(&(p->gw_reservoir_storage_deficit_m), pos);
        pos = copy_to_double(&(p->timestep_h), pos);
        serialize_reservoir(&(p->soil_reservoir), pos);
        pos += sizeof_reservoir(&(p->soil_reservoir));
        
        serialize_reservoir(&(p->gw_reservoir), pos);
        pos += sizeof_reservoir(&(p->gw_reservoir));
        
        serialize_soil_params(&(p->NWM_soil_params), pos);
        pos += sizeof_soil_params();
        
        serialize_et(&(p->et_struct), pos);
        pos += sizeof_et();

        serialize_massbal(&(p->vol_struct), pos);
        pos += sizeof_massbal();

        serialize_runoff_params(&(p->direct_runoff_params_struct), pos);
        pos += sizeof_runoff_params();

        pos = copy_to_long(&(p->epoch_start_time), pos);
        pos = copy_to_int(&(p->num_timesteps), pos);
        pos = copy_to_int(&(p->current_time_step), pos);
        pos = copy_to_int(&(p->time_step_size), pos);
        pos = copy_to_double(&(p->time_step_fraction), pos);
        pos = copy_to_int(&(p->is_forcing_from_bmi), pos);
        
        memcpy(pos, p->forcing_file, strlen(p->forcing_file)+1); //copy string plus null terminator
        pos += strlen(p->forcing_file)+1;

        pos = copy_to_double(&(p->K_lf), pos);
        pos = copy_to_double(&(p->K_nash), pos);
        pos = copy_to_int(&(p->N_nash), pos);
        pos = copy_to_int(&(p->num_giuh_ordinates), pos);

        serialize_aorc_data(&(p->aorc), pos);
        pos += sizeof_aorc_data();

        //FIXME copy dynamic pointers
        if(p->is_forcing_from_bmi){
            pos = copy_to_double(p->forcing_data_precip_kg_per_m2, pos);
            pos = copy_to_long(p->forcing_data_time, pos);
        }
        else{
            //Error??? Or try to serialize the num_time_steps+1 data???
            assert(0); //not using bmi!!!
        }

        pos = copy_to_double_array(p->giuh_ordinates, pos, p->num_giuh_ordinates);
        pos = copy_to_double_array(p->nash_storage, pos, p->N_nash);
        pos = copy_to_double_array(p->runoff_queue_m_per_timestep, pos, p->num_giuh_ordinates);

        //printf("BEFORE: %lf\n", p->flux_output_direct_runoff_m);
        pos = copy_to_double(p->flux_output_direct_runoff_m, pos);
        //printf("AFTER: %lf\n", (double*) pos);
        pos = copy_to_double(p->flux_giuh_runoff_m, pos);
        pos = copy_to_double(p->flux_nash_lateral_runoff_m, pos);
        pos = copy_to_double(p->flux_from_deep_gw_to_chan_m, pos);
        pos = copy_to_double(p->flux_perc_m, pos);
        pos = copy_to_double(p->flux_lat_m, pos);
        pos = copy_to_double(p->flux_Qout_m, pos);

        pos = copy_to_int(&(p->verbosity), pos);
    }
}

static void deserialize_cfe_state(struct cfe_state_struct* p, char* buffer){
    if(p!= NULL && buffer != NULL){
        char* pos = buffer;
        //printf("\n%p %p\n\n", pos, buffer);
        // printf("Putting %lf from buff\n", *(double*) buffer);
        // double tmp = -1;
        // copy_from_double(buffer, &tmp);
        // printf("Reading %lf from var\n", tmp);
        pos = copy_from_double(pos, &(p->timestep_rainfall_input_m));
        //printf("Reading %lf from struct\n", p->timestep_rainfall_input_m);
        //printf("\n%p %p\n\n", pos, buffer+sizeof(double));
        pos = copy_from_double(pos, &(p->soil_reservoir_storage_deficit_m));
        //printf("Reading next %lf \n", p->soil_reservoir_storage_deficit_m);
        //printf("From buff offset %lf, \n", *(double*)(buffer+sizeof(double)));
        pos = copy_from_double(pos, &(p->infiltration_depth_m));
        pos = copy_from_double(pos, &(p->gw_reservoir_storage_deficit_m));
        pos = copy_from_double(pos, &(p->timestep_h));
        deserialize_reservoir(&(p->soil_reservoir), pos);
        pos += sizeof_reservoir(&(p->soil_reservoir));
        
        deserialize_reservoir(&(p->gw_reservoir), pos);
        pos += sizeof_reservoir(&(p->gw_reservoir));
        
        deserialize_soil_params(&(p->NWM_soil_params), pos);
        pos += sizeof_soil_params();
        
        deserialize_et(&(p->et_struct), pos);
        pos += sizeof_et();

        deserialize_massbal(&(p->vol_struct), pos);
        pos += sizeof_massbal();

        deserialize_runoff_params(&(p->direct_runoff_params_struct), pos);
        pos += sizeof_runoff_params();

        pos = copy_from_long(pos, &(p->epoch_start_time));
        pos = copy_from_int(pos, &(p->num_timesteps));
        pos = copy_from_int(pos, &(p->current_time_step));
        pos = copy_from_int(pos, &(p->time_step_size));
        pos = copy_from_double(pos, &(p->time_step_fraction));
        pos = copy_from_int(pos, &(p->is_forcing_from_bmi));
        
        memcpy(p->forcing_file, pos, strlen(pos)+1); //copy string plus null terminator
        pos += strlen(p->forcing_file)+1;

        pos = copy_from_double(pos, &(p->K_lf));
        pos = copy_from_double(pos, &(p->K_nash));
        pos = copy_from_int(pos, &(p->N_nash));
        pos = copy_from_int(pos, &(p->num_giuh_ordinates));

        deserialize_aorc_data(&(p->aorc), pos);
        pos += sizeof_aorc_data();

        //FIXME copy dynamic pointers
        if(p->is_forcing_from_bmi){
            pos = copy_from_double(pos, p->forcing_data_precip_kg_per_m2);
            pos = copy_from_long(pos, p->forcing_data_time);
        }
        else{
            //Error??? Or try to serialize the num_time_steps+1 data???
            assert(0); //not using bmi!!!
        }

        pos = copy_from_double_array(pos, p->giuh_ordinates, p->num_giuh_ordinates);
        pos = copy_from_double_array(pos, p->nash_storage, p->N_nash);
        pos = copy_from_double_array(pos, p->runoff_queue_m_per_timestep, p->num_giuh_ordinates);

        //printf("BEFORE: %lf\n", (double*)pos);
        pos = copy_from_double(pos, p->flux_output_direct_runoff_m);
        //printf("AFTER: %lf\n", p->flux_output_direct_runoff_m);

        pos = copy_from_double(pos, p->flux_giuh_runoff_m);
        pos = copy_from_double(pos, p->flux_nash_lateral_runoff_m);
        pos = copy_from_double(pos, p->flux_from_deep_gw_to_chan_m);
        pos = copy_from_double(pos, p->flux_perc_m);
        pos = copy_from_double(pos, p->flux_lat_m);
        pos = copy_from_double(pos, p->flux_Qout_m);

        pos = copy_from_int(pos, &(p->verbosity));
    }
    else{
        assert(0);
    }
}

typedef struct cfe_state_struct cfe_state_struct;

extern double greg_2_jul(long year, long mon, long day, long h, long mi,
                         double se);
extern void calc_date_cfe(double jd, long *y, long *m, long *d, long *h, long *mi,
                      double *sec);

extern void itwo_alloc_cfe( int ***ptr, int x, int y);
extern void dtwo_alloc_cfe( double ***ptr, int x, int y);
extern void d_alloc_cfe(double **var,int size);
extern void i_alloc_cfe(int **var,int size);

extern void parse_aorc_line_cfe(char *theString,long *year,long *month, long *day,long *hour,
                            long *minute, double *dsec, struct aorc_forcing_data_cfe *aorc);

extern void get_word_cfe(char *theString,int *start,int *end,char *theWord,int *wordlen);

/*int read_init_config_cfe(const char* config_file, cfe_state_struct* model, double* alpha_fc, double* soil_storage,
                     int* is_soil_storage_ratio);*/
//LKC removed double alpha_fc since it has been added to the soil parameter structure                     
int read_init_config_cfe(const char* config_file, cfe_state_struct* model);


/*extern void init_soil_reservoir(cfe_state_struct* cfe_ptr, double alpha_fc, double max_storage, double storage,
                     int is_storage_ratios);*/
//LKC removed double alpha_fc since it has been added to the soil parameter structure                 
extern void init_soil_reservoir(cfe_state_struct* cfe_ptr);

//extern double init_reservoir_storage(int is_ratio, double amount, double max_amount);

extern void initialize_volume_trackers(cfe_state_struct* cfe_ptr);

extern void print_cfe_flux_header();
extern void print_cfe_flux_at_timestep(cfe_state_struct* cfe_ptr);
extern void mass_balance_check(cfe_state_struct* cfe_ptr);

// Bmi* register_bmi(Bmi *model);

Bmi* register_bmi_cfe(Bmi *model);

extern void run_cfe(cfe_state_struct* model);

cfe_state_struct * new_bmi_cfe(void);

#if defined(__cplusplus)
}
#endif

#endif
