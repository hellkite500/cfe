#ifndef _CONCEPTUAL_RESERVOIR_H
#define _CONCEPTUAL_RESERVOIR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "serialize.h"

#define TRUE 1
#define FALSE 0

struct conceptual_reservoir {
  // this data structure describes a nonlinear reservoir having two outlets, one primary with an activation
  // threshold that may be zero, and a secondary outlet with a threshold that may be zero
  // this will also simulate a linear reservoir by setting the exponent parameter to 1.0 iff is_exponential==FALSE
  // iff is_exponential==TRUE, then it uses the exponential discharge function from the NWM V2.0 forumulation
  // as the primary discharge with a zero threshold, and does not calculate a secondary discharge.
  //--------------------------------------------------------------------------------------------------
  int is_exponential;  // set this true TRUE to use the exponential form of the discharge equation
  int is_sft_coupled; // boolean - true if SFT is ON otherwise OFF (default is OFF)
  int n_soil_layers; // number of soil layers
  int aet_root_zone; // boolean - true if aet_root_zone is ON otherwise OFF (default is OFF)
  int max_root_zone_layer;  // maximum root zone layer is used to identify the maximum layer to remove water from for AET
  
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
  double soil_water_content_field_capacity;  // water content [m/m] at field capacity.  Used in AET routine 
  
  //---Root zone adjusted AET development -rlm -ahmad -------------
  double *smc_profile; //soil moisture content profile, array of size n_soil_layers
  double *soil_layer_depths_m; // soil layer depths defined in the config file in units of [m], array of size n_soil_layers
  double *delta_soil_layer_depth_m; // used to calculate the total soil moisture in each layer, array of size n_soil_layers
  //---------------------------------------------------------------
};

static u_int32_t sizeof_reservoir(struct conceptual_reservoir* reservoir){
  u_int32_t bytes = 0;
  bytes += 5*sizeof(int); //the 5 ints in the struct
  bytes += 13*sizeof(double); //the 13 doubles in the struct
  bytes += reservoir->n_soil_layers*sizeof(double); //bytes in soil_layer_depths_m
  bytes += reservoir->n_soil_layers*sizeof(double); //bytes in smc_profile
  bytes += (reservoir->n_soil_layers+1)*sizeof(double); //bytes in delta_soil_layer_depth_m
  return bytes;
}

static void serialize_reservoir(struct conceptual_reservoir* reservoir, char* buffer){
  if(reservoir != NULL && buffer != NULL){
    char* pos = buffer;
    //serialize is_exponential
    pos = copy_to_int(&(reservoir->is_exponential), pos);
    //serialize is_sft_coupled
    pos = copy_to_int(&(reservoir->is_sft_coupled), pos);
    //serialize n_soil_layers
    pos = copy_to_int(&(reservoir->n_soil_layers), pos);
    //serialize aet_root_zone
    pos = copy_to_int(&(reservoir->aet_root_zone), pos);
    //serialize max_root_zone_layer
    pos = copy_to_int(&(reservoir->max_root_zone_layer), pos);

    //serialize gw_storage
    pos = copy_to_double(&(reservoir->gw_storage), pos);
    //serialize storage_max_m
    pos = copy_to_double(&(reservoir->storage_max_m), pos);
    //serialize storage_m
    pos = copy_to_double(&(reservoir->storage_m), pos);
    //serialize storage_change_m
    pos = copy_to_double(&(reservoir->storage_change_m), pos);
    //serialize coeff_primary
    pos = copy_to_double(&(reservoir->coeff_primary), pos);
    //serialize exponent_primary
    pos = copy_to_double(&(reservoir->exponent_primary), pos);
    //serialize storage_threshold_primary_m
    pos = copy_to_double(&(reservoir->storage_threshold_primary_m), pos);
    //serialize storage_threshold_secondary_m
    pos = copy_to_double(&(reservoir->storage_threshold_secondary_m), pos);
    //serialize coeff_secondary
    pos = copy_to_double(&(reservoir->coeff_secondary), pos);
    //serialize exponent_secondary
    pos = copy_to_double(&(reservoir->exponent_secondary), pos);
    //serialize ice_fraction_schaake
    pos = copy_to_double(&(reservoir->ice_fraction_schaake), pos);
    //serialize ice_fraction_xinanjiang
    pos = copy_to_double(&(reservoir->ice_fraction_xinanjiang), pos);
    //serialize soil_water_content_field_capacity
    pos = copy_to_double(&(reservoir->soil_water_content_field_capacity), pos);

    //Handle the dynamic parts
    //smc_profile
    pos = copy_to_double_array(reservoir->smc_profile, pos, reservoir->n_soil_layers);
    //soil_layer_depths_m
    pos = copy_to_double_array(reservoir->soil_layer_depths_m, pos, reservoir->n_soil_layers);
    //soil_layer_depths_m
    pos = copy_to_double_array(reservoir->soil_layer_depths_m, pos, reservoir->n_soil_layers);

  }
}

static void deserialize_reservoir(struct conceptual_reservoir* reservoir, char* buffer)
{
  if(reservoir != NULL && buffer != NULL){
    char* pos = buffer;
    //serialize is_exponential
    pos = copy_from_int(pos, &(reservoir->is_exponential));
    //serialize is_sft_coupled
    pos = copy_from_int(pos, &(reservoir->is_sft_coupled));
    //serialize n_soil_layers
    pos = copy_from_int(pos, &(reservoir->n_soil_layers));
    //serialize aet_root_zone
    pos = copy_from_int(pos, &(reservoir->aet_root_zone));
    //serialize max_root_zone_layer
    pos = copy_from_int(pos, &(reservoir->max_root_zone_layer));

    //serialize gw_storage
    pos = copy_from_double(pos, &(reservoir->gw_storage));
    //serialize storage_max_m
    pos = copy_from_double(pos, &(reservoir->storage_max_m));
    //serialize storage_m
    pos = copy_from_double(pos, &(reservoir->storage_m));
    //serialize storage_change_m
    pos = copy_from_double(pos, &(reservoir->storage_change_m));
    //serialize coeff_primary
    pos = copy_from_double(pos, &(reservoir->coeff_primary));
    //serialize exponent_primary
    pos = copy_from_double(pos, &(reservoir->exponent_primary));
    //serialize storage_threshold_primary_m
    pos = copy_from_double(pos, &(reservoir->storage_threshold_primary_m));
    //serialize storage_threshold_secondary_m
    pos = copy_from_double(pos, &(reservoir->storage_threshold_secondary_m));
    //serialize coeff_secondary
    pos = copy_from_double(pos, &(reservoir->coeff_secondary));
    //serialize exponent_secondary
    pos = copy_from_double(pos, &(reservoir->exponent_secondary));
    //serialize ice_fraction_schaake
    pos = copy_from_double(pos, &(reservoir->ice_fraction_schaake));
    //serialize ice_fraction_xinanjiang
    pos = copy_from_double(pos, &(reservoir->ice_fraction_xinanjiang));
    //serialize soil_water_content_field_capacity
    pos = copy_from_double(pos, &(reservoir->soil_water_content_field_capacity));

    //Handle the dynamic parts
    //smc_profile
    pos = copy_from_double_array(pos, reservoir->smc_profile, reservoir->n_soil_layers);
    //soil_layer_depths_m
    pos = copy_from_double_array(pos, reservoir->soil_layer_depths_m, reservoir->n_soil_layers);
    //soil_layer_depths_m
    pos = copy_from_double_array(pos, reservoir->soil_layer_depths_m, reservoir->n_soil_layers);
  }
}

extern void conceptual_reservoir_flux_calc(struct conceptual_reservoir *da_reservoir,
                                           double *primary_flux_m, double *secondary_flux_m);

#endif
