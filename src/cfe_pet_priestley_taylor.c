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


// A SIMPLE IMPLEMENTATION OF PRIESTLEY-TAYLOR ET FOR TESTING PURPOSES
// Uses OARC met variables (incoming shortwave and longwave radiation)
// to calculate radiation balance et, then applies P-T method.  Note
// this function is NOT intended ffor operational water prediction.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
//                  ################################ 
//                  ###FOR TESTING PURPOSES ONLY.###
//                  ################################ 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This function simplistically assumes constant land-surface temperature,
// ignores the influence of canopy, and assumes a constan atmospheric 
// longwave emmissivity, which it isn't.  FLO 9/2025

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cfe_pet_priestley_taylor.h"

// Priestley-Taylor potential evapotranspiration calculation
// Returns PET in meters per timestep.  Included by FLO for testing purposes.
// Priestley-Taylor potential evapotranspiration calculation
// Returns PET in meters per timestep
//###################################
double calculate_pet_priestley_taylor(const cfe_forcing_struct* forcing, int dt_seconds, double alpha_pt)
{
    if (forcing == NULL || dt_seconds <= 0) return 0.0;
    
    // Constants
    // const double ALPHA_PT = 1.26;          // Priestley-Taylor coefficient (dimensionless) from config file.
    const double STEFAN_BOLTZMANN = 5.67e-8; // Stefan-Boltzmann constant (W/m^2/K^4)
    const double LATENT_HEAT_VAPORIZATION = 2.45e6; // J/kg at 20C
    const double PSYCHROMETRIC_CONSTANT = 0.0665;   // kPa/K (approximate at sea level)
    const double LIQUID_WATER_DENSITY = 998.0;
    
    // Extract values from forcing structure
    double temp_K = forcing->TMP_2maboveground;      // Temperature in Kelvin
    
    // Calculate net radiation properly
    double surface_temperature_K = 282.0;   // THIS IS A MAJOR (bad) ASSUMPTION.  LAND SURFACE TEMPERATURE IS NOT CONSTANT.
    double surface_albedo = 0.25;
    double surface_longwave_rad_W_m_2 = STEFAN_BOLTZMANN * pow(surface_temperature_K, 4.0);
    double atmos_emissivity = 0.757; // assumed
    double downwelling_longwave = atmos_emissivity *  STEFAN_BOLTZMANN * pow(forcing->TMP_2maboveground,4.0);
    double net_radiation = (1.0 - surface_albedo) * forcing->DSWRF_surface + downwelling_longwave - surface_longwave_rad_W_m_2;
    
    // Convert temperature to Celsius for some calculations
    double temp_C = temp_K - 273.15;
    
    // Calculate saturation vapor pressure (Tetens equation) in kPa
    if (temp_C <= -237.3) {
        fprintf(stderr, "ERROR: air_temperature_C <= -237.3 C in calculate_pet_priestley_taylor().\n");
        exit(EXIT_FAILURE);
    }
    double es_kPa = 0.6108 * exp((17.27 * temp_C) / (temp_C + 237.3));

    // Calculate slope of saturation vapor pressure curve (kPa/K)
    double delta = (4098.0 * es_kPa) / pow(temp_C + 237.3, 2.0);
    
    // For Priestley-Taylor, we need net radiation minus soil heat flux
    // Assume soil heat flux is 10% of net radiation (typical approximation)
    double available_energy = net_radiation * 0.9; // W/m^2
    
    // Handle negative available energy (nighttime/winter)
    if (available_energy <= 0.0) {
        return 0.0; // No evapotranspiration
    }
    
    // Priestley-Taylor equation: PET = alpha * (Delta/(Delta+gamma)) * (Rn-G) / lambda
    // Where:
    // alpha = Priestley-Taylor coefficient (1.26)
    // Delta = slope of saturation vapor pressure curve (kPa/K)
    // gamma = psychrometric constant (kPa/K)
    // Rn-G = available energy (W/m^2)
    // lambda = latent heat of vaporization (J/kg)
    
    double pet_rate_m_per_s = alpha_pt * 
                              (delta / (delta + PSYCHROMETRIC_CONSTANT)) * 
                              (available_energy / LATENT_HEAT_VAPORIZATION) / LIQUID_WATER_DENSITY;
    
    // Convert from m/s to m per timestep
    double calc_pet_m_per_timestep = pet_rate_m_per_s * dt_seconds;
    
    // Sanity check: limit to reasonable values (max ~10mm/day = 0.01m/day)
    double max_daily_pet_m_per_day = 0.01; // 10mm/day in meters
    double hours_sunshine_per_day =12.0;   // assumed
    double max_pet_m_per_h = max_daily_pet_m_per_day / hours_sunshine_per_day;
    double max_pet_this_timestep = max_pet_m_per_h;
    
    if (calc_pet_m_per_timestep > max_pet_this_timestep) {
        calc_pet_m_per_timestep = max_pet_this_timestep;
    }
    
    return calc_pet_m_per_timestep;
}
