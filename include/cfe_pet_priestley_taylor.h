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

#ifndef CFE_PET_H
#define CFE_PET_H

#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Calculate potential evapotranspiration using Priestley-Taylor method
 * 
 * Calculates PET in meters per timestep using the Priestley-Taylor equation.
 * This is an energy-based method that requires temperature and radiation data.
 * 
 * @param forcing Pointer to forcing data structure containing meteorological variables
 * @param dt_seconds Time step duration in seconds
 * @param alpha_pt Priestley-Taylor coefficient (typically 1.26, dimensionless)
 * @return PET in meters per timestep, or 0.0 on error/invalid input
 */
double calculate_pet_priestley_taylor(const cfe_forcing_struct* forcing, 
                                     int dt_seconds, 
                                     double alpha_pt);

#ifdef __cplusplus
}
#endif

#endif /* CFE_PET_H */
