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
 

// ======================================================================
// Stateless soil step (BMI kernel) NDISC generic
// Positive vertical flux is downward.
// Units:
//   theta: m3/m3
//   psi:   m
//   K:     m/h
//   rainfall, PET inputs: mm/h (converted to m/h inside)
//   step-integrated totals: m
//   rates: m/h
//
// Author: Fred L. Ogden, July, 2025, NOAA/NWS Office of Water Prediction
//
// This code solves the soil moisture evolution given rainfall input with
// percolation output to groundwater plus output to a lateral flow routine
// from a homogeneous 2 m thick soil, discretized into 0.1, 0.3, 0.6, and 
// 1.0 m thick discretizations (discs) (top-down), as in Noah-MP.  It uses
// Darcy-Buckingham flux calculations and the arithmetic average of the
// unsaturated hydraulic conductivity in the discs on either side of the
// interface between them.   Like Noah-MP is uses the field capacity as
// the threshold to activate the percolation and lateral flow fluxes.  It
// extracts AET from the wettest root zone disc.  It uses substeps to keep
// the solution stable, based on a fraction of the available pore space 
// filled by rainfall during a sub-time-step to less than 20%, or the 
// distance that the Darcy flux moves in a sub-time-step to be less than
// 10 percent of the disc thickness. These ratios were determined by trial
// and error, and may not be optimal.  The code can optionally use look-
// up tables of pre-calculated Clapp-Hornberger K(theta) and psi(theta)
// consisting of 5 values (defined in cfe_soil_discrete.h) to reduce 
// computation.  The Clapp-Hornberger functions are linear after logarithmic
// transform so the fit is fantastic, even with only 5 points.
// ======================================================================

#include <math.h>
#include "cfe.h"
#include "discrete_soil_moisture.h"
#include "soil_helpers.h"

#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a < _b ? _a : _b; })

//####################################
static inline double storage_sum_ndisc(const double *theta, const double *dz_m)
{
    double s = 0.0;
    for (int i = 0; i < NDISC; i++) s += theta[i] * dz_m[i];
    return s;
}

//##############################
static inline int any_disc_above(const double *theta, double thresh, int ndisc)
{
    for (int i = 0; i < ndisc; i++) if (theta[i] > thresh) return 1;
    return 0;
}

// Pick a conservative substep count based on initial fluxes and rainfall demand.
// Generic across NDISC.  Note: discretization is abbreviated here as disc.
//
//
//#############################
static int choose_n_substeps_generic(double dt_hours,
                                double rain_mm_per_h,
                                const double *dz_m,
                                const double *theta, double theta_sat,
                                const double *q0_m_per_h, // [0..NDISC-2]
                                int ndisc)
{
    int nintf = ndisc - 1;  // the number of interfaces between discs

    // Max interface magnitude
    double qmax = 0.0;
    for (int i = 0; i < nintf; i++) {
        double a = fabs(q0_m_per_h[i]);
        if (a > qmax) qmax = a;
    }

    // Thinnest discretization
    double dzmin = dz_m[0];  // Uses Noah-MP discs of 0.1, 0.3, 0.6, 1.0 m

    // Flux criterion: how much of the water would move in dt_hours?
    double move_potential = qmax * dt_hours;                 // m
    double flux_ratio     = move_potential / (0.10 * dzmin);

    // Rain criterion: fraction of top-cell storage capacity asked for this hour
    double rain_rate_m_per_h = rain_mm_per_h / 1000.0;
    double rain_hour_m       = rain_rate_m_per_h * dt_hours;
    double cap1_m            = (theta_sat - theta[0]) * dz_m[0];
    if (cap1_m < 1e-12) cap1_m = 1e-12;
    double rain_ratio        = rain_hour_m / (0.20 * cap1_m);

    double severity = 0.0;
    if (flux_ratio > rain_ratio) { 
        severity = flux_ratio;
    } else {
        severity = rain_ratio;
    }

    int n_substeps = 1;
    if (severity <= 1.0 && rain_mm_per_h > 0.0) n_substeps = 2;
    else if (severity <= 2.0)  n_substeps = 4;
    else if (severity <= 3.0)  n_substeps = 6;
    else if (severity <= 5.0)  n_substeps = 8;
    else                       n_substeps = 12;

    if (n_substeps < 1)  n_substeps = 1;
    if (n_substeps > 12) n_substeps = 12;
    return n_substeps;
}

//##############################
int DSBM_step_one_hour_stateless(
    double                   available_gw_storage_m,
    const SoilControl        *control,
    const SoilGeometry       *geom,
    const SoilParameters     *params,
    const SoilLookupTables   *lut,                     // may be NULL
    const SoilStateIn        *state_in,
    struct EVAPOTRANSPIRATION_STRUCTURE *evap_struct,  // Pass the actual ET struct
    const SoilForcing        *forcing,
    SoilStateOut             *state_out,
    SoilFluxes               *flux,
    TimestepSoilVolbal       *volbal,
    FILE                     *debug_fptr)              // may be NULL too
{
    (void)debug_fptr;

    available_gw_storage_m = fmax(0.0, available_gw_storage_m);
    
    // ---- local working state ------------------------------------------------
    double theta[NDISC];
    for (int i = 0; i < NDISC; i++) theta[i] = state_in->theta_in[i];

    // outputs zeroed
    for (int i = 0; i < NDISC; i++) {
        flux->AET_by_disc_m[i] = 0.0;
        flux->lateral_by_disc_m[i] = 0.0;
        flux->interface_vol_m[i] = 0.0;
        flux->interface_rate_m_per_h[i] = 0.0;
        state_out->ch_lut_hint_out[i] = state_in->ch_lut_hint_in[i]; // start with input hints
    }
    flux->percolation_to_gw_m = 0.0;
    flux->rain_into_soil_m    = 0.0;
    flux->rain_excess_m       = 0.0;
    flux->n_substeps_used          = 0;

    volbal->in_rain_m       = 0.0;
    volbal->excess_m        = 0.0;
    volbal->perc_m          = 0.0;
    volbal->AET_m           = 0.0;
    volbal->lateral_m       = 0.0;
    volbal->delta_storage_m = 0.0;
    volbal->residual_m      = 0.0;

    const int ndisc  = control->ndisc;
    const int nintf  = ndisc - 1;
    const double delta_t_h = control->dt_hours;

    const double rain_rate_m_per_h = forcing->rain_mm_per_h / 1000.0;

    const double theta_floor = fmax(THETA_MIN, params->theta_r);

    // initial storage
    const double storage_start = storage_sum_ndisc(theta, geom->dz_m);

    //-- NEW
    // Create a working copy of the input soil state, and populate theta in discs

    SoilStateIn temp_soil_state = *state_in;

    for (int i = 0; i < NDISC; i++) {
        temp_soil_state.theta_in[i] = theta[i];
    }

    // Call ET calculation once ffor the full timestep, before substep loop, remove water from wettest disc
    et_from_soil_discrete(control, geom, params, &temp_soil_state, evap_struct); //calc. AET for this timestep

    // Update theta array with post-ET values
    for (int i = 0; i < NDISC; i++) {
        theta[i] = temp_soil_state.theta_in[i];
    }

    // Track ET removal by disc for flux accounting
    double et_removed = 0.0;
    for (int i = 0; i < NDISC; i++) {
        et_removed += (state_in->theta_in[i] - theta[i]) * geom->dz_m[i]; //FIXME- could just return value...
        flux->AET_by_disc_m[i] = et_removed;
    }

    volbal->AET_m = evap_struct->actual_et_from_soil_m_per_timestep;  // was updated in et_from_soil_discrete()
    //-- END NEW

    // Calculate initial fluxes at disc interfaces
    double q0[NDISC-1];  
    for (int i = 0; i < nintf; i++) {
        // calculate the Darcy-Buckingham (DB) flux from disc 0-1, disc 1-2, disc 2-3.
        q0[i] = flux_DB_pair(state_in->psi_in[i], state_in->K_in[i],
                             state_in->psi_in[i+1], state_in->K_in[i+1],
                             geom->dz_m[i], geom->dz_m[i+1]);
    }

    // Determine the number of substeps (needed in case of very wet soils in disc1)
    int n_substeps = choose_n_substeps_generic(delta_t_h, forcing->rain_mm_per_h,
                                     geom->dz_m, theta, params->theta_sat, q0, ndisc);
    if (n_substeps < 1) n_substeps = 1;
    flux->n_substeps_used = n_substeps;

    const double dt_sub = delta_t_h / (double)n_substeps;

    // per-substep work arrays
    double psi[NDISC], K[NDISC];
    double store_cap[NDISC];
    double pot_downflux[NDISC > 1 ? NDISC-1 : 1];  // >=0
    double V_if[NDISC > 1 ? NDISC - 1 : 1];          // signed desired
    double Accept[NDISC + 1];                      // [0..ndisc], nd = bottom
    for (int i = 0; i <= ndisc; i++) Accept[i] = 0.0;

    // External inflow (incident) is available to the caller via forcing and delta_t_h.
    // Here we only track infiltrated vs excess ffor the step.
    // Loop over substeps

    //  <----------------------------------------------------------- Start of substep loop
    for (int substep = 0; substep < n_substeps; substep++) {
        // See if partitioned soil moisture from Schaake/Xinanjiang fits into disc 0
        double rain_sub = rain_rate_m_per_h * dt_sub;  // m
        double cap1     = (params->theta_sat - theta[0]) * geom->dz_m[0];
        if (cap1 < 0.0) cap1 = 0.0;

        double used   = rain_sub;
        if (used > cap1) used = cap1;

        double excess = rain_sub - used;
        if (excess < 0.0) excess = 0.0;

        theta[0] += used / geom->dz_m[0];
        if (theta[0] > params->theta_sat) theta[0] = params->theta_sat;

        flux->rain_into_soil_m += used;
        flux->rain_excess_m    += excess;

        // Properties after rainfall addition (LUT or analytic)
        int hint_local[NDISC];
        for (int i = 0; i < NDISC; i++) hint_local[i] = state_out->ch_lut_hint_out[i];

        compute_props_with_option_stateless(
            (control->use_ch_lookup_table ? lut : NULL),
            theta, psi, K,
            params->theta_r, params->theta_sat,
            params->K_sat_cm_per_h, params->phi_sat_cm, params->b_exp,
            hint_local);

        // keep updated hints
        for (int i = 0; i < NDISC; i++) state_out->ch_lut_hint_out[i] = hint_local[i];

        // Interface fluxes and desired substep volumes
        for (int i = 0; i < nintf; i++) {
            const double q = flux_DB_pair(psi[i], K[i], psi[i+1], K[i+1],
                                          geom->dz_m[i], geom->dz_m[i+1]);
            V_if[i] = q * dt_sub;
            pot_downflux[i] = (q > 0.0) ? (q * dt_sub) : 0.0;
        }

        // Per-disc free storage up to saturation
        for (int d = 0; d < ndisc; d++) {
            double s = (params->theta_sat - theta[d]) * geom->dz_m[d];
            if (s < 0.0) s = 0.0;
            store_cap[d] = s;
        }

       // Modified to consider the situation where the available storage in the groundwater reservoior
       // is insufficient to accept all the percolation this time step. 
       // Key changes marked with 

        //  bottom potential percolation (K(theta4) * limiter) | GW storage limit
        //  This is precisely how they doo it in Noah-MP
        
        double bottom_potential = 0.0;
        if (theta[ndisc-1] > params->theta_fc) {
            double K_now = K_from_theta(theta[ndisc-1], params->theta_sat,
                                        params->K_sat_cm_per_h, params->b_exp);
            double percolation_rate = params->perc_limiter_0_to_1 * K_now;     // m/h
            if (percolation_rate > 0.0) {
                bottom_potential = percolation_rate * dt_sub;
                
                // Limited by available groundwater storage
                // More conservative: limit each substep to remaining storage / remaining substeps
                double remaining_gw_storage = available_gw_storage_m - flux->percolation_to_gw_m;
                double remaining_gw_storage_sub = remaining_gw_storage / (double)(n_substeps - substep);
                
                if (remaining_gw_storage_sub < bottom_potential) {
                    bottom_potential = fmax(0.0, remaining_gw_storage_sub);
                }
            }
        }

        // Downstream acceptance (bottom up) 
        // Accept[ndisc] = last store + bottom_potential (now GW-storage-limited)
        Accept[ndisc] = store_cap[ndisc-1] + bottom_potential;

        for (int i = nintf-1; i >= 0; i--) {
            double pass = pot_downflux[i];
            int down_index = i + 2;                 // downstream Accept slot; last interface -> ndisc
            if (down_index >= ndisc) down_index = ndisc;
            if (pass > Accept[down_index]) pass = Accept[down_index];
            Accept[i+1] = store_cap[i] + pass;
        }

        // Apply capped transfers across all interfaces
        for (int i = 0; i < nintf; i++) {
            double V = V_if[i];

            double accept_down =
                (i + 2 <= ndisc-1) ? Accept[i+2] : Accept[ndisc];

            if (V > 0.0) {
                double donor_avail = (theta[i] - theta_floor) * geom->dz_m[i];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i+1];

                double chain_pass = pot_downflux[i];
                if (chain_pass > accept_down) chain_pass = accept_down;

                double max_out = store_cap[i] + chain_pass;

                if (V > max_out) V = max_out;
                if (V > donor_avail) V = donor_avail;
                if (V > recv_space + chain_pass) V = recv_space + chain_pass;
                if (V < 0.0) V = 0.0;

                theta[i]   -= V / geom->dz_m[i];
                theta[i+1] += V / geom->dz_m[i+1];

                if (theta[i]   < theta_floor)      theta[i]   = theta_floor;
                if (theta[i+1] > params->theta_sat)   theta[i+1] = params->theta_sat;
            }
            else if (V < 0.0) {
                double need = -V;

                double donor_avail = (theta[i+1] - theta_floor) * geom->dz_m[i+1];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i];

                double move = need;
                if (move > donor_avail) move = donor_avail;
                if (move > recv_space)  move = recv_space;

                V = -move;

                theta[i+1] -= move / geom->dz_m[i+1];
                theta[i]   += move / geom->dz_m[i];

                if (theta[i+1] < theta_floor)      theta[i+1] = theta_floor;
                if (theta[i]   > params->theta_sat)   theta[i]   = params->theta_sat;
            }

            // accumulate internal interface volume
            flux->interface_vol_m[i] += V;
        }

        // Apply bottom percolation 
        double perc_vol = 0.0;
        if (theta[ndisc-1] > params->theta_fc && bottom_potential > 0.0) {
            double avail = (theta[ndisc-1] - params->theta_fc) * geom->dz_m[ndisc-1];
            if (avail < 0.0) avail = 0.0;

            perc_vol = bottom_potential;  // This is already GW-storage-limited
            if (perc_vol > avail) perc_vol = avail;

            // Additional safety check against total GW storage
            double total_perc_after = flux->percolation_to_gw_m + perc_vol;
            if (total_perc_after > available_gw_storage_m) {
                perc_vol = available_gw_storage_m - flux->percolation_to_gw_m;
                perc_vol = fmax(0.0, perc_vol);
            }

            if (perc_vol > 0.0) {
                theta[ndisc-1] -= perc_vol / geom->dz_m[ndisc-1];
                if (theta[ndisc-1] < params->theta_fc) theta[ndisc-1] = params->theta_fc;
            }
        }
        flux->percolation_to_gw_m += perc_vol;
        flux->interface_vol_m[ndisc-1] += perc_vol;


        // Remove calculated lateral flow to subsurface Nash cascade but only iff any disc > FC
        if (params->klf_per_h > 0.0 && any_disc_above(theta, params->theta_fc, ndisc)) {
            double lat_removed =
                remove_lateral_to_subsurface_nash_substep(
                    theta, geom->dz_m, params->theta_fc, params->theta_sat,
                    params->klf_per_h, dt_sub, flux->lateral_by_disc_m);

            volbal->lateral_m += lat_removed;
        }


        // safety clamp
        for (int i = 0; i < ndisc; i++) {
            if (theta[i] < theta_floor)    theta[i] = theta_floor;
            if (theta[i] > params->theta_sat) theta[i] = params->theta_sat;
        }
    } // <------------------------------------------------------------------end substep loop

    // finalize rates
    state_out->total_storage_m = 0.0;
    for (int i = 0; i < NDISC; i++) {
        flux->interface_rate_m_per_h[i] = flux->interface_vol_m[i] / delta_t_h;
    }

    // write state out
    for (int i = 0; i < NDISC; i++) state_out->theta_out[i] = theta[i];
    state_out->total_storage_m = storage_sum_ndisc(theta, geom->dz_m); 
    state_out->storage_deficit_m = geom->depth_m * params->theta_sat - state_out->total_storage_m;
    
    // volume balances
    double storage_end = storage_sum_ndisc(theta, geom->dz_m);

    volbal->in_rain_m  = flux->rain_into_soil_m;
    volbal->excess_m   = flux->rain_excess_m;
    volbal->perc_m     = flux->percolation_to_gw_m;

    // total lateral already accumulated in volbal->lateral_m inside loop
    // Ensure it matches sum of per-disc laterals (defensive)
    double lat_sum = 0.0; 
    for (int i = 0; i < NDISC; i++) lat_sum += flux->lateral_by_disc_m[i];
    volbal->lateral_m = lat_sum;

    volbal->delta_storage_m = (storage_end - storage_start);

    volbal->residual_m = (volbal->in_rain_m)
                   - (volbal->perc_m + volbal->AET_m + volbal->lateral_m)
                   -  volbal->delta_storage_m;

    return 0;
}


//##############################################################
//############   ET FROM SOIL DISCRETE  ########################
//##############################################################
void et_from_soil_discrete
    (
    const SoilControl*    soil_control,
    const SoilGeometry*   soil_geometry, 
    const SoilParameters* soil_parameters, 
    SoilStateIn*    soil_state,
    struct EVAPOTRANSPIRATION_STRUCTURE* evap_struct
    )
{
    //    FLO August, 2025
    //    Take AET from wettest root zone discretization (disc) in a discretized soil
    //    using Budyko type function to limit PET iff wilting_point < theta < field_capacity
    //    Assumes AET=PET when theta >= field_capacity moisture content
    //            AET=0  when theta <= wilting point moisture content
    //            0<=AET/PET<=1 between WP and FC (linear)
    //    Returns AET and modified theta in wettest root zone discretization
   

    // assume that root zone disc with highest moisture content satisfies all AET 
    
    int    wettest_disc            = -1;
    double wettest_theta           = 0.0;
    double dz_m_wettest_disc         = 0.0;
    double PET                     = evap_struct->reduced_potential_et_m_per_timestep;
    double AET                     = 0.0;
    
    for(int i = 0; i < soil_control->deepest_root_disc; i++) {       // find wettest disc
      if(soil_state->theta_in[i] > wettest_theta) {
          wettest_theta = soil_state->theta_in[i];
          wettest_disc = i;
      }
    }
    dz_m_wettest_disc = soil_geometry->dz_m[wettest_disc];    
    double wettest_disc_storage_m = wettest_theta * dz_m_wettest_disc;
      
    if(wettest_theta <= soil_parameters->theta_wp) {
        AET    = 0.0;
    } else if(wettest_theta >= soil_parameters->theta_fc) {
        AET    = min(PET, wettest_disc_storage_m);
    } else {
        double Budyko_numerator   = wettest_theta - soil_parameters->theta_wp;
        double Budyko_denominator = soil_parameters->theta_fc - soil_parameters->theta_wp;
        double Budyko_multiplier  = Budyko_numerator / Budyko_denominator;       
        AET = min(PET * Budyko_multiplier, wettest_disc_storage_m);
    }

    // remove AET from storage in wettest disc:
    if(AET > 0.0) {
      double delta_theta = AET / dz_m_wettest_disc;          // because: delta_theta = delta_Vw/Vt
      soil_state->theta_in[wettest_disc] -= delta_theta;   // negative by definition
    }
    evap_struct->actual_et_from_soil_m_per_timestep = AET;
    return;
    
}

