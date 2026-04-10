/*
 * bmi_cfe.c — CFE v3 BMI implementation
 *
 * Based on cfe3-project/src/cfe_bmi.c by FLO, migrated into the v2
 * BMI function structure.  Model state lives in CFE_Model_Context
 * stored in self->data (accessed via the CONTEXT macro).
 *
 * Key differences from v2 bmi_cfe.c:
 *  - State is CFE_Model_Context* in self->data (was cfe_state_struct*)
 *  - Initialize delegates to cfe_context_create_from_config()
 *  - Update delegates to cfe_context_update()
 *  - Finalize delegates to cfe_context_destroy()
 *  - Variable names are v3 names (discharge_m, surface_runoff_m, etc.)
 *  - get_value_ptr supports all variables + ngen mass balance protocol
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "bmi.h"
#include "bmi_cfe.h"
#include "cfe_context.h"
#include "cfe_helpers.h"
#include "parser_helpers.h"
#include "ngen_utilities.h"

/* ------------------------------------------------------------------ */
/* Cast helper — extract CFE_Model_Context from the BMI data pointer  */
/* ------------------------------------------------------------------ */
#define CONTEXT(self) ((CFE_Model_Context*)(self)->data)

/* ================================================================== */
/*  Variable name tables                                               */
/* ================================================================== */

/* --- inputs --- */
static const char* input_var_names[] = {
    "rainfall_depth_m",
    "et_potential_m",
    "verbosity",
    "forcing_file_path"
};
static const int INPUT_VAR_NAME_COUNT = 4;

/* --- outputs --- */
static const char* output_var_names[] = {
    /* Primary fluxes (per timestep) */
    "discharge_m",                              /*  0 */
    "surface_runoff_m",                         /*  1 */
    "lateral_flow_m",                           /*  2 */
    "baseflow_m",                               /*  3 */
    "actual_et_m",                              /*  4 */
    "vol_balance_residual_m",                   /*  5 */

    /* State scalars (checkpointing) */
    "state_soil_storage_m",                     /*  6 */
    "state_gw_storage_m",                       /*  7 */
    "state_current_timestep",                   /*  8 */

    /* State arrays (checkpointing) */
    "state_soil_moisture_theta",                /*  9  — NDISC elements */
    "state_nash_surface_storage",               /* 10  — N_nash elements */
    "state_nash_subsurface_storage",            /* 11  — 2 elements */
    "state_giuh_queue",                         /* 12  — num_giuh elements */

    /* Config / parameters for interpretation */
    "config_simulate_discrete_soil_moisture",   /* 13 */
    "param_catchment_area_km2",                 /* 14 */
    "param_soil_depth_m",                       /* 15 */
    "param_soil_porosity",                      /* 16 */

    /* Per-timestep volume balance */
    "timestep_storage_start_m",                 /* 17 */
    "timestep_input_m",                         /* 18 */
    "timestep_output_m",                        /* 19 */
    "timestep_storage_end_m"                    /* 20 */
};
static const int OUTPUT_VAR_NAME_COUNT = 21;

/* --- calibration parameters (get_value / set_value / get_value_ptr) --- */
/* v3 canonical names are listed first; v2 aliases in comments */
static const char* param_var_names[] = {
    "soil_effective_porosity",                   /*  0  v2: maxsmc */
    "soil_saturated_hydraulic_conductivity",     /*  1  v2: satdk — m/s */
    "soil_percolation_rate_limiter",             /*  2  v2: slope — 0-1 */
    "soil_Clapp_Hornberger_b",                  /*  3  v2: b */
    "soil_lateral_flow_K",                       /*  4  v2: Klf — per h */
    "subsurface_nash_K",                         /*  5  v2: Kn — per h */
    "gw_discharge_coefficient",                  /*  6  v2: Cgw — m/s */
    "gw_discharge_exponent",                     /*  7  v2: expon */
    "gw_max_storage_m",                          /*  8  v2: max_gw_storage */
    "soil_saturated_capillary_head",             /*  9  v2: satpsi — m */
    "soil_wilting_point",                        /* 10  v2: wltsmc */
    "soil_field_capacity_fraction",              /* 11  v2: alpha_fc — Pcap/Patm */
    "refkdt",                                   /* 12  v2: refkdt */
    "Xinanjiang_inflection_a",                   /* 13  v2: a_Xinanjiang_inflection_point_parameter */
    "Xinanjiang_shape_b",                        /* 14  v2: b_Xinanjiang_shape_parameter */
    "Xinanjiang_shape_x",                        /* 15  v2: x_Xinanjiang_shape_parameter */
    "surface_nash_Kinf",                         /* 16  v2: Kinf_nash_surface — per h */
    "surface_nash_retention_depth_m",            /* 17  v2: retention_depth_nash_surface */
    "Priestley_Taylor_alpha",                    /* 18  (new in v3) */
    "soil_ice_imperv_threshold"                  /* 19  (new in v3) */
};
static const int PARAM_VAR_NAME_COUNT = 20;

/* Resolve a parameter name (v3 or v2 alias) to a pointer into ctx->parameters.
 * Returns NULL if unrecognized. */
static double* param_field_ptr(CFE_Model_Context *ctx, const char *name) {
    cfe_parameters_struct *p = &ctx->parameters;

    if (strcmp(name, "soil_effective_porosity") == 0               || strcmp(name, "maxsmc") == 0)
        return &p->effective_porosity;
    if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0 || strcmp(name, "satdk") == 0)
        return &p->ksat_m_per_s;
    if (strcmp(name, "soil_percolation_rate_limiter") == 0         || strcmp(name, "slope") == 0)
        return &p->soil_to_gw_percolation_rate_limiter_0_1;
    if (strcmp(name, "soil_Clapp_Hornberger_b") == 0              || strcmp(name, "b") == 0)
        return &p->soil_b;
    if (strcmp(name, "soil_lateral_flow_K") == 0                   || strcmp(name, "Klf") == 0)
        return &p->soil_k_lateral_per_h;
    if (strcmp(name, "subsurface_nash_K") == 0                     || strcmp(name, "Kn") == 0)
        return &p->nash_subsurface_K_per_h;
    if (strcmp(name, "gw_discharge_coefficient") == 0              || strcmp(name, "Cgw") == 0)
        return &p->gw_discharge_coeff_m_per_s;
    if (strcmp(name, "gw_discharge_exponent") == 0                 || strcmp(name, "expon") == 0)
        return &p->gw_discharge_exponent;
    if (strcmp(name, "gw_max_storage_m") == 0                      || strcmp(name, "max_gw_storage") == 0)
        return &p->gw_max_storage_m;
    if (strcmp(name, "soil_saturated_capillary_head") == 0         || strcmp(name, "satpsi") == 0)
        return &p->sat_capillary_head_m;
    if (strcmp(name, "soil_wilting_point") == 0                    || strcmp(name, "wltsmc") == 0)
        return &p->wilting_point;
    if (strcmp(name, "soil_field_capacity_fraction") == 0          || strcmp(name, "alpha_fc") == 0)
        return &p->field_capacity_Pcap_over_Patm;
    if (strcmp(name, "refkdt") == 0)
        return &p->refkdt;
    if (strcmp(name, "Xinanjiang_inflection_a") == 0               || strcmp(name, "a_Xinanjiang_inflection_point_parameter") == 0)
        return &p->xj_tension_inflection_0_1;
    if (strcmp(name, "Xinanjiang_shape_b") == 0                    || strcmp(name, "b_Xinanjiang_shape_parameter") == 0)
        return &p->xj_tension_b;
    if (strcmp(name, "Xinanjiang_shape_x") == 0                    || strcmp(name, "x_Xinanjiang_shape_parameter") == 0)
        return &p->xj_free_b;
    if (strcmp(name, "surface_nash_Kinf") == 0                     || strcmp(name, "Kinf_nash_surface") == 0)
        return &p->surface_Kinf_per_h;
    if (strcmp(name, "surface_nash_retention_depth_m") == 0        || strcmp(name, "retention_depth_nash_surface") == 0)
        return &p->surface_retention_depth_m;
    if (strcmp(name, "Priestley_Taylor_alpha") == 0                || strcmp(name, "alpha_pt") == 0)
        return &p->alpha_pt;
    if (strcmp(name, "soil_ice_imperv_threshold") == 0)
        return &p->soil_ice_imperv_threshold;

    return NULL;
}

/* ================================================================== */
/*  Lifecycle: Initialize / Update / Finalize                          */
/* ================================================================== */

static int Initialize(Bmi *self, const char *cfg_file) {
    if (cfg_file == NULL || self == NULL) return BMI_FAILURE;

    CFE_Model_Context* ctx = NULL;
    int result = cfe_context_create_from_config(cfg_file, &ctx);
    if (result != 0 || ctx == NULL) return BMI_FAILURE;

    self->data = (void*)ctx;
    return BMI_SUCCESS;
}

static int Update(Bmi *self) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    CFE_Model_Context *ctx = CONTEXT(self);
    int result = cfe_context_update(ctx);
    if (result != 0) return BMI_FAILURE;

    /* cache the volume balance residual for get_value_ptr */
    ctx->vol_balance_residual_m = ctx->volbal.volstart + ctx->volbal.volin
                                - ctx->volbal.volout   - ctx->volbal.volend;
    return BMI_SUCCESS;
}

static int Update_until(Bmi *self, double then) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;

    double current_time = 0.0;
    int dt_seconds;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_seconds);
    int step;
    cfe_context_get_current_step(CONTEXT(self), &step);
    current_time = (double)step * (double)dt_seconds;

    while (current_time < then) {
        int result = Update(self);
        if (result != BMI_SUCCESS) return result;
        cfe_context_get_current_step(CONTEXT(self), &step);
        current_time = (double)step * (double)dt_seconds;
    }
    return BMI_SUCCESS;
}

static int Finalize(Bmi *self) {
    if (CONTEXT(self) != NULL) {
        cfe_context_destroy(CONTEXT(self));
        self->data = NULL;
    }
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Component / exchange item metadata                                 */
/* ================================================================== */

static int Get_component_name(Bmi *self, char *name) {
    strcpy(name, "CFE - Conceptual Functional Equivalent");
    return BMI_SUCCESS;
}

static int Get_input_item_count(Bmi *self, int *count) {
    *count = INPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}

static int Get_output_item_count(Bmi *self, int *count) {
    *count = OUTPUT_VAR_NAME_COUNT;
    return BMI_SUCCESS;
}

static int Get_input_var_names(Bmi *self, char **names) {
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++)
        strcpy(names[i], input_var_names[i]);
    return BMI_SUCCESS;
}

static int Get_output_var_names(Bmi *self, char **names) {
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++)
        strcpy(names[i], output_var_names[i]);
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Variable metadata                                                  */
/* ================================================================== */

static int Get_var_grid(Bmi *self, const char *name, int *grid) {
    if      (strcmp(name, "state_soil_moisture_theta")    == 0) *grid = 1;
    else if (strcmp(name, "state_nash_surface_storage")   == 0) *grid = 2;
    else if (strcmp(name, "state_nash_subsurface_storage")== 0) *grid = 3;
    else if (strcmp(name, "state_giuh_queue")             == 0) *grid = 4;
    else *grid = 0;
    return BMI_SUCCESS;
}

static int Get_var_type(Bmi *self, const char *name, char *type) {
    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0) {
        strcpy(type, "int");
        return BMI_SUCCESS;
    }
    if (strcmp(name, "forcing_file_path") == 0) {
        strcpy(type, "string");
        return BMI_SUCCESS;
    }
    /* Check output and input variable names (all double) */
    for (int i = 0; i < OUTPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, output_var_names[i]) == 0) { strcpy(type, "double"); return BMI_SUCCESS; }
    for (int i = 0; i < INPUT_VAR_NAME_COUNT; i++)
        if (strcmp(name, input_var_names[i]) == 0) { strcpy(type, "double"); return BMI_SUCCESS; }
    /* Check calibration parameter names (all double) */
    if (CONTEXT(self) && param_field_ptr(CONTEXT(self), name) != NULL) {
        strcpy(type, "double");
        return BMI_SUCCESS;
    }
    return BMI_FAILURE;
}

static int Get_var_units(Bmi *self, const char *name, char *units) {
    if (strcmp(name, "rainfall_depth_m")         == 0 ||
        strcmp(name, "et_potential_m")           == 0 ||
        strcmp(name, "discharge_m")              == 0 ||
        strcmp(name, "surface_runoff_m")         == 0 ||
        strcmp(name, "lateral_flow_m")           == 0 ||
        strcmp(name, "baseflow_m")               == 0 ||
        strcmp(name, "actual_et_m")              == 0 ||
        strcmp(name, "vol_balance_residual_m")   == 0 ||
        strcmp(name, "state_soil_storage_m")     == 0 ||
        strcmp(name, "state_gw_storage_m")       == 0 ||
        strcmp(name, "param_soil_depth_m")       == 0 ||
        strcmp(name, "timestep_storage_start_m") == 0 ||
        strcmp(name, "timestep_input_m")         == 0 ||
        strcmp(name, "timestep_output_m")        == 0 ||
        strcmp(name, "timestep_storage_end_m")   == 0) {
        strcpy(units, "m");
    }
    else if (strcmp(name, "state_soil_moisture_theta") == 0 ||
             strcmp(name, "param_soil_porosity") == 0) {
        strcpy(units, "-");
    }
    else if (strcmp(name, "state_nash_surface_storage")    == 0 ||
             strcmp(name, "state_nash_subsurface_storage") == 0 ||
             strcmp(name, "state_giuh_queue")              == 0) {
        strcpy(units, "m");
    }
    else if (strcmp(name, "param_catchment_area_km2") == 0) {
        strcpy(units, "km2");
    }
    else if (strcmp(name, "verbosity") == 0 ||
             strcmp(name, "state_current_timestep") == 0 ||
             strcmp(name, "config_simulate_discrete_soil_moisture") == 0 ||
             strcmp(name, "forcing_file_path") == 0) {
        strcpy(units, "1");
    }
    /* --- calibration parameter units ---
     * Units are the internal representation in cfe_parameters_struct.
     * These match the units documented in bmi_config_cat87_v3.cf3 and
     * the unit conversions applied in cfe_helpers.c map_config_to_parameters_and_options().
     * Note: v2 param names (satdk, Cgw, etc.) resolve to the same fields. */
    else if (strcmp(name, "soil_effective_porosity") == 0     || strcmp(name, "maxsmc") == 0 ||
             strcmp(name, "soil_wilting_point") == 0          || strcmp(name, "wltsmc") == 0 ||
             strcmp(name, "soil_field_capacity_fraction") == 0|| strcmp(name, "alpha_fc") == 0 ||
             strcmp(name, "soil_percolation_rate_limiter") == 0 || strcmp(name, "slope") == 0 ||
             strcmp(name, "Xinanjiang_inflection_a") == 0     || strcmp(name, "a_Xinanjiang_inflection_point_parameter") == 0 ||
             strcmp(name, "soil_ice_imperv_threshold") == 0) {
        strcpy(units, "-");  /* dimensionless fractions (V/V or 0-1) */
    }
    else if (strcmp(name, "soil_saturated_hydraulic_conductivity") == 0 || strcmp(name, "satdk") == 0 ||
             strcmp(name, "gw_discharge_coefficient") == 0               || strcmp(name, "Cgw") == 0) {
        strcpy(units, "m s-1");  /* stored as m/s in cfe_parameters_struct */
    }
    else if (strcmp(name, "soil_lateral_flow_K") == 0        || strcmp(name, "Klf") == 0 ||
             strcmp(name, "subsurface_nash_K") == 0           || strcmp(name, "Kn") == 0 ||
             strcmp(name, "surface_nash_Kinf") == 0           || strcmp(name, "Kinf_nash_surface") == 0) {
        strcpy(units, "h-1");  /* per-hour rate constants */
    }
    else if (strcmp(name, "soil_saturated_capillary_head") == 0 || strcmp(name, "satpsi") == 0 ||
             strcmp(name, "gw_max_storage_m") == 0               || strcmp(name, "max_gw_storage") == 0 ||
             strcmp(name, "surface_nash_retention_depth_m") == 0  || strcmp(name, "retention_depth_nash_surface") == 0) {
        strcpy(units, "m");  /* meters */
    }
    else if (strcmp(name, "soil_Clapp_Hornberger_b") == 0    || strcmp(name, "b") == 0 ||
             strcmp(name, "gw_discharge_exponent") == 0       || strcmp(name, "expon") == 0 ||
             strcmp(name, "Xinanjiang_shape_b") == 0           || strcmp(name, "b_Xinanjiang_shape_parameter") == 0 ||
             strcmp(name, "Xinanjiang_shape_x") == 0           || strcmp(name, "x_Xinanjiang_shape_parameter") == 0 ||
             strcmp(name, "refkdt") == 0 ||
             strcmp(name, "Priestley_Taylor_alpha") == 0       || strcmp(name, "alpha_pt") == 0) {
        strcpy(units, "-");  /* dimensionless exponents and coefficients */
    }
    else {
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_var_itemsize(Bmi *self, const char *name, int *size) {
    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0) {
        *size = sizeof(int);
    }
    else if (strcmp(name, "forcing_file_path") == 0) {
        *size = sizeof(char);
    }
    else {
        /* all remaining recognized variables (outputs, inputs, params) are double */
        char type[BMI_MAX_TYPE_NAME];
        if (Get_var_type(self, name, type) != BMI_SUCCESS) return BMI_FAILURE;
        if (strcmp(type, "double") == 0) *size = sizeof(double);
        else return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_var_nbytes(Bmi *self, const char *name, int *nbytes) {
    if (strcmp(name, "verbosity") == 0 ||
        strcmp(name, "state_current_timestep") == 0 ||
        strcmp(name, "config_simulate_discrete_soil_moisture") == 0) {
        *nbytes = sizeof(int);
    }
    else if (strcmp(name, "forcing_file_path") == 0) {
        *nbytes = PATH_FILENAME_STRING_LENGTH;
    }
    else if (strcmp(name, "state_soil_moisture_theta") == 0) {
        *nbytes = NDISC * sizeof(double);
    }
    else if (strcmp(name, "state_nash_surface_storage") == 0) {
        if (CONTEXT(self) && CONTEXT(self)->options.surface_routing_scheme == SURF_ROUTE_NASH_CASCADE)
            *nbytes = CONTEXT(self)->parameters.nash_surface_N * sizeof(double);
        else
            *nbytes = 0;
    }
    else if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        *nbytes = 2 * sizeof(double);
    }
    else if (strcmp(name, "state_giuh_queue") == 0) {
        if (CONTEXT(self) && CONTEXT(self)->options.surface_routing_scheme == SURF_ROUTE_GIUH)
            *nbytes = CONTEXT(self)->parameters.giuh_num_ordinates * sizeof(double);
        else
            *nbytes = 0;
    }
    else {
        /* remaining recognized variables (scalar outputs, inputs, params) are single doubles */
        int itemsize;
        if (Get_var_itemsize(self, name, &itemsize) != BMI_SUCCESS) return BMI_FAILURE;
        *nbytes = itemsize;
    }
    return BMI_SUCCESS;
}

static int Get_var_location(Bmi *self, const char *name, char *location) {
    strcpy(location, "node");
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Time information                                                   */
/* ================================================================== */

static int Get_current_time(Bmi *self, double *time) {
    if (CONTEXT(self) == NULL) { *time = 0.0; return BMI_SUCCESS; }
    int step;
    cfe_context_get_current_step(CONTEXT(self), &step);
    int dt_s;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
    *time = (double)step * (double)dt_s;
    return BMI_SUCCESS;
}

static int Get_start_time(Bmi *self, double *time) {
    *time = 0.0;
    return BMI_SUCCESS;
}

static int Get_end_time(Bmi *self, double *time) {
    if (CONTEXT(self) != NULL && CONTEXT(self)->options.num_timesteps > 0) {
        int dt_s;
        cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
        *time = (double)CONTEXT(self)->options.num_timesteps * (double)dt_s;
    } else {
        *time = -1.0;  /* unknown — forcings arrive via BMI */
    }
    return BMI_SUCCESS;
}

static int Get_time_units(Bmi *self, char *units) {
    strcpy(units, "s");
    return BMI_SUCCESS;
}

static int Get_time_step(Bmi *self, double *dt) {
    if (CONTEXT(self) == NULL) return BMI_FAILURE;
    int dt_s;
    cfe_context_get_time_step_seconds(CONTEXT(self), &dt_s);
    *dt = (double)dt_s;
    return BMI_SUCCESS;
}

/* ================================================================== */
/*  Get value / Get value ptr / Get value at indices                   */
/* ================================================================== */

static int Get_value(Bmi *self, const char *name, void *dest) {
    CFE_Model_Context *ctx = CONTEXT(self);
    if (ctx == NULL) return BMI_FAILURE;
    double *d = (double*)dest;

    /* --- primary fluxes --- */
    if (strcmp(name, "discharge_m") == 0) {
        *d = ctx->last_outputs.qout_m;
    }
    else if (strcmp(name, "surface_runoff_m") == 0) {
        *d = ctx->last_outputs.surface_runoff_generated_m;
    }
    else if (strcmp(name, "lateral_flow_m") == 0) {
        *d = ctx->last_outputs.lateral_flow_m;
    }
    else if (strcmp(name, "baseflow_m") == 0) {
        *d = ctx->last_outputs.baseflow_m;
    }
    else if (strcmp(name, "actual_et_m") == 0) {
        *d = ctx->last_outputs.actual_et_m;
    }
    else if (strcmp(name, "vol_balance_residual_m") == 0) {
        *d = ctx->vol_balance_residual_m;
    }

    /* --- state scalars --- */
    else if (strcmp(name, "state_soil_storage_m") == 0) {
        *d = ctx->state.soil_storage_m;
    }
    else if (strcmp(name, "state_gw_storage_m") == 0) {
        *d = ctx->state.gw_storage_m;
    }
    else if (strcmp(name, "state_current_timestep") == 0) {
        *(int*)dest = ctx->state.current_time_step;
    }

    /* --- state arrays --- */
    else if (strcmp(name, "state_soil_moisture_theta") == 0) {
        for (int i = 0; i < NDISC; i++)
            d[i] = ctx->state.soil_discrete_storage_theta[i];
    }
    else if (strcmp(name, "state_nash_surface_storage") == 0) {
        if (ctx->options.surface_routing_scheme != SURF_ROUTE_NASH_CASCADE)
            return BMI_FAILURE;
        for (int i = 0; i < ctx->parameters.nash_surface_N; i++)
            d[i] = ctx->state.nash_surface_storage_m[i];
    }
    else if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        for (int i = 0; i < 2; i++)
            d[i] = ctx->state.nash_subsurface_storage_m[i];
    }
    else if (strcmp(name, "state_giuh_queue") == 0) {
        if (ctx->options.surface_routing_scheme != SURF_ROUTE_GIUH)
            return BMI_FAILURE;
        for (int i = 0; i < ctx->parameters.giuh_num_ordinates; i++)
            d[i] = ctx->state.giuh_queue_m[i];
    }

    /* --- config / parameters --- */
    else if (strcmp(name, "config_simulate_discrete_soil_moisture") == 0) {
        *(int*)dest = ctx->options.simulate_discrete_soil_moisture;
    }
    else if (strcmp(name, "param_catchment_area_km2") == 0) {
        *d = ctx->parameters.catchment_area_km2;
    }
    else if (strcmp(name, "param_soil_depth_m") == 0) {
        *d = ctx->parameters.soil_depth_m;
    }
    else if (strcmp(name, "param_soil_porosity") == 0) {
        *d = ctx->parameters.effective_porosity;
    }

    /* --- per-timestep volume balance --- */
    else if (strcmp(name, "timestep_storage_start_m") == 0) {
        *d = ctx->timestep_storage_start_m;
    }
    else if (strcmp(name, "timestep_input_m") == 0) {
        *d = ctx->timestep_input_m;
    }
    else if (strcmp(name, "timestep_output_m") == 0) {
        *d = ctx->timestep_output_m;
    }
    else if (strcmp(name, "timestep_storage_end_m") == 0) {
        *d = ctx->timestep_storage_end_m;
    }

    /* --- inputs (read-back) --- */
    else if (strcmp(name, "rainfall_depth_m") == 0) {
        *d = ctx->forcing.rainfall_depth_m;
    }
    else if (strcmp(name, "et_potential_m") == 0) {
        *d = ctx->forcing.et_potential_m;
    }
    else if (strcmp(name, "verbosity") == 0) {
        *(int*)dest = ctx->options.verbosity;
    }
    else if (strcmp(name, "forcing_file_path") == 0) {
        strncpy((char*)dest, ctx->options.input_forcing_filename,
                PATH_FILENAME_STRING_LENGTH);
    }

    else {
        /* check calibration parameters */
        double *pp = param_field_ptr(ctx, name);
        if (pp != NULL) { *d = *pp; }
        else return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_value_ptr(Bmi *self, const char *name, void **dest) {
    CFE_Model_Context *ctx = CONTEXT(self);
    if (ctx == NULL) return BMI_FAILURE;

    /* --- ngen mass balance protocol --- */
    if (strcmp(name, NGEN_MASS_IN) == 0)     { *dest = &ctx->volbal.cumulative_vol;   return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_OUT) == 0)    { *dest = &ctx->volbal.volout;           return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_STORED) == 0) { *dest = &ctx->volbal.volume_in_domain; return BMI_SUCCESS; }
    if (strcmp(name, NGEN_MASS_LEAKED) == 0) { *dest = &ctx->volbal.leakage;          return BMI_SUCCESS; }

    /* --- output scalars --- */
    if (strcmp(name, "discharge_m") == 0)              { *dest = &ctx->last_outputs.qout_m;                    return BMI_SUCCESS; }
    if (strcmp(name, "surface_runoff_m") == 0)         { *dest = &ctx->last_outputs.surface_runoff_generated_m; return BMI_SUCCESS; }
    if (strcmp(name, "lateral_flow_m") == 0)           { *dest = &ctx->last_outputs.lateral_flow_m;            return BMI_SUCCESS; }
    if (strcmp(name, "baseflow_m") == 0)               { *dest = &ctx->last_outputs.baseflow_m;                return BMI_SUCCESS; }
    if (strcmp(name, "actual_et_m") == 0)              { *dest = &ctx->last_outputs.actual_et_m;               return BMI_SUCCESS; }
    if (strcmp(name, "state_soil_storage_m") == 0)     { *dest = &ctx->state.soil_storage_m;                   return BMI_SUCCESS; }
    if (strcmp(name, "state_gw_storage_m") == 0)       { *dest = &ctx->state.gw_storage_m;                     return BMI_SUCCESS; }
    if (strcmp(name, "state_current_timestep") == 0)   { *dest = &ctx->state.current_time_step;                return BMI_SUCCESS; }
    if (strcmp(name, "config_simulate_discrete_soil_moisture") == 0) { *dest = &ctx->options.simulate_discrete_soil_moisture; return BMI_SUCCESS; }
    if (strcmp(name, "param_catchment_area_km2") == 0) { *dest = &ctx->parameters.catchment_area_km2;         return BMI_SUCCESS; }
    if (strcmp(name, "param_soil_depth_m") == 0)       { *dest = &ctx->parameters.soil_depth_m;                return BMI_SUCCESS; }
    if (strcmp(name, "param_soil_porosity") == 0)      { *dest = &ctx->parameters.effective_porosity;          return BMI_SUCCESS; }
    if (strcmp(name, "timestep_storage_start_m") == 0) { *dest = &ctx->timestep_storage_start_m;               return BMI_SUCCESS; }
    if (strcmp(name, "timestep_input_m") == 0)         { *dest = &ctx->timestep_input_m;                       return BMI_SUCCESS; }
    if (strcmp(name, "timestep_output_m") == 0)        { *dest = &ctx->timestep_output_m;                      return BMI_SUCCESS; }
    if (strcmp(name, "timestep_storage_end_m") == 0)   { *dest = &ctx->timestep_storage_end_m;                 return BMI_SUCCESS; }
    if (strcmp(name, "vol_balance_residual_m") == 0) { *dest = &ctx->vol_balance_residual_m;                return BMI_SUCCESS; }

    /* --- output arrays --- */
    if (strcmp(name, "state_soil_moisture_theta") == 0)     { *dest = ctx->state.soil_discrete_storage_theta;  return BMI_SUCCESS; }
    if (strcmp(name, "state_nash_surface_storage") == 0)    { *dest = ctx->state.nash_surface_storage_m;       return BMI_SUCCESS; }
    if (strcmp(name, "state_nash_subsurface_storage") == 0) { *dest = ctx->state.nash_subsurface_storage_m;    return BMI_SUCCESS; }
    if (strcmp(name, "state_giuh_queue") == 0)              { *dest = ctx->state.giuh_queue_m;                 return BMI_SUCCESS; }

    /* --- inputs --- */
    if (strcmp(name, "rainfall_depth_m") == 0)  { *dest = &ctx->forcing.rainfall_depth_m;          return BMI_SUCCESS; }
    if (strcmp(name, "et_potential_m") == 0)    { *dest = &ctx->forcing.et_potential_m;             return BMI_SUCCESS; }
    if (strcmp(name, "verbosity") == 0)         { *dest = &ctx->options.verbosity;                  return BMI_SUCCESS; }
    if (strcmp(name, "forcing_file_path") == 0) { *dest = ctx->options.input_forcing_filename;      return BMI_SUCCESS; }

    /* --- calibration parameters --- */
    double *pp = param_field_ptr(ctx, name);
    if(pp != NULL) { 
        *dest = pp; 
        return BMI_SUCCESS; 
    }

    return BMI_FAILURE;
}

static int Get_value_at_indices(Bmi *self, const char *name, void *dest, int *inds, int count) {
    if (strcmp(name, "state_soil_moisture_theta") == 0 && CONTEXT(self)) {
        double *d = (double*)dest;
        for (int i = 0; i < count; i++) {
            if (inds[i] < 0 || inds[i] >= NDISC) return BMI_FAILURE;
            d[i] = CONTEXT(self)->state.soil_discrete_storage_theta[inds[i]];
        }
        return BMI_SUCCESS;
    }
    if (count == 1) return Get_value(self, name, dest);
    return BMI_FAILURE;
}

/* ================================================================== */
/*  Set value / Set value at indices                                   */
/* ================================================================== */

static int Set_value(Bmi *self, const char *name, void *src) {
    CFE_Model_Context *ctx = CONTEXT(self);
    if (ctx == NULL) return BMI_FAILURE;

    /* --- forcing inputs --- */
    if (strcmp(name, "rainfall_depth_m") == 0) {
        ctx->forcing.rainfall_depth_m = *(double*)src;
    }
    else if (strcmp(name, "et_potential_m") == 0) {
        ctx->forcing.et_potential_m = *(double*)src;
    }
    else if (strcmp(name, "verbosity") == 0) {
        ctx->options.verbosity = *(int*)src;
    }
    else if (strcmp(name, "forcing_file_path") == 0) {
        strncpy(ctx->options.input_forcing_filename, (char*)src,
                sizeof(ctx->options.input_forcing_filename) - 1);
        ctx->options.input_forcing_filename[sizeof(ctx->options.input_forcing_filename) - 1] = '\0';
    }

    /* --- state arrays (hotstart) --- */
    else if (strcmp(name, "state_soil_moisture_theta") == 0) {
        double *s = (double*)src;
        for (int i = 0; i < NDISC; i++)
            ctx->state.soil_discrete_storage_theta[i] = s[i];
    }
    else if (strcmp(name, "state_nash_surface_storage") == 0) {
        if (ctx->options.surface_routing_scheme != SURF_ROUTE_NASH_CASCADE)
            return BMI_FAILURE;
        double *s = (double*)src;
        for (int i = 0; i < ctx->parameters.nash_surface_N; i++)
            ctx->state.nash_surface_storage_m[i] = s[i];
    }
    else if (strcmp(name, "state_nash_subsurface_storage") == 0) {
        double *s = (double*)src;
        for (int i = 0; i < 2; i++)
            ctx->state.nash_subsurface_storage_m[i] = s[i];
    }
    else if (strcmp(name, "state_giuh_queue") == 0) {
        if (ctx->options.surface_routing_scheme != SURF_ROUTE_GIUH)
            return BMI_FAILURE;
        double *s = (double*)src;
        for (int i = 0; i < ctx->parameters.giuh_num_ordinates; i++)
            ctx->state.giuh_queue_m[i] = s[i];
    }

    /* --- state scalars (hotstart) --- */
    else if (strcmp(name, "state_soil_storage_m") == 0) {
        ctx->state.soil_storage_m = *(double*)src;
    }
    else if (strcmp(name, "state_gw_storage_m") == 0) {
        ctx->state.gw_storage_m = *(double*)src;
    }
    else if (strcmp(name, "state_current_timestep") == 0) {
        ctx->state.current_time_step = *(int*)src;
    }

    else {
        /* calibration parameters — all are double type */
        double *pp = param_field_ptr(ctx, name);
        if (pp != NULL) {
            *pp = *(double*)src;
        } else {
            return BMI_FAILURE;
        }
    }
    return BMI_SUCCESS;
}

static int Set_value_at_indices(Bmi *self, const char *name, int *inds, int count, void *src) {
    if (strcmp(name, "state_soil_moisture_theta") == 0 && CONTEXT(self)) {
        double *s = (double*)src;
        for (int i = 0; i < count; i++) {
            if (inds[i] < 0 || inds[i] >= NDISC) return BMI_FAILURE;
            CONTEXT(self)->state.soil_discrete_storage_theta[inds[i]] = s[i];
        }
        return BMI_SUCCESS;
    }
    if (count == 1) return Set_value(self, name, src);
    return BMI_FAILURE;
}

/* ================================================================== */
/*  Grid information                                                   */
/* ================================================================== */

static int Get_grid_rank(Bmi *self, int grid, int *rank) {
    *rank = (grid == 0) ? 0 : 1;
    return BMI_SUCCESS;
}

static int Get_grid_size(Bmi *self, int grid, int *size) {
    CFE_Model_Context *ctx = CONTEXT(self);
    switch (grid) {
        case 0: *size = 1; break;
        case 1: *size = NDISC; break;
        case 2:
            if (ctx && ctx->options.surface_routing_scheme == SURF_ROUTE_NASH_CASCADE)
                *size = ctx->parameters.nash_surface_N;
            else
                *size = 0;
            break;
        case 3: *size = 2; break;
        case 4:
            if (ctx && ctx->options.surface_routing_scheme == SURF_ROUTE_GIUH)
                *size = ctx->parameters.giuh_num_ordinates;
            else
                *size = 0;
            break;
        default: return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

static int Get_grid_type(Bmi *self, int grid, char *type) {
    strcpy(type, (grid == 0) ? "scalar" : "vector");
    return BMI_SUCCESS;
}

static int Get_grid_shape(Bmi *self, int grid, int *shape) {
    if (grid > 0) return Get_grid_size(self, grid, shape);
    return BMI_FAILURE;
}

static int Get_grid_spacing  (Bmi *self, int g, double *s)       { return BMI_FAILURE; }
static int Get_grid_origin   (Bmi *self, int g, double *o)       { return BMI_FAILURE; }
static int Get_grid_x        (Bmi *self, int g, double *x)       { return BMI_FAILURE; }
static int Get_grid_y        (Bmi *self, int g, double *y)       { return BMI_FAILURE; }
static int Get_grid_z        (Bmi *self, int g, double *z)       { return BMI_FAILURE; }
static int Get_grid_node_count(Bmi *self, int g, int *c)         { return Get_grid_size(self, g, c); }
static int Get_grid_edge_count(Bmi *self, int g, int *c)         { return BMI_FAILURE; }
static int Get_grid_face_count(Bmi *self, int g, int *c)         { return BMI_FAILURE; }
static int Get_grid_edge_nodes(Bmi *self, int g, int *e)         { return BMI_FAILURE; }
static int Get_grid_face_edges(Bmi *self, int g, int *e)         { return BMI_FAILURE; }
static int Get_grid_face_nodes(Bmi *self, int g, int *n)         { return BMI_FAILURE; }
static int Get_grid_nodes_per_face(Bmi *self, int g, int *n)     { return BMI_FAILURE; }

/* ================================================================== */
/*  Registration                                                       */
/* ================================================================== */

Bmi* register_bmi_cfe(Bmi *model) {
    if (model == NULL) return NULL;

    model->data = NULL;

    model->initialize          = Initialize;
    model->update              = Update;
    model->update_until        = Update_until;
    model->finalize            = Finalize;

    model->get_component_name  = Get_component_name;
    model->get_input_item_count  = Get_input_item_count;
    model->get_output_item_count = Get_output_item_count;
    model->get_input_var_names   = Get_input_var_names;
    model->get_output_var_names  = Get_output_var_names;

    model->get_var_grid        = Get_var_grid;
    model->get_var_type        = Get_var_type;
    model->get_var_units       = Get_var_units;
    model->get_var_itemsize    = Get_var_itemsize;
    model->get_var_nbytes      = Get_var_nbytes;
    model->get_var_location    = Get_var_location;

    model->get_current_time    = Get_current_time;
    model->get_start_time      = Get_start_time;
    model->get_end_time        = Get_end_time;
    model->get_time_units      = Get_time_units;
    model->get_time_step       = Get_time_step;

    model->get_value           = Get_value;
    model->get_value_ptr       = Get_value_ptr;
    model->get_value_at_indices = Get_value_at_indices;

    model->set_value           = Set_value;
    model->set_value_at_indices = Set_value_at_indices;

    model->get_grid_rank       = Get_grid_rank;
    model->get_grid_size       = Get_grid_size;
    model->get_grid_type       = Get_grid_type;
    model->get_grid_shape      = Get_grid_shape;
    model->get_grid_spacing    = Get_grid_spacing;
    model->get_grid_origin     = Get_grid_origin;
    model->get_grid_x          = Get_grid_x;
    model->get_grid_y          = Get_grid_y;
    model->get_grid_z          = Get_grid_z;
    model->get_grid_node_count = Get_grid_node_count;
    model->get_grid_edge_count = Get_grid_edge_count;
    model->get_grid_face_count = Get_grid_face_count;
    model->get_grid_edge_nodes = Get_grid_edge_nodes;
    model->get_grid_face_edges = Get_grid_face_edges;
    model->get_grid_face_nodes = Get_grid_face_nodes;
    model->get_grid_nodes_per_face = Get_grid_nodes_per_face;

    return model;
}

CFE_Model_Context* new_bmi_cfe(void) {
    /* Allocate a zeroed context — Initialize will populate it */
    return (CFE_Model_Context*)calloc(1, sizeof(CFE_Model_Context));
}

void delete_bmi_cfe(Bmi *model) {
    if (model != NULL) {
        if (model->data != NULL) {
            Finalize(model);
        }
        free(model);
    }
}
