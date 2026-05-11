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
 
 // This code contains helper functions ffor the CFE model to parse config files.
 // Note that it treats legacy (pre version 2.1) config files very differently from
 // cfe 2.1 and greater.  This is because the keywords, file structure, and variety
 // of features supported were all changed (hopefully ffor the better) in version
 // 2.1.   Keywords were also added to support the Discrete Moisture Balance Module
 // (DSBM) added with CFE 3.0.   FLO 9/2025

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser_helpers.h"
#include "cfe_config.h"
#include "cfe.h"  /* physics constants: GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2, etc. */

// Function to trim whitespace from a string
//###################
char* trim_whitespace(char* str) {
    if (!str) return NULL;

    char* end;

    while (isspace((unsigned char)*str)) str++;

    if (*str == 0) return str;

    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;

    end[1] = '\0';
    return str;
}


// Function to extract units from bracketed section
//################
void extract_units(const char* line, char* units, size_t units_size) {
    const char* start = strchr(line, '[');
    const char* end = strchr(line, ']');
    
    if (start && end && end > start) {
        start++; // Move past the '['
        size_t len = end - start;
        if (len < units_size - 1) {
            strncpy(units, start, len);
            units[len] = '\0';
            trim_whitespace(units);
        }
    } else {
        units[0] = '\0'; // No units found
    }
}

// Function to parse array values (comma-separated)
//####################
int parse_double_array(const char* value_str, double* array, int max_elements) {
    size_t len = strlen(value_str);
    char* str_copy = (char*)malloc(len + 1);
    if (!str_copy) return 0;

    memcpy(str_copy, value_str, len);
    str_copy[len] = '\0';

    char* token = strtok(str_copy, ",");
    int count = 0;

    while (token && count < max_elements) {
        char* trimmed = trim_whitespace(token);
        char* endptr;
        double value = strtod(trimmed, &endptr);
        if (endptr == trimmed || *endptr != '\0') {
            fprintf(stderr, "ERROR: Invalid number format: '%s'\n", trimmed);
            free(str_copy);
            return -1; // Return error code
        }
        array[count] = value; 
        count++;
        token = strtok(NULL, ",");
    }

    if (token != NULL) {
        fprintf(stderr, "Warning: parse_double_array truncated input after %d elements.\n", max_elements);
    }

    free(str_copy);
    return count;
}

// Portable case-insensitive string comparison
//############################
int string_compare_ignore_case(const char* str1, const char* str2) {
    while (*str1 && *str2) {
        char c1 = (*str1 >= 'A' && *str1 <= 'Z') ? *str1 + 32 : *str1;
        char c2 = (*str2 >= 'A' && *str2 <= 'Z') ? *str2 + 32 : *str2;
        if (c1 != c2) return c1 - c2;
        str1++;
        str2++;
    }
    return *str1 - *str2;
}

// Function to parse boolean values
//###############
int parse_boolean(const char* value_str) {
    char* trimmed = trim_whitespace((char*)value_str);
    if (string_compare_ignore_case(trimmed, "TRUE") == 0 || 
        string_compare_ignore_case(trimmed, "1") == 0) {
        return 1;
    }
    return 0;
}

// Function to clean quoted strings (remove quotes iff present)
//######################
void clean_quoted_string(char* str) {
    size_t len = strlen(str);
    if (len >= 2 && str[0] == '"' && str[len-1] == '"') {
        // Remove quotes
        memmove(str, str + 1, len - 2);
        str[len - 2] = '\0';
    }
}

// Helper function to make sure that N values provided when N required.
//######################
int validate_array_count(const char* keyword, int expected_count, int actual_count, 
                        double* array, int allow_single_zero) {
    
    // Special case: if only one value provided and it's 0.0, fill all with 0.0
    if (allow_single_zero && actual_count == 1 && array[0] == 0.0) {
        for (int i = 1; i < expected_count; i++) {
            array[i] = 0.0;
        }
        return expected_count; // Return corrected count
    }
    
    // Otherwise, require exact match
    if (actual_count != expected_count) {
        fprintf(stderr, "ERROR: %s requires exactly %d values, got %d\n", 
                keyword, expected_count, actual_count);
        return -1;
    }
    
    return actual_count;
}

// Returns cfe version iff CFE_CONFIG_VERSION= keyword string is found, ellse returns 0.0
//############################
double read_cfe_config_version(const char* cfg_path)
{
    FILE* f = fopen(cfg_path, "r");
    if (!f) {
        perror(cfg_path);
        return 0.0;
    }

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "cfe_config_version")) {
            double version;
            if (sscanf(line, "cfe_config_version = %lf", &version) == 1) {
                fclose(f);
                return version;
            }
        }
    }

    fclose(f);
    return 0.0;  // No version key found \u2192 assume legacy
}

// a helper ffor the legacy parsing function
//##########################
int split_legacy_config_line(const char* line, char* key, char* value, char* units)
{
    char* eq = strchr(line, '=');
    if (!eq) return 0;

    char* unit_start = strchr(eq + 1, '[');
    char* unit_end   = strchr(eq + 1, ']');

    strncpy(key, line, eq - line);
    key[eq - line] = '\0';
    trim_whitespace(key);

    if (unit_start && unit_end && unit_end > unit_start) {
        strncpy(value, eq + 1, unit_start - (eq + 1));
        value[unit_start - (eq + 1)] = '\0';
        strncpy(units, unit_start + 1, unit_end - unit_start - 1);
        units[unit_end - unit_start - 1] = '\0';
    } else {
        strncpy(value, eq + 1, 255);
        value[255] = '\0';
        units[0] = '\0';
    }

    trim_whitespace(value);
    return 1;
}
// Parser ffor the pre version 2.1 (legacy) config files
//############################
int parse_config_legacy_format(const char* filename, CFE_CONFIG* config)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open legacy config file: %s\n", filename);
        return -1;
    }
    
    // Flags to track which required inputs were found
    int yes_has_forcing_filename            = FALSE;
    int yes_has_soil_depth                  = FALSE;
    int yes_has_soil_params_bb              = FALSE;
    int yes_has_Ksat                        = FALSE;
    int yes_has_satpsi                      = FALSE;
    int yes_has_porosity                    = FALSE;
    int yes_has_wilting_point               = FALSE;  // retained ffor backward compatibility with CFE versions < 3 -FLO
    int yes_has_slop                        = FALSE;
    int yes_has_alpha_fc                    = FALSE;
    int yes_has_gw_max_storage              = FALSE;
    int yes_has_gw_init_storage             = FALSE;
    int yes_has_gw_discharge_coeff          = FALSE;
    int yes_has_gw_exponent                 = FALSE;
    int yes_has_subsurface_nash_K           = FALSE;
    int yes_has_nash_storage_subsurface     = FALSE;
    int yes_has_soil_K_to_lateral_flow      = FALSE;

// These checks were deleted in CFE 3 because the added parameters did not produce added model skill
// partly because of strong parameter interaction, but mostly because they aren't useful -FLO
//    int yes_has_surface_routing_scheme      = FALSE;
//    int yes_has_num_surface_nash_reservoirs = FALSE;
//    int yes_has_surface_routing_nash_K      = FALSE;
//    int yes_has_nash_surface_storage        = FALSE;
//    int yes_has_surface_nash_Kinf           = FALSE;
//    int yes_has_surf_retention_depth        = FALSE;

    int yes_has_partitioning_scheme_name    = FALSE;
    int yes_has_xinanjiang_a                = FALSE;
    int yes_has_xinanjiang_b                = FALSE;
    int yes_has_xinanjiang_x                = FALSE;
    int yes_has_urban_fraction              = FALSE;
    int yes_has_total_timesteps             = FALSE;
    int yes_has_verbosity                   = FALSE;

    // Initialize defaults
    
    memset(config, 0, sizeof(CFE_CONFIG));
    config->version    = 2.0;  // legacy family that includes GIUH or Nash w/ retention depth
    config->timestep_h = 1.0;  // defaultt timestep ffor CFE versions < 3.0
    strncpy(config->output_path_name, "./", sizeof(config->output_path_name)); // write output files to ./ by default

    int yes_has_surface_routing_scheme = TRUE;  // default
    snprintf(config->surface_routing_scheme_name, sizeof(config->surface_routing_scheme_name), "%s", "giuh");  // default  

    // Temporary cache: legacy "soil_storage" line gives theta (dimensionless)
    // We convert to absolute storage later after soil depth is known.
    int    has_init_soil_theta = FALSE;
    double init_soil_theta     = 0.0;

    char line[1024];
    char key[256], value[256], units[64];
    int errors = 0;

    // NOTE DEFAULTS FOR LEGACY CONFIG FILES ARE SET IN FUNCTION normalize_config_unit() in cfe_helpers.c

    //============================
    // 1) Parsing loop (no bounds or presence logic beyond setting flags)
    //============================
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        if (!split_legacy_config_line(line, key, value, units)) {
            fprintf(stderr, "WARNING: Skipping malformed config line: %s", line);
            errors++;
            continue;
        }

        if (string_compare_ignore_case(key, "forcing_file") == 0) {
            snprintf(config->control_input_forcing_filename,
                     sizeof(config->control_input_forcing_filename), "%s", value);
            yes_has_forcing_filename = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.depth") == 0 || string_compare_ignore_case(key, "soil_params.D") == 0) {
            config->soil_depth_m = atof(value);
            yes_has_soil_depth   = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.bb") == 0 || string_compare_ignore_case(key, "soil_params.b") == 0) {
            config->soil_Clapp_Hornberger_exponent_b = atof(value);
            yes_has_soil_params_bb = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.satdk") == 0) {
            // Expecting m/s in legacy, convert to cm/h (m/s * 100 cm/m * 3600 s/h)
            config->soil_sat_hydraulic_conductivity_cm_per_h = atof(value) * 100.0 * 3600.0;
            yes_has_Ksat = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.satpsi") == 0) {
            // Expecting meters in legacy, convert to cm
            config->soil_sat_capillary_head_cm = atof(value) * 100.0;
            yes_has_satpsi = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.smcmax") == 0) {
            config->soil_effective_porosity = atof(value);
            yes_has_porosity = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.wltsmc") == 0) {   // should retain backward compatibility with CFE versions < 3 -FLO
            config->soil_wilting_point_moisture_content = atof(value);
            yes_has_wilting_point = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_params.slop") == 0) {
            config->soil_to_gw_percolation_rate_limiter_0_to_1 = atof(value);
            yes_has_slop = TRUE;
        }
        else if (string_compare_ignore_case(key, "alpha_fc") == 0 || string_compare_ignore_case(key, "soil_params.alpha_fc") == 0) {
            config->soil_field_capacity_Pcap_over_Patm_0_1 = atof(value);
            yes_has_alpha_fc = TRUE;
        }
        else if (string_compare_ignore_case(key, "soil_storage") == 0) {
            // Legacy provides theta (dimensionless); convert later using soil_depth_m
            init_soil_theta     = atof(value);
            has_init_soil_theta = TRUE;
        }
        else if (string_compare_ignore_case(key, "max_gw_storage") == 0) {
            config->gw_reservoir_max_storage_m = atof(value);
            yes_has_gw_max_storage = TRUE;
        }
        else if (string_compare_ignore_case(key, "gw_storage") == 0) {
            config->gw_reservoir_init_storage_m = atof(value);
            yes_has_gw_init_storage = TRUE;
        }
        else if (string_compare_ignore_case(key, "Cgw") == 0) {
            config->gw_discharge_coeff_m_per_timestep = atof(value);
            yes_has_gw_discharge_coeff = TRUE;
        }
        else if (string_compare_ignore_case(key, "expon") == 0) {
            config->gw_discharge_exponent = atof(value);
            yes_has_gw_exponent = TRUE;
        }
        else if (string_compare_ignore_case(key, "K_nash_subsurface") == 0) {
            config->subsurface_routing_nash_K = atof(value);
            yes_has_subsurface_nash_K = TRUE;
        }
        else if (string_compare_ignore_case(key, "nash_storage_subsurface") == 0) {
            parse_double_array(value, config->subsurface_routing_nash_cascade_init_storage_m, 2);
            yes_has_nash_storage_subsurface = TRUE;
        }
        else if (string_compare_ignore_case(key, "K_lf") == 0) {
            config->soil_reservoir_rate_const_to_subsurface_lateral_flow = atof(value);
            yes_has_soil_K_to_lateral_flow = TRUE;
        }
//#####################################
//

        else if (string_compare_ignore_case(key, "surface_water_partitioning_scheme") == 0) {
            strncpy(config->partitioning_scheme_name, value,
                    sizeof(config->partitioning_scheme_name) - 1);
            config->partitioning_scheme_name[sizeof(config->partitioning_scheme_name) - 1] = '\0';
            yes_has_partitioning_scheme_name = TRUE;
        }
        else if (string_compare_ignore_case(key, "a_Xinanjiang_inflection_point_parameter") == 0) {
            config->soil_Xinanjiang_tension_water_inflection_point = atof(value);
            yes_has_xinanjiang_a = TRUE;
        }
        else if (string_compare_ignore_case(key, "b_Xinanjiang_shape_parameter") == 0) {
            config->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent = atof(value);
            yes_has_xinanjiang_b = TRUE;
        }
        else if (string_compare_ignore_case(key, "x_Xinanjiang_shape_parameter") == 0) {
            config->soil_Xinanjiang_free_water_soil_moist_distrib_exponent = atof(value);
            yes_has_xinanjiang_x = TRUE;
        }
        else if (string_compare_ignore_case(key, "urban_decimal_fraction") == 0) {
            config->cat_impervious_fraction = atof(value);
            yes_has_urban_fraction = TRUE;
        }
        else if (string_compare_ignore_case(key, "num_timesteps") == 0) {
            config->total_timesteps = atoi(value);
            yes_has_total_timesteps = TRUE;
        }
        else if (string_compare_ignore_case(key, "verbosity") == 0) {
            config->verbosity = atoi(value);
            yes_has_verbosity = TRUE;
        }
        else if (string_compare_ignore_case(key, "giuh_ordinates") == 0) {
            int count = parse_double_array(value, config->surface_routing_giuh_ordinates,
                                           MAX_NUM_GIUH_ORDINATES);
            if (count > 0) {
                config->surface_routing_num_giuh_ordinates = count;
                for (int i = 0; i < count; i++)
                    config->surface_routing_init_giuh_convolution_queue_m[i] = 0.0;
            }
        }
        // v2-only keys accepted for backward compatibility (warn but don't fail)
        else if (string_compare_ignore_case(key, "soil_params.expon") == 0 ||
                 string_compare_ignore_case(key, "soil_params.expon_secondary") == 0 ||
                 string_compare_ignore_case(key, "refkdt") == 0 ||
                 string_compare_ignore_case(key, "debug") == 0 ||
                 string_compare_ignore_case(key, "nsubsteps_nash_surface") == 0 ||
                 string_compare_ignore_case(key, "surface_runoff_scheme") == 0 ||
                 string_compare_ignore_case(key, "N_nash_surface") == 0 ||
                 string_compare_ignore_case(key, "K_nash_surface") == 0 ||
                 string_compare_ignore_case(key, "nash_storage_surface") == 0 ||
                 string_compare_ignore_case(key, "Kinf_nash_surface") == 0 ||
                 string_compare_ignore_case(key, "retention_depth_nash_surface") == 0) {
            fprintf(stderr, "WARNING: Ignoring deprecated v2 config key: %s (not supported in CFE v3)\n", key);
        }
        else {
            fprintf(stderr, "WARNING: Unknown legacy keyword: %s\n", key);
            errors++;
        }
    }

    fclose(fp);

    //============================
    // 2) Presence and option checks ONLY (no bounds checking here)
    //============================

    // Must have either forcing_file or num_timesteps
    if (!yes_has_forcing_filename && config->total_timesteps <= 0) {
        fprintf(stderr, "ERROR: Missing both forcing_file and num_timesteps. At least one is required.\n");
        return -1;
    }

    // Soil core requirements for computing initial soil storage from theta
    if (!yes_has_soil_depth) {
        fprintf(stderr, "WARNING: Missing soil_params.depth, using default 2.0 m\n");
        config->soil_depth_m = 2.0;
        yes_has_soil_depth = TRUE;
    }
    if (!yes_has_soil_params_bb) {
        fprintf(stderr, "ERROR: Missing soil_params.bb\n");
        return -1;
    }
    if (!yes_has_Ksat) {
        fprintf(stderr, "ERROR: Missing soil_params.satdk\n");
        return -1;
    }
    if (!yes_has_satpsi) {
        fprintf(stderr, "ERROR: Missing soil_params.satpsi\n");
        return -1;
    }
    if (!yes_has_porosity) {
        fprintf(stderr, "ERROR: Missing soil_params.smcmax (porosity)\n");
        return -1;
    }
    if (!yes_has_wilting_point) {    // check and make sure we've got the parameters needed to calculate it -FLO new ffor CFE 3
        if (yes_has_soil_params_bb && yes_has_satpsi && yes_has_porosity) {
        
            // new with CFE version 3 - Calculate using 15 atm capillary pressure -FLO
            double capillary_pressure_wilting_atm = 15.0;  // atmospheres (positive)
        
            // Convert to meters of water
            double g_m_per_s2 = GRAVITATIONAL_ACCELERATION_EARTH_m_per_s2;
            double rho_lw_kg_per_m3 = WATER_LIQUID_DENSITY_kg_per_m3;
            double std_atm_press_Pa = STANDARD_ATM_PRESS_Pa;
            double psi_atm_m = std_atm_press_Pa / (g_m_per_s2 * rho_lw_kg_per_m3);
            double capillary_pressure_wilting_m = capillary_pressure_wilting_atm * psi_atm_m;  // positive, meters of water
        
            // Clapp-Hornberger: theta = theta_sat * (psi_sat/psi)^(1/b)
            config->soil_wilting_point_moisture_content = config->soil_effective_porosity * 
                pow(config->soil_sat_capillary_head_cm / 100.0 / capillary_pressure_wilting_m, 
                    (1.0 / config->soil_Clapp_Hornberger_exponent_b));
            yes_has_wilting_point = TRUE;
        } else {
            fprintf(stderr, "ERROR: Unable to calculate wilting point - missing CH soil parameters\n");
        return -1;
        }
    }
    if (!yes_has_alpha_fc) {
        fprintf(stderr, "ERROR: Missing alpha_fc\n");
        return -1;
    }
    if (!yes_has_gw_max_storage) {
        fprintf(stderr, "ERROR: Missing max_gw_storage\n");
        return -1;
    }
    if (!yes_has_gw_init_storage) {
        fprintf(stderr, "ERROR: Missing gw_storage\n");
        return -1;
    }
    if (!yes_has_gw_discharge_coeff) {
        fprintf(stderr, "ERROR: Missing Cgw (GW discharge coefficient)\n");
        return -1;
    }
    if (!yes_has_gw_exponent) {
        fprintf(stderr, "ERROR: Missing exponent (GW discharge exponent)\n");
        return -1;
    }
    if (!yes_has_subsurface_nash_K) {
        fprintf(stderr, "ERROR: Missing K_nash_subsurface\n");
        return -1;
    }
    if (!yes_has_soil_K_to_lateral_flow) {
        fprintf(stderr, "ERROR: Missing K_lf\n");
        return -1;
    }
    if (!yes_has_partitioning_scheme_name) {
        fprintf(stderr, "ERROR: Missing surface_water_partitioning_scheme\n");
        return -1;
    }

    // Partitioning-specific presence checks
    if (string_compare_ignore_case(config->partitioning_scheme_name, "Xinanjiang") == 0) {
        if (!yes_has_xinanjiang_a || !yes_has_xinanjiang_b || !yes_has_xinanjiang_x) {
            fprintf(stderr, "ERROR: Missing Xinanjiang parameters (a, b, x)\n");
            return -1;
        }
    }

//
//##################################

    // Need these three ffor converting legacy soil_storage (theta) to absolute storage
    if (has_init_soil_theta) {
        if (!yes_has_soil_depth ) {
            fprintf(stderr,
                    "ERROR: soil_params.depth required to compute initial soil_storage from theta.\n");
            return -1;
        }
        config->soil_reservoir_init_storage_m = init_soil_theta * config->soil_depth_m * config->soil_effective_porosity;
        int fred_debug = FALSE;
        if(fred_debug) {
            printf("init_soil_theta =%f, soil_depth = %f m, porosity=%f init_storage=%f m\n",init_soil_theta, config->soil_depth_m,
                  config->soil_effective_porosity, config->soil_reservoir_init_storage_m);
        }
    }

   //
   // Note: parameter bounds checking is done in validate_required_parameters() that lives in cfe_helpers.c

    // Final unknown-keyword check
    if (errors > 0) {
        fprintf(stderr, "ERROR: %d unknown or malformed keyword(s) encountered.\n", errors);
        return -1;
    }

    return 0;
}


// Parsing function ffor CFE version >= 2.1
//##########################   
int parse_cfe_config_ge_v2_1(const char* filename, CFE_CONFIG* config) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return -1;
    }

    PARSER_ARRAY_COUNTS array_counts = {0};
    
    memset(config, 0, sizeof(CFE_CONFIG));
    for (int i = 0; i < MAX_NUM_SUBSURFACE_NASH_CASCADE; i++) {  // or whatever max for subsurface
        config->subsurface_routing_nash_cascade_init_storage_m[i] = 0.0;
    }
    char line[1024];
    char keyword[256];
    char value_part[768];
    char units[64];

    snprintf(config->surface_routing_scheme_name, sizeof(config->surface_routing_scheme_name), "%s", "giuh");  // default and only option
    
    while (fgets(line, sizeof(line), file)) {
        // Skip empty lines and comments
        if (line[0] == '\n' || line[0] == '#') continue;
        char* line_start = line;
        while (*line_start && isspace(*line_start)) line_start++;   // skip leading white spaces   
        char* comment_pos = strchr(line_start, '#');  // strip comments that might exist at the end of a keyword=val line.
        char* cpp_comment_pos = strstr(line_start, "//");
        if (comment_pos && cpp_comment_pos) {
            // Use whichever comes first
            if (cpp_comment_pos < comment_pos) comment_pos = cpp_comment_pos;
        } else if (cpp_comment_pos) {
            comment_pos = cpp_comment_pos;
        }
        if (comment_pos) *comment_pos = '\0';
        // Find the equals sign
        char* equals_pos = strchr(line_start, '=');
        if (!equals_pos) continue;
        // Extract keyword: start from beginning, stop at equals, trim whitespace
        char* keyword_start = line_start;
        char* keyword_end = equals_pos;
        // Skip leading whitespace in keyword
        while (keyword_start < keyword_end && isspace(*keyword_start)) {
            keyword_start++;
        }
        // Skip trailing whitespace in keyword
        while (keyword_end > keyword_start && isspace(*(keyword_end - 1))) {
            keyword_end--;
        }

        // Copy keyword
        size_t keyword_len = keyword_end - keyword_start;
        if (keyword_len >= sizeof(keyword)) keyword_len = sizeof(keyword) - 1;
        strncpy(keyword, keyword_start, keyword_len);
        keyword[keyword_len] = '\0';

        // Extract value part (everything after '=')
        char* value_start = equals_pos + 1;
        while (*value_start && isspace(*value_start)) value_start++;  

        // Handle units in brackets [unit]
        char* bracket_pos = strchr(value_start, '[');
        if (bracket_pos) {
            size_t value_len = bracket_pos - value_start;
            strncpy(value_part, value_start, value_len);
            value_part[value_len] = '\0';
        } else {
            strncpy(value_part, value_start, sizeof(value_part) - 1);
            value_part[sizeof(value_part) - 1] = '\0';
            // Remove newline if present
            char* newline = strchr(value_part, '\n');
            if (newline) *newline = '\0';
        }

        // Trim whitespace from value
        trim_whitespace(value_part);

        // Extract units
        extract_units(line, units, sizeof(units));

        // Process keywords (same as before)
        if (string_compare_ignore_case(keyword, "cfe_config_version") == 0) {
            config->version = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_model_timestep_h") == 0) {
            config->timestep_h = atof(value_part);
            snprintf(config->timestep_units, sizeof(config->timestep_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_id") == 0) {
            snprintf(config->cat_id, sizeof(config->cat_id), "%s", value_part);
        }
        else if (string_compare_ignore_case(keyword, "catchment_latitude_decimal_degree") == 0) {
            config->cat_latitude = atof(value_part);
            snprintf(config->cat_latitude_units, sizeof(config->cat_latitude_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_longitude_decimal_degree") == 0) {
            config->cat_longitude = atof(value_part);
            snprintf(config->cat_longitude_units, sizeof(config->cat_longitude_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_elevation") == 0) {
            config->cat_elev = atof(value_part);
            snprintf(config->cat_elev_units, sizeof(config->cat_elev_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_area_km2") == 0) {
            config->cat_area_km2 = atof(value_part);
            snprintf(config->cat_area_units, sizeof(config->cat_area_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "catchment_impervious_fraction_0-1") == 0) {
            config->cat_impervious_fraction = atof(value_part);
            snprintf(config->cat_impervious_units, sizeof(config->cat_impervious_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "control_input_forcing_filename") == 0) {
            snprintf(config->control_input_forcing_filename, sizeof(config->control_input_forcing_filename), "%s", value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_total_num_simulation_timesteps") == 0) {
            config->total_timesteps = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_verbosity") == 0) {
            config->verbosity = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_depth_m") == 0) {
            config->soil_depth_m = atof(value_part);
            snprintf(config->soil_depth_units, sizeof(config->soil_depth_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_Clapp_Hornberger_exponent_b") == 0) {
            config->soil_Clapp_Hornberger_exponent_b = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_sat_hydraulic_conductivity_cm_per_h") == 0) {
            config->soil_sat_hydraulic_conductivity_cm_per_h = atof(value_part);
            snprintf(config->soil_sat_hydraulic_conductivity_units, sizeof(config->soil_sat_hydraulic_conductivity_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_sat_capillary_head_cm") == 0) {
            config->soil_sat_capillary_head_cm = atof(value_part);
            snprintf(config->soil_sat_capillary_head_units, sizeof(config->soil_sat_capillary_head_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_to_gw_percolation_rate_limiter_0_to_1") == 0) {
            config->soil_to_gw_percolation_rate_limiter_0_to_1 = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_effective_porosity") == 0) {
            config->soil_effective_porosity = atof(value_part);
        }
//   Because this is not an input parameter ffor CFE versions >= 3, it is appropriately calculated from Clapp-Hornberger relation
//        else if (string_compare_ignore_case(keyword, "soil_wilting_point_moisture_content") == 0) {
//            config->soil_wilting_point_moisture_content = atof(value_part);
//        }
        else if (string_compare_ignore_case(keyword, "soil_field_capacity_Pcap_over_Patm_0_1") == 0) {
            config->soil_field_capacity_Pcap_over_Patm_0_1 = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_soil_reservoir_init_storage_m") == 0) {
            config->soil_reservoir_init_storage_m = atof(value_part);
            snprintf(config->soil_reservoir_init_storage_units, sizeof(config->soil_reservoir_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_soil_reservoir_init_discrete_storage_theta") == 0) {
            int parsed_count = parse_double_array(value_part, config->soil_reservoir_init_discrete_storage_theta, NDISCS);

            if (parsed_count == -1) {
                fprintf(stderr, "ERROR: Failed to parse soil_reservoir_init_discrete_storage_theta\n");
                return -1; // parsing error
            }

            if (parsed_count != NDISCS) {
                fprintf(stderr, "ERROR: soil_reservoir_init_discrete_storage_theta requires exactly %d values, got %d\n", 
                        NDISCS, parsed_count);
                return -1; // wrong number of elements
            }

            // Optional: validate theta values are reasonable
            for (int i = 0; i < NDISCS; i++) {
                if (config->soil_reservoir_init_discrete_storage_theta[i] < 0.0 || 
                    config->soil_reservoir_init_discrete_storage_theta[i] > 1.0) {
                    fprintf(stderr, "WARNING: soil_reservoir_init_discrete_storage_theta[%d] = %.3f is outside typical range [0.0, 1.0]\n", 
                            i, config->soil_reservoir_init_discrete_storage_theta[i]);
                }
            }

            snprintf(config->soil_reservoir_init_discrete_storage_theta_units, 
                     sizeof(config->soil_reservoir_init_discrete_storage_theta_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "soil_reservoir_rate_const_to_subsurface_lateral_flow") == 0) {
            config->soil_reservoir_rate_const_to_subsurface_lateral_flow = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "soil_ice_content_impervious_threshold") == 0) {
            config->soil_ice_content_impervious_threshold = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "gw_reservoir_max_storage_m") == 0) {
            config->gw_reservoir_max_storage_m = atof(value_part);
            snprintf(config->gw_reservoir_max_storage_units, sizeof(config->gw_reservoir_max_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_gw_reservoir_init_storage_m") == 0) {
            config->gw_reservoir_init_storage_m = atof(value_part);
            snprintf(config->gw_reservoir_init_storage_units, sizeof(config->gw_reservoir_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "gw_discharge_coeff_m_per_timestep") == 0) {
            config->gw_discharge_coeff_m_per_timestep = atof(value_part);
            snprintf(config->gw_discharge_coeff_m_per_timestep_units, sizeof(config->gw_discharge_coeff_m_per_timestep_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "gw_discharge_exponent") == 0) {
            config->gw_discharge_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "subsurface_routing_nash_reservoir_time_constant_k") == 0) {
            config->subsurface_routing_nash_K = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "state_subsurface_routing_init_nash_cascade_storage_m") == 0) {
            array_counts.num_subsurf_nash_storages_read = parse_double_array(value_part, 
                                                                             config->subsurface_routing_nash_cascade_init_storage_m, 
                                                                             MAX_NUM_SUBSURFACE_NASH_CASCADE);
            snprintf(config->subsurface_routing_nash_cascade_init_storage_units, sizeof(config->subsurface_routing_nash_cascade_init_storage_units), "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_simulate_freeze_thaw_true_false") == 0) {
            config->control_soil_simulate_freeze_thaw_true_false = parse_boolean(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_simulate_discrete_soil_moisture_true_false") == 0) {
            config->control_soil_simulate_discrete_soil_moisture_true_false = parse_boolean(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_soil_use_lookup_table_num_points") == 0) {
            config->control_soil_use_lookup_table_num_points = atoi(value_part);
            if(config->verbosity > 1) printf("DEBUG: Parsed control_soil_use_lookup_table_num_points = %s -> %d\n", 
                                             value_part, config->control_soil_use_lookup_table_num_points);
        }
        else if (string_compare_ignore_case(keyword, "control_ET_simulate_Priestley_Taylor") == 0) {
            config->et_alpha_pt = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "control_ET_deepest_root_zone_discretization") == 0) {
            config->control_ET_deepest_root_zone_discretization = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_scheme_name") == 0) {
            // Convert to lowercase for consistent comparison
            char temp_value[sizeof(config->partitioning_scheme_name)];
            snprintf(temp_value, sizeof(temp_value), "%s", value_part);
            for (int i = 0; temp_value[i]; i++) {
                temp_value[i] = tolower(temp_value[i]);
            }
            snprintf(config->partitioning_scheme_name, sizeof(config->partitioning_scheme_name), "%s", temp_value);
        }
        else if (string_compare_ignore_case(keyword, "surface_routing_num_giuh_ordinates") == 0) {
            config->surface_routing_num_giuh_ordinates = atoi(value_part);
        }
        else if (string_compare_ignore_case(keyword, "surface_routing_giuh_ordinates") == 0) {
            array_counts.num_giuh_ordinates_read = parse_double_array(value_part,
                                                                      config->surface_routing_giuh_ordinates,
                                                                      MAX_NUM_GIUH_ORDINATES);
            snprintf(config->surface_routing_giuh_ordinates_units,
                     sizeof(config->surface_routing_giuh_ordinates_units),
                     "%s", units);
        }
        else if (string_compare_ignore_case(keyword, "state_surface_routing_init_giuh_convolution_queue_m") == 0) {
            array_counts.num_giuh_convolution_read = parse_double_array(value_part,
                               config->surface_routing_init_giuh_convolution_queue_m,
                               MAX_NUM_GIUH_ORDINATES);
            snprintf(config->surface_routing_init_giuh_convolution_queue_units,
                     sizeof(config->surface_routing_init_giuh_convolution_queue_units),
                     "%s", units);
        }


        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_tension_water_inflection_point") == 0) {
            config->soil_Xinanjiang_tension_water_inflection_point = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent") == 0) {
            config->soil_Xinanjiang_tension_water_soil_moist_distrib_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent") == 0) {
            config->soil_Xinanjiang_free_water_soil_moist_distrib_exponent = atof(value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_status_warnings_filename") == 0) {
            snprintf(config->output_status_warnings_filename, sizeof(config->output_status_warnings_filename), "%s", value_part);
            clean_quoted_string(config->output_status_warnings_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_internal_fluxes_m_per_timestep_filename") == 0) {
            snprintf(config->output_internal_fluxes_filename, sizeof(config->output_internal_fluxes_filename), "%s", value_part);
            clean_quoted_string(config->output_internal_fluxes_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_internal_storages_m_per_timestep_filename") == 0) {
            snprintf(config->output_internal_storages_filename, sizeof(config->output_internal_storages_filename), "%s", value_part);
            clean_quoted_string(config->output_internal_storages_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_volume_balance_filename") == 0) {
            snprintf(config->output_volume_balance_filename, sizeof(config->output_volume_balance_filename), "%s", value_part);
            clean_quoted_string(config->output_volume_balance_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_soil_moisture_theta_filename") == 0) {
            snprintf(config->output_soil_moisture_theta_filename, sizeof(config->output_soil_moisture_theta_filename), "%s", value_part);
            clean_quoted_string(config->output_soil_moisture_theta_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_discharge_m_per_timestep_filename") == 0) {
            snprintf(config->output_discharge_filename, sizeof(config->output_discharge_filename), "%s", value_part);
            clean_quoted_string(config->output_discharge_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_total_discharge_m3_per_sec_filename") == 0) {
            snprintf(config->output_total_discharge_m3_per_sec_filename, 
                     sizeof(config->output_total_discharge_m3_per_sec_filename), "%s", value_part);
            clean_quoted_string(config->output_total_discharge_m3_per_sec_filename);
        }
        else if (string_compare_ignore_case(keyword, "output_path_name") == 0) {
            snprintf(config->output_path_name, sizeof(config->output_path_name), "%s", value_part);
            clean_quoted_string(config->output_path_name);
        }
        else if (string_compare_ignore_case(keyword, "output_time_standard_format") == 0) {
            snprintf(config->output_time_standard_format, sizeof(config->output_time_standard_format), "%s", value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_file_delimiter") == 0) {
            snprintf(config->output_file_delimiter, sizeof(config->output_file_delimiter), "%s", value_part);
        }
        else if (string_compare_ignore_case(keyword, "output_value_format") == 0) {
            snprintf(config->output_value_format, sizeof(config->output_value_format), "%s", value_part);
            clean_quoted_string(config->output_value_format);

            // ensure that the specified output value format is valid- "%.Ne" or "%.Nf" where 0<N<16
            // iff not, set to defaultt "%.8f"
            if (!validate_and_fix_output_format(config->output_value_format, sizeof(config->output_value_format))) {
                if (config->verbosity > 0) {
                    fprintf(stderr, "WARNING: Invalid output_value_format, using default %.8f\n");
                }
            }
        }
        else if (string_compare_ignore_case(keyword, "output_new_config_filename_prefix") == 0) {
            snprintf(config->output_new_config_filename, sizeof(config->output_new_config_filename), "%s", value_part);
            clean_quoted_string(config->output_new_config_filename);
            
        }

    }

    fclose(file);
    return 0;
}  //  <--------------------- END OF VERSION > 2.0 parser

// This function concatenates the output filename to the output path
void build_full_output_path(const char* base_path, const char* filename, 
                              char* full_path, size_t full_path_size) {
      if (base_path && strlen(base_path) > 0) {
          snprintf(full_path, full_path_size, "%s%s", base_path, filename);
      } else {
          snprintf(full_path, full_path_size, "%s", filename);
      }
  }
  
// Note: output_path_name string must end in "/".  If not, add the "/".
//###############################
void validate_and_fix_output_path(char* path_str, size_t buffer_size) {
      if (!path_str || strlen(path_str) == 0) {
          snprintf(path_str, buffer_size, "./");
          return;
      }
      
      size_t len = strlen(path_str);
      if (path_str[len-1] != '/') {
          if (len < buffer_size - 1) {
              strcat(path_str, "/");
          }
      }
  }
  
// Validate and fix output format string iff necessary
// Returns 1 iff format was valid, 0 iff it was corrected to defaultt
//################################
int validate_and_fix_output_format(char* format_str, size_t buffer_size) {
    if (!format_str || strlen(format_str) == 0) {
        snprintf(format_str, buffer_size, "%%.8f");
        return 0;
    }
    
    // Check if format matches %.Nf or %.Ne pattern
    int precision;
    char format_type;
    
    // Try to parse the format string
    int matches = sscanf(format_str, "%%.%d%c", &precision, &format_type);
    
    if (matches == 2 && 
        (format_type == 'f' || format_type == 'e') && 
        precision >= 1 && precision <= 15) {
        // Valid format
        return 1;
    }
    
    // Invalid format - set to default
    snprintf(format_str, buffer_size, "%%.8f");
    return 0;
}

// For each of the initial storages read ensure that the correct number of values was read and applied
//###########################
int validate_giuh_arrays(CFE_CONFIG* config, const PARSER_ARRAY_COUNTS* counts) {

    
    // GIUH validation----------------------
    int expected_num = config->surface_routing_num_giuh_ordinates;
    int ordinates_read = counts->num_giuh_ordinates_read;
    int queue_read = counts->num_giuh_convolution_read;
    
    if (ordinates_read != expected_num) {
        fprintf(stderr, "ERROR: Expected %d GIUH ordinates, read: %d\n", 
                expected_num, ordinates_read);
        return -1;
    }
    
    if (queue_read != expected_num) {
        fprintf(stderr, "ERROR: Expected %d GIUH queue values, read: %d\n", 
                expected_num, queue_read);
        return -1;
    }
    
    // Validate GIUH ordinates sum to 1.0 and are non-negative
    double sum = 0.0;
    int is_negative = FALSE;
    for (int i = 0; i < expected_num; i++) {
        sum += config->surface_routing_giuh_ordinates[i];
        if(config->surface_routing_giuh_ordinates[i] < 0.0) is_negative = TRUE;
    }
    if (fabs(sum - 1.0) > 1e-6) {
        fprintf(stderr, "ERROR: GIUH ordinates sum to %.6f, must sum to 1.0\n", sum);
        return -1;
    }
    if(is_negative) {
        fprintf(stderr,"ERROR: negative GIUH ordinate encountered. All ordinates must be positive and sum to 1.0.\n");
        return -1;
    }
    
    // Check for negative GIUH queue values
    for (int i = 0; i < expected_num; i++) {
        if (config->surface_routing_init_giuh_convolution_queue_m[i] < 0.0) {
            fprintf(stderr, "ERROR: GIUH queue[%d] cannot be negative: %.6f\n", 
                    i, config->surface_routing_init_giuh_convolution_queue_m[i]);
            return -1;
        }
    }
    
    
    // Subsurface Nash validation (always expects 2)-----------------
    int subsurf_read = counts->num_subsurf_nash_storages_read;
    if (subsurf_read != 2) {
        fprintf(stderr, "ERROR: Expected 2 subsurface Nash storage values, read: %d\n", 
                subsurf_read);
        return -1;
    }
    
    // Check for negative subsurface Nash storage values
    for (int i = 0; i < 2; i++) {
        if (config->subsurface_routing_nash_cascade_init_storage_m[i] < 0.0) {
            fprintf(stderr, "ERROR: Subsurface Nash storage[%d] cannot be negative: %.6f\n", 
                    i, config->subsurface_routing_nash_cascade_init_storage_m[i]);
            return -1;
        }
    }
    
    return 0;
}


// Function to print the configuration (ffor testing) - Updated to include new fields
//###############
void print_config(const CFE_CONFIG* config) {
    printf("CFE Configuration:\n");
    printf("Version: %f\n", config->version);
    printf("Timestep: %.1f %s\n", config->timestep_h, config->timestep_units);
    printf("Catchment ID: %s\n", config->cat_id);
    printf("Latitude: %.6f %s\n", config->cat_latitude, config->cat_latitude_units);
    printf("Longitude: %.6f %s\n", config->cat_longitude, config->cat_longitude_units);
    printf("Elevation: %.2f %s\n", config->cat_elev, config->cat_elev_units);
    printf("Area: %.2f %s\n", config->cat_area_km2, config->cat_area_units);
    printf("Impervious fraction: %.3f %s\n", config->cat_impervious_fraction, config->cat_impervious_units);
    printf("Forcing file: %s\n", config->control_input_forcing_filename);
    printf("Total timesteps: %d\n", config->total_timesteps);
    printf("Verbosity: %d\n", config->verbosity);
    printf("Soil depth: %.2f %s\n", config->soil_depth_m, config->soil_depth_units);
    printf("Soil Clapp-Hornberger exponent: %.3f\n", config->soil_Clapp_Hornberger_exponent_b);
    printf("Soil sat hydraulic conductivity: %.4f %s\n", 
           config->soil_sat_hydraulic_conductivity_cm_per_h, 
           config->soil_sat_hydraulic_conductivity_units);
    printf("Soil discrete storage theta: %.2f, %.2f, %.2f, %.2f %s\n",
           config->soil_reservoir_init_discrete_storage_theta[0],
           config->soil_reservoir_init_discrete_storage_theta[1],
           config->soil_reservoir_init_discrete_storage_theta[2],
           config->soil_reservoir_init_discrete_storage_theta[3],
           config->soil_reservoir_init_discrete_storage_theta_units);
    printf("GW max storage: %.3f %s\n", config->gw_reservoir_max_storage_m, config->gw_reservoir_max_storage_units);
    printf("Rainfall partitioning scheme: %s\n", config->partitioning_scheme_name);
    printf("Surface routing scheme: %s\n", config->surface_routing_scheme_name);
    printf("Control discrete soil moisture: %s\n", 
           config->control_soil_simulate_discrete_soil_moisture_true_false ? "TRUE" : "FALSE");
    printf("Control use soil lookup table: %s\n", 
           config->control_soil_use_lookup_table_num_points ? "TRUE" : "FALSE");
    printf("Number of look-up table points to calculate: %d\n",config->control_soil_use_lookup_table_num_points);
    printf("\nOutput Configuration:\n");
    printf("Status warnings: %s\n", config->output_status_warnings_filename);
    printf("Internal fluxes: %s\n", config->output_internal_fluxes_filename);
    printf("Internal storages: %s\n", config->output_internal_storages_filename);
    printf("Volume balance: %s\n", config->output_volume_balance_filename);
    printf("Soil moisture theta: %s\n", config->output_soil_moisture_theta_filename);
    printf("Discharge: %s\n", config->output_discharge_filename);
    printf("Time format: %s\n", config->output_time_standard_format);
    printf("File delimiter: %s\n", config->output_file_delimiter);
}

// Add these functions to parser_helpers.c:

// Convert delimiter name to actual delimiter character
//##############################
const char* get_delimiter_string(const char* delimiter_name) {
    if (!delimiter_name) return ",";  // Default
    
    if (string_compare_ignore_case(delimiter_name, "comma") == 0) {
        return ",";
    } else if (string_compare_ignore_case(delimiter_name, "space") == 0) {
        return " ";
    } else if (string_compare_ignore_case(delimiter_name, "tab") == 0) {
        return "\t";
    } else if (string_compare_ignore_case(delimiter_name, "pipe") == 0) {
        return "|";
    } else {
        return ",";  // Default fallback
    }
}


// Validate time format
//######################
int validate_time_format(const char* time_format) {
    if (!time_format) return 0;
    
    if (string_compare_ignore_case(time_format, "timestep") == 0 ||
        string_compare_ignore_case(time_format, "datetime") == 0 ||
        string_compare_ignore_case(time_format, "juliandate") == 0) {
        return 1;  // Valid
    }
    
    
    return 0;  // Invalid
}
