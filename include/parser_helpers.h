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


#ifndef PARSER_HELPERS_H
#define PARSER_HELPERS_H

#include <stddef.h>
#include "cfe_config.h"


// This structt is used to validate giuh, giuh convolution, and nash init storage arrays.
typedef struct {
    int num_giuh_ordinates_read;
    int num_giuh_convolution_read;
    int num_surf_nash_storages_read;
    int num_subsurf_nash_storages_read;
} PARSER_ARRAY_COUNTS;


// Function prototypes for parser helper functions

/**
 * Trim whitespace from the beginning and end of a string
 * @param str String to trim (modified in place)
 * @return Pointer to trimmed string
 */
char* trim_whitespace(char* str);

/**
 * Extract units from bracketed section of a line
 * @param line Input line containing bracketed units
 * @param units Output buffer for units string
 * @param units_size Size of units buffer
 */
void extract_units(const char* line, char* units, size_t units_size);

/**
 * Parse comma-separated array of double values
 * @param value_str String containing comma-separated values
 * @param array Output array to populate
 * @param max_elements Maximum number of elements to parse
 * @return Number of elements actually parsed
 */
int parse_double_array(const char* value_str, double* array, int max_elements);

/**
 * Portable case-insensitive string comparison
 * @param str1 First string to compare
 * @param str2 Second string to compare
 * @return 0 if equal, <0 if str1 < str2, >0 if str1 > str2
 */
int string_compare_ignore_case(const char* str1, const char* str2);

/**
 * Parse boolean values from string (TRUE/FALSE, 1/0)
 * @param value_str String containing boolean value
 * @return 1 for true, 0 for false
 */
int parse_boolean(const char* value_str);

/**
 * Clean quoted strings (remove quotes if present)
 * @param str String to clean (modified in place)
 */
void clean_quoted_string(char* str);

/*
 * Find the config file format version (version 2.1 or >= version 2.1)
 * @return version number
 */
double read_cfe_config_version(const char* cfg_path);

/*
 * This function parses the pre version 2.1 config file format
 */
int parse_config_legacy_format(const char* filename, CFE_CONFIG* config);

// helper function ffor parse_config_legacy_format()
int split_legacy_config_line(const char* line, char* key, char* value, char* units);

/**
 * Main function to parse CFE configuration file
 * @param filename Path to configuration file
 * @param config Pointer to CFE_CONFIG structure to populate
 * @return 0 on success, -1 on error
 */
int parse_cfe_config_ge_v2_1(const char* filename, CFE_CONFIG* config);

/**
 * The cfe2.1 parser calls this function after parsing to ensure that for each storage array read,
 *  the correct number of values was read, and properly applied.
 */
int validate_giuh_nash_arrays(CFE_CONFIG* config, const PARSER_ARRAY_COUNTS* counts);
/**
 * Print configuration values for testing/debugging
 * @param config Pointer to CFE_CONFIG structure to print
 */
void print_config(const CFE_CONFIG* config);

/**
 * Convert delimiter name to actual delimiter string
 * @param delimiter_name Name from config ("comma", "space", "tab", "pipe")
 * @return Actual delimiter string
 */
const char* get_delimiter_string(const char* delimiter_name);

/**
 * Validate time format option
 * @param time_format Format name ("timestep", "datetime", "juliandate")
 * @return 1 if valid, 0 if invalid
 */
int validate_time_format(const char* time_format);

/**
 * output path name string must end in "/".  If not, add it.
 */
void validate_and_fix_output_path(char* path_str, size_t buffer_size);

/*
 * Validate and fix output format string iff necessary
 * Returns 1 iff format was valid, 0 iff it was corrected to defaultt
 */
int validate_and_fix_output_format(char* format_str, size_t buffer_size);

/*
 * function concatenates the output filename to the output path
 */
void build_full_output_path(const char* base_path, const char* filename, 
                              char* full_path, size_t full_path_size);

#endif // PARSER_HELPERS_H
