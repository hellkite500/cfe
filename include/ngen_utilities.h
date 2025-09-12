/*
Author: Nels Frazier
Copyright (C) 2025 Lynker
------------------------------------------------------------------------
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
------------------------------------------------------------------------

Version 0.1
Protocol utilities for use with the ngen model engine
*/
#define NGEN_MASS_IN "ngen::mass_in"
#define NGEN_MASS_OUT "ngen::mass_out"
#define NGEN_MASS_STORED "ngen::mass_stored"
#define NGEN_MASS_LEAKED "ngen::mass_leaked"

#define MASS_BALANCE_VAR_NAME_COUNT 4

#include <bmi.h>
#include <string.h>

//These are the bare minimum required to interface with ngen
static const char *mass_balance_var_names[MASS_BALANCE_VAR_NAME_COUNT] = { NGEN_MASS_IN, NGEN_MASS_OUT, NGEN_MASS_STORED, NGEN_MASS_LEAKED };
static const char *mass_balance_var_types[MASS_BALANCE_VAR_NAME_COUNT] = { "double", "double", "double", "double" };
static const int   mass_balance_var_item_count[MASS_BALANCE_VAR_NAME_COUNT] = { 1, 1, 1, 1};
static const char *mass_balance_var_units[MASS_BALANCE_VAR_NAME_COUNT] = { "m", "m", "m", "m" };
//These may be useful in the future, but aren't strictly necessary for the protocol at the moment
// static const char *mass_balance_var_grids[MASS_BALANCE_VAR_NAME_COUNT] = { 0, 0, 0 };
// static const char *mass_balance_var_locations[MASS_BALANCE_VAR_NAME_COUNT] = { "node", "node", "node" };

static inline int get_mass_balance_var_type(const char *name, char *type)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            strncpy(type, mass_balance_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}

static inline int get_mass_balance_item_count(const char *name)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            return mass_balance_var_item_count[i];
        }
    }
    return -1; // Default to -1 if not found
}

static inline int get_mass_balance_unit(const char *name, char* unit)
{
    for (int i = 0; i < MASS_BALANCE_VAR_NAME_COUNT; ++i) {
        if (strcmp(name, mass_balance_var_names[i]) == 0) {
            strncpy(unit, mass_balance_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    return BMI_FAILURE;
}