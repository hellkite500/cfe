# CFE v3 Tests

Tests are configured for use with CTest directly within CMakeLists.txt.

## Test Categories

### BMI Unit Tests (44 tests)

- **combined_bmi_funcs_test**: Legacy batch test exercising all BMI functions
  in a single routine, including a 20-parameter calibration round-trip
  (`set_value` → `get_value` → `get_value_ptr` for each parameter).
- **test_bmi_model**: Individual unit tests for each BMI function, invoked
  by name (e.g. `test_bmi_model test_initialize`).
- **test_mass_balance_protocol**: Validates the ngen mass balance protocol
  identity (`mass_in = mass_out + mass_stored + mass_leaked`) by running
  10 timesteps with known forcing and checking closure to machine epsilon.

### Integration Tests (2 tests)

- **integration_v2_legacy**: Runs `cfe_bmi_driver` with the legacy v2 config
  (`legacy_cfe_config_cat87.cf2`) and compares discharge, fluxes, and storage
  against golden reference outputs at 1e-10 tolerance.
- **integration_v3_dsbm**: Same, with the v3 DSBM config
  (`bmi_config_cat87_v3.cf3`), including soil moisture theta comparison.

Golden reference outputs are in `test/golden/v2/` and `test/golden/v3/`.

## Running Tests

```bash
# Build
cmake -B build -S .
cmake --build build

# All tests
ctest --test-dir build

# Unit tests only
ctest --test-dir build -E integration

# Integration tests only
ctest --test-dir build -L integration

# Verbose output on failure
ctest --test-dir build --output-on-failure

# Single test by name
ctest --test-dir build -R test_mass_balance_protocol
```

## Test Code Organization

### combined_bmi_funcs_test.c
Legacy batch test. Exercises all BMI functions in sequence, including
calibration parameter set/get/ptr round-trip for all 20 parameters.

### test_bmi_model.c
Individual BMI function unit tests. The `main` function dispatches to
specific test functions based on command-line argument. Also contains
`test_mass_balance_protocol`.

### bmi_test_utils.[c,h]
Test fixture struct, setup/teardown, BMI variable name arrays, and helper
functions for getting/setting values. Expected variable names and counts
are defined here.

### general_test_utils.[c,h]
Generic comparison utilities (doubles, ints, strings).

### compare_outputs.sh
Numerical comparison script for golden output testing. Compares two
output files line by line, skipping comments, with configurable relative
error tolerance.

### run_integration_test.sh
Wrapper script that runs `cfe_bmi_driver` and compares all output files
against a golden directory.

## Adding New Tests

1. Add test code (new function in `test_bmi_model.c` or new source file).
2. If new source file: add `add_executable` in CMakeLists.txt.
3. Add `add_test(NAME ... COMMAND ...)` in CMakeLists.txt.
4. Rebuild: `cmake -B build -S . && cmake --build build`
5. Run: `ctest --test-dir build -R <test_name>`
