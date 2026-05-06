# CFE v3 — Build, Test, and Run Instructions

Building CFE requires a C compiler (GCC or Clang) and [CMake](https://cmake.org/) >= 3.10.

## Quick Start

```bash
git clone https://github.com/NOAA-OWP/cfe
cd cfe
cmake -B build -S .
cmake --build build
ctest --test-dir build
```

This builds the shared library (`libcfebmi`), the BMI driver executable
(`cfe_bmi_driver`), and all test executables, then runs the full test suite
(unit tests + integration tests with golden output comparison).

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `STANDALONE` | `OFF` | Build the standalone non-BMI driver (`cfe_main_driver`) |
| `NGEN` | `ON` | Accepted for backward compatibility (no effect — always builds) |
| `CMAKE_BUILD_TYPE` | (none) | Set to `Debug` for debug symbols and verbosity=1 |

Example with standalone driver:
```bash
cmake -B build -S . -DSTANDALONE=ON
cmake --build build
```

## Running Tests

### All tests
```bash
ctest --test-dir build
```

### Unit tests only
```bash
ctest --test-dir build -E integration
```

### Integration tests only
```bash
ctest --test-dir build -L integration
```

The integration tests run the `cfe_bmi_driver` with v2 legacy and v3 DSBM
configs against golden reference outputs at 1e-10 tolerance (exact match).

### Verbose output on failure
```bash
ctest --test-dir build --output-on-failure
```

## Running the BMI Driver

The `cfe_bmi_driver` executable reads a config file and forcing data, runs
the model through the BMI interface, and writes output files.

### With v2 legacy config
```bash
build/cfe_bmi_driver \
    -c configs/legacy_cfe_config_cat87.cf2 \
    -f forcings/cat87_01Dec2015.csv \
    -q output/q.out \
    -x output/fluxes.out \
    -s output/storage.out \
    -v 1
```

### With v3 DSBM config
```bash
build/cfe_bmi_driver \
    -c configs/bmi_config_cat87_v3.cf3 \
    -f forcings/cat87_01Dec2015.csv \
    -q output/q.out \
    -x output/fluxes.out \
    -s output/storage.out \
    -t output/thetas.out \
    -v 1
```

### Driver options
```
Usage:
  cfe_bmi_driver -c <config> -f <forcing> [OPTIONS]

Required:
  -c <file>     Configuration file (.cf2 or .cf3)

Optional:
  -f <file>     Forcing data file (overrides config)
  -q <file>     Discharge output file (m/timestep)
  -b <file>     Volume balance summary file
  -x <file>     Internal fluxes output file
  -s <file>     Internal storages output file
  -t <file>     Soil moisture theta output (DSBM only)
  -v <level>    Verbosity (0=quiet, 1=normal, 2=verbose)
  -dryrun       Run 120 timesteps without forcing data
```

## Config File Formats

CFE v3 supports two config formats:

### Legacy v2 format (`.cf2` / `.txt`)
Key=value pairs, one per line. Units in brackets. Example:
```
forcing_file=BMI
soil_params.depth=2.0[m]
soil_params.b=4.05[]
soil_params.satdk=0.00000338[m s-1]
...
```

### v3 keyword format (`.cf3`)
Comments with `#` or `//`, richer keyword names, explicit units. Example:
```
cfe_config_version=3.0[]
soil_depth_m=2.0[m]
soil_Clapp_Hornberger_exponent_b=4.05[]
soil_sat_hydraulic_conductivity_cm_per_h=1.2168[cm h-1]
control_soil_simulate_discrete_soil_moisture_true_false=TRUE
...
```

See `configs/clean_config.cf3` for a fully annotated template of all v3 options.

## ngen Framework Integration

CFE v3 builds a shared library (`libcfebmi.dylib` / `libcfebmi.so`) that
can be loaded by the [ngen](https://github.com/NOAA-OWP/ngen) framework.

### Building within ngen
```bash
# From ngen directory
git submodule update --remote extern/cfe/cfe
cmake -B extern/cfe/cfe/cmake_build -S extern/cfe/cfe/
make -C extern/cfe/cfe/cmake_build
```

The realization config should reference the shared library and a CFE config
file (either v2 or v3 format):

```json
{
    "model_type_name": "bmi_c",
    "library_file": "./extern/cfe/cfe/cmake_build/libcfebmi.so",
    "init_config": "./extern/cfe/cfe/configs/legacy_cfe_config_cat87.cf2"
}
```

### BMI Variable Names

See `CHANGELOG.md` for the complete v2→v3 variable name mapping.

**Inputs:** `rainfall_depth_m`, `et_potential_m`, `verbosity`, `forcing_file_path`

**Key outputs:** `discharge_m`, `surface_runoff_m`, `lateral_flow_m`, `baseflow_m`, `actual_et_m`

**Calibration parameters** (20 total, accessible via `set_value` / `get_value_ptr`):
`soil_effective_porosity`, `soil_saturated_hydraulic_conductivity`,
`soil_Clapp_Hornberger_b`, `gw_discharge_coefficient`, `gw_discharge_exponent`,
`gw_max_storage_m`, etc. See `bmi_cfe.c` `param_var_names[]` for the full list.
v2 parameter names (e.g. `maxsmc`, `satdk`, `Cgw`) are accepted as aliases.

### Notes

- `include/bmi.h` must be the canonical [CSDMS BMI-C](https://github.com/csdms/bmi-c)
  release header. ABI compatibility with ngen cannot be guaranteed if this file
  is modified or replaced with a non-standard version.
- The `-DNGEN=ON` flag is accepted for backward compatibility but is no longer
  required. A plain `cmake -B build -S .` builds everything.
