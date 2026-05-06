# Conceptual Functional Equivalent (CFE) Model

CFE is a simplified conceptual hydrological model designed to be functionally
equivalent to the stormflow generation components of the NOAA National Water
Model (versions 3.1 and earlier). Originally conceived by Fred L. Ogden
(NOAA/NWS Office of Water Prediction).

For the conceptual basis and hypotheses underpinning CFE, see
[MODEL.md](MODEL.md).

## Version 3 Highlights

- **Discrete Soil Moisture Balance Model (DSBM)**: optional 4-layer
  Noah-MP-style soil discretization with Darcy-Buckingham vertical fluxes
  and Clapp-Hornberger hydraulic properties.
- **Unified Nash Cascade routing** for both surface and subsurface flow,
  with optional retention depth and runon infiltration.
- **Priestley-Taylor PET** estimation from AORC radiation data (testing only).
- **New v3 config format** (`.cf3`) with comments and explicit units, plus
  full backward compatibility with v2 legacy format (`.cf2`/`.txt`).
- **BMI compliance** with the CSDMS BMI-C standard, including `get_value_ptr`
  for all variables and the ngen mass balance protocol.
- **20 calibration parameters** accessible via BMI `set_value`/`get_value_ptr`,
  with v2 parameter name aliases for backward compatibility.

## Build and Run

```bash
cmake -B build -S .
cmake --build build
ctest --test-dir build       # 46 tests: 44 unit + 2 integration
```

See [INSTALL.md](INSTALL.md) for detailed build options, driver usage, and
ngen framework integration instructions.

## Configuration

CFE supports both legacy v2 and new v3 config formats. Example configs are
in the `configs/` directory. See [configs/README.md](configs/README.md) for
parameter reference, and `configs/clean_config.cf3` for a fully annotated
v3 template.

## Testing

- **Unit tests**: 44 CTest tests covering all BMI functions, calibration
  parameter round-trips, and mass balance protocol validation.
- **Integration tests**: 2 golden output comparisons (v2 legacy + v3 DSBM)
  at 1e-10 tolerance (exact numerical match).

```bash
ctest --test-dir build -L integration   # integration only
ctest --test-dir build -E integration   # unit tests only
```

See [test/README.md](test/README.md) for test code organization.

## Migration from v2

See [CHANGELOG.md](CHANGELOG.md) for the complete v2-to-v3 variable name
mapping, calibration parameter reference, and list of removed items.

## Getting Help

For questions, please open a GitHub Issue. See [CONTRIBUTING.md](CONTRIBUTING.md)
for development guidelines.
