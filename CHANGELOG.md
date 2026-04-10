All notable changes to this project will be documented in this file.
We follow the [Semantic Versioning 2.0.0](http://semver.org/) format.


## 3.0.0 - 2026-04-09

CFE v3 migration from cfe3-project reference implementation by
Fred L. Ogden (NOAA/NWS Office of Water Prediction).

### Added

#### Model
- Discrete Soil Moisture Balance Model (DSBM): 4-layer Noah-MP-style soil
  discretization (0.1, 0.3, 0.6, 1.0 m) with Darcy-Buckingham fluxes,
  Clapp-Hornberger hydraulic properties, and optional lookup tables.
  Activated via `control_soil_simulate_discrete_soil_moisture_true_false=TRUE`.
- Priestley-Taylor PET estimation from AORC radiation data (testing only).
  Activated via `control_ET_simulate_Priestley_Taylor=<alpha>`.
- New v3.0 config format (`.cf3`) with keyword=arg(s) [units] and comments.
  Full backward compatibility with legacy v2 format (`.cf2` / `.txt`).

#### BMI
- `CONTEXT(self)` macro pattern: model state (`CFE_Model_Context`) stored in
  `Bmi.data`, accessed cleanly throughout bmi_cfe.c.
- 21 output variables (was 15), 4 input variables (was 5) — see tables below.
- Array state variables on dedicated grids: `state_soil_moisture_theta` (grid 1,
  NDISC), `state_nash_surface_storage` (grid 2), `state_nash_subsurface_storage`
  (grid 3), `state_giuh_queue` (grid 4).
- Per-timestep volume balance outputs: `timestep_storage_start_m`,
  `timestep_input_m`, `timestep_output_m`, `timestep_storage_end_m`.
- `vol_balance_residual_m` output: cached volstart + volin - volout - volend.
- 20 calibration parameters via `get_value`, `set_value`, and `get_value_ptr`,
  with v2 legacy name aliases for backward compatibility with ngen calibration
  configs (e.g. `maxsmc` resolves to the same field as `soil_effective_porosity`).
- ngen mass balance protocol: `ngen::mass_in`, `ngen::mass_out`,
  `ngen::mass_stored`, `ngen::mass_leaked` via `get_value_ptr`.
- `STANDALONE` CMake option to build `cfe_main_driver` (non-BMI executable).
- 2 CTest integration tests with golden output comparison at 1e-10 tolerance
  (v2 legacy config + v3 DSBM config).

#### Core Model Improvements
- Conceptual reservoir rewritten with epsilon-based comparisons, linear
  reservoir special case (exponent ≈ 1.0), and improved flux priority logic.
- Unified `nash_cascade_routing()` replaces separate surface/subsurface
  functions, with retention depth and runon infiltration support.
- CSDMS ABI compatibility notice added to `include/bmi.h`.

### Changed

#### BMI Variable Names

Output variables (v2 → v3):

| v2 | v3 |
|----|-----|
| `Q_OUT` | `discharge_m` |
| `GIUH_RUNOFF` | `surface_runoff_m` |
| `NASH_LATERAL_RUNOFF` | `lateral_flow_m` |
| `DEEP_GW_TO_CHANNEL_FLUX` | `baseflow_m` |
| `ACTUAL_ET` | `actual_et_m` |
| `SOIL_STORAGE` | `state_soil_storage_m` |
| `GW_STORAGE` | `state_gw_storage_m` |

Input variables (v2 → v3):

| v2 | v3 |
|----|-----|
| `atmosphere_water__liquid_equivalent_precipitation_rate` | `rainfall_depth_m` (m/timestep) |
| `water_potential_evaporation_flux` | `et_potential_m` (m/timestep) |

Calibration parameters (v2 alias → v3 canonical name):

| v2 | v3 | Units |
|----|-----|-------|
| `maxsmc` | `soil_effective_porosity` | - |
| `satdk` | `soil_saturated_hydraulic_conductivity` | m s-1 |
| `slope` | `soil_percolation_rate_limiter` | - |
| `b` | `soil_Clapp_Hornberger_b` | - |
| `Klf` | `soil_lateral_flow_K` | h-1 |
| `Kn` | `subsurface_nash_K` | h-1 |
| `Cgw` | `gw_discharge_coefficient` | m s-1 |
| `expon` | `gw_discharge_exponent` | - |
| `max_gw_storage` | `gw_max_storage_m` | m |
| `satpsi` | `soil_saturated_capillary_head` | m |
| `wltsmc` | `soil_wilting_point` | - |
| `alpha_fc` | `soil_field_capacity_fraction` | - |
| `a_Xinanjiang_inflection_point_parameter` | `Xinanjiang_inflection_a` | - |
| `b_Xinanjiang_shape_parameter` | `Xinanjiang_shape_b` | - |
| `x_Xinanjiang_shape_parameter` | `Xinanjiang_shape_x` | - |
| `Kinf_nash_surface` | `surface_nash_Kinf` | h-1 |
| `retention_depth_nash_surface` | `surface_nash_retention_depth_m` | m |
| — | `Priestley_Taylor_alpha` | - |
| — | `soil_ice_imperv_threshold` | - |

#### Build System
- CMake project version bumped to 3.0.0.
- `-DNGEN=ON` defaults ON and is no longer required; `cmake -B build -S .` works.
- `Get_var_type`, `Get_var_itemsize`, `Get_var_nbytes` return `BMI_FAILURE` for
  unrecognized variable names (was silent fallthrough to double).
- Component name: `"CFE - Conceptual Functional Equivalent"` (was `"The CFE Model"`).

### Removed

- v2 standalone executables: `main.c` (BASE), `main_pass_forcings.c` (FORCING),
  `main_cfe_aorc_pet.c` (FORCINGPET), `main_cfe_aorc_pet_rz_aet.cxx` (AETROOTZONE).
- `run_cfe.sh` launcher script.
- `bmi/bmi.h` duplicate (canonical copy is `include/bmi.h`).
- BASE, FORCING, FORCINGPET, AETROOTZONE CMake options and build targets.
- `BMI_ACTIVE` compile definition (unused by v3 code).
- extern submodules: `aorc_bmi`, `evapotranspiration` (only used by removed targets;
  v3 handles PET internally via `cfe_pet_priestley_taylor.c`).
- v2 output variables: `RAIN_RATE`, `INFILTRATION_EXCESS`, `DIRECT_RUNOFF`,
  `SOIL_TO_GW_FLUX`, `POTENTIAL_ET`, `SOIL_STORAGE_CHANGE`, `SURF_RUNOFF_SCHEME`,
  `NWM_PONDED_DEPTH`.
- v2 input variables: `ice_fraction_schaake`, `ice_fraction_xinanjiang`,
  `soil_moisture_profile`.

### Fixed

- `conceptual_reservoir_flux_calc` now uses epsilon guards to avoid
  floating-point edge cases near storage thresholds.
- Legacy config parser accepts v2-only keys (`soil_params.expon`,
  `soil_params.expon_secondary`, `nsubsteps_nash_surface`) without error.
