# CFE v3 Configuration Files

CFE v3 uses a keyword-based config format (`.cf3`). Legacy v2 configs are not
supported directly — use `cfe_migrate_config` to convert them. For the v2
codebase, see the [v2.1.0 tag](https://github.com/NOAA-OWP/cfe/tree/v2.1.0).

## Config Format

Comments with `#` or `//`, one keyword per line, optional units in brackets:
```
cfe_config_version=3.0[]
control_model_timestep_h=1.0[h]
soil_depth_m=2.0[m]
soil_sat_hydraulic_conductivity_cm_per_h=1.2168[cm h-1]
control_soil_simulate_discrete_soil_moisture_true_false=TRUE
```

See `clean_config.cf3` for a fully annotated template of all options.

## Example Configs

| File | Partitioning | Soil | Notes |
|------|-------------|------|-------|
| `cfe_config_cat_87_pass.cf3` | Schaake | Linear reservoir | Used by unit tests |
| `migrated_legacy_cat87.cf3` | Xinanjiang | Linear reservoir | Integration test (migrated from v2) |
| `bmi_config_cat87_v3.cf3` | Schaake | DSBM (4-layer) | Integration test, DSBM enabled |
| `cfe_config_cat87_v3.cf3` | Schaake | DSBM (4-layer) | Standalone driver config |
| `clean_config.cf3` | — | — | Annotated template |

## Migrating from v2

```bash
build/cfe_migrate_config old_config.txt new_config.cf3
```

The migration utility converts key names, units (satdk m/s to cm/h, satpsi m
to cm), and storage representations (theta to meters). Configs that relied on
Nash Cascade surface routing will error — GIUH ordinates are required in v3.

## Parameter Reference

Parameters marked with * are calibratable via BMI `set_value` /
`get_value_ptr`.

### Soil Parameters

| Config Key | BMI Name | Units | Calibration Range | Description |
|-----------|----------|-------|-------------------|-------------|
| `soil_depth_m` | `param_soil_depth_m` | m | — | Soil column depth |
| `soil_Clapp_Hornberger_exponent_b`* | `soil_Clapp_Hornberger_b` | - | 2.0–15.0 | Clapp-Hornberger exponent |
| `soil_sat_hydraulic_conductivity_cm_per_h`* | `soil_saturated_hydraulic_conductivity` | m s-1 (internal) | 0.07–510.0 cm/h | Saturated hydraulic conductivity |
| `soil_sat_capillary_head_cm`* | `soil_saturated_capillary_head` | m (internal) | — | Saturated capillary head |
| `soil_effective_porosity`* | `soil_effective_porosity` | - | 0.16–0.58 (optional) | Effective porosity |
| `soil_wilting_point_moisture_content` | — | - | — | Auto-calculated from Clapp-Hornberger at 15 atm |
| `soil_to_gw_percolation_rate_limiter_0_to_1`* | `soil_percolation_rate_limiter` | - | 0.0–1.0 | Percolation rate limiter |
| `soil_field_capacity_Pcap_over_Patm_0_1`* | `soil_field_capacity_fraction` | - | 0.15 (sand) / 0.33 (other) | Field capacity; typically not calibrated |
| `soil_reservoir_rate_const_to_subsurface_lateral_flow`* | `soil_lateral_flow_K` | h-1 | 0.0–1.0 | Lateral flow rate constant |
| `state_soil_reservoir_init_storage_m` | — | m | — | Initial soil storage |

Note: `refkdt` is a constant (=3.0, per Schaake et al. 1996) hardcoded in the
model. It is not a calibration parameter — calibrating both `refkdt` and `satdk`
is degenerate (a change in one can be fully compensated by the other).

### Groundwater Parameters

| Config Key | BMI Name | Units | Calibration Range | Description |
|-----------|----------|-------|-------------------|-------------|
| `gw_reservoir_max_storage_m`* | `gw_max_storage_m` | m | 0.1–0.25 | Maximum GW storage |
| `gw_discharge_coeff_m_per_timestep`* | `gw_discharge_coefficient` | m s-1 (internal) | 1.8e-6–1.8e-3 | GW discharge coefficient |
| `gw_discharge_exponent`* | `gw_discharge_exponent` | - | 1.0–8.0 | GW discharge exponent |
| `state_gw_reservoir_init_storage_m` | — | m | — | Initial GW storage |

### GIUH Surface Routing

| Config Key | Units | Description |
|-----------|-------|-------------|
| `surface_routing_num_giuh_ordinates` | — | Number of GIUH ordinates (required) |
| `surface_routing_giuh_ordinates` | - | Comma-separated ordinates summing to 1.0 |
| `state_surface_routing_init_giuh_convolution_queue_m` | m | Initial convolution queue (one per ordinate) |

GIUH is the only surface routing method in CFE v3. Nash Cascade surface routing
was removed as it did not produce added model skill.

### Subsurface Routing (Nash Cascade)

| Config Key | BMI Name | Units | Calibration Range | Description |
|-----------|----------|-------|-------------------|-------------|
| `subsurface_routing_nash_reservoir_time_constant_k`* | `subsurface_nash_K` | h-1 | 0.0–1.0 | Subsurface Nash time constant |
| `state_subsurface_routing_init_nash_cascade_storage_m` | — | m | — | Initial Nash cascade storage (2 values) |

### Xinanjiang Parameters (when `partitioning_scheme_name=XINANJIANG`)

| Config Key | BMI Name | Units | Calibration Range |
|-----------|----------|-------|-------------------|
| `partitioning_Xinanjiang_tension_water_inflection_point`* | `Xinanjiang_inflection_a` | - | -0.493–0.493 |
| `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent`* | `Xinanjiang_shape_b` | - | 0.0–1.0 |
| `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent`* | `Xinanjiang_shape_x` | - | 0.0–1.0 |

### Additional Calibration Parameters

| Config Key | BMI Name | Units | Calibration Range | Description |
|-----------|----------|-------|-------------------|-------------|
| `control_ET_simulate_Priestley_Taylor`* | `Priestley_Taylor_alpha` | - | 0.75–1.6 | P-T alpha coefficient (0 disables PT PET) |
| `soil_ice_content_impervious_threshold`* | `soil_ice_imperv_threshold` | - | — | Ice fraction threshold for Schaake |

### DSBM and v3 Feature Controls

| Config Key | Description |
|-----------|-------------|
| `control_soil_simulate_discrete_soil_moisture_true_false` | `TRUE` to enable DSBM |
| `control_soil_use_lookup_table_num_points` | LUT points (0 = analytic) |
| `control_ET_deepest_root_zone_discretization` | Deepest root zone layer (1–NDISC) |

Calibration ranges are from NWM calibration experience (OWP/FLO). Parameters
marked `—` have no established calibration range or are not typically calibrated.

## Infiltration Excess Partitioning

1. **Schaake** (`partitioning_scheme_name=SCHAAKE`) — no additional parameters
2. **Xinanjiang** (`partitioning_scheme_name=XINANJIANG`) — requires the three Xinanjiang parameters above
