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

Parameters marked with <sup>\*</sup> are calibratable via BMI `set_value` /
`get_value_ptr`.

### Soil Parameters

| Config Key | BMI Name | Units | Description |
|-----------|----------|-------|-------------|
| `soil_depth_m` | `param_soil_depth_m` | m | Soil column depth |
| `soil_Clapp_Hornberger_exponent_b`<sup>\*</sup> | `soil_Clapp_Hornberger_b` | - | Clapp-Hornberger exponent |
| `soil_sat_hydraulic_conductivity_cm_per_h`<sup>\*</sup> | `soil_saturated_hydraulic_conductivity` | m s-1 (internal) | Saturated hydraulic conductivity |
| `soil_sat_capillary_head_cm`<sup>\*</sup> | `soil_saturated_capillary_head` | m (internal) | Saturated capillary head |
| `soil_effective_porosity`<sup>\*</sup> | `soil_effective_porosity` | - | Effective porosity |
| `soil_wilting_point_moisture_content`<sup>\*</sup> | `soil_wilting_point` | - | Wilting point |
| `soil_to_gw_percolation_rate_limiter_0_to_1`<sup>\*</sup> | `soil_percolation_rate_limiter` | - | Percolation rate limiter (0-1) |
| `soil_field_capacity_Pcap_over_Patm_0_1`<sup>\*</sup> | `soil_field_capacity_fraction` | - | Field capacity (Pcap/Patm) |
| `soil_reservoir_rate_const_to_subsurface_lateral_flow`<sup>\*</sup> | `soil_lateral_flow_K` | h-1 | Lateral flow rate constant |
| `state_soil_reservoir_init_storage_m` | — | m | Initial soil storage |

### Groundwater Parameters

| Config Key | BMI Name | Units | Description |
|-----------|----------|-------|-------------|
| `gw_reservoir_max_storage_m`<sup>\*</sup> | `gw_max_storage_m` | m | Maximum GW storage |
| `gw_discharge_coeff_m_per_timestep`<sup>\*</sup> | `gw_discharge_coefficient` | m s-1 (internal) | GW discharge coefficient |
| `gw_discharge_exponent`<sup>\*</sup> | `gw_discharge_exponent` | - | GW discharge exponent |
| `state_gw_reservoir_init_storage_m` | — | m | Initial GW storage |

### GIUH Surface Routing

| Config Key | Units | Description |
|-----------|-------|-------------|
| `surface_routing_num_giuh_ordinates` | — | Number of GIUH ordinates (required) |
| `surface_routing_giuh_ordinates` | - | Comma-separated ordinates summing to 1.0 |
| `state_surface_routing_init_giuh_convolution_queue_m` | m | Initial convolution queue (one per ordinate) |

GIUH is the only surface routing method in CFE v3. Nash Cascade surface routing
was removed as it did not produce added model skill.

### Subsurface Routing (Nash Cascade)

| Config Key | BMI Name | Units | Description |
|-----------|----------|-------|-------------|
| `subsurface_routing_nash_reservoir_time_constant_k`<sup>\*</sup> | `subsurface_nash_K` | h-1 | Subsurface Nash time constant |
| `state_subsurface_routing_init_nash_cascade_storage_m` | — | m | Initial Nash cascade storage (2 values) |

### Xinanjiang Parameters (when `partitioning_scheme_name=XINANJIANG`)

| Config Key | BMI Name | Units |
|-----------|----------|-------|
| `partitioning_Xinanjiang_tension_water_inflection_point`<sup>\*</sup> | `Xinanjiang_inflection_a` | - |
| `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent`<sup>\*</sup> | `Xinanjiang_shape_b` | - |
| `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent`<sup>\*</sup> | `Xinanjiang_shape_x` | - |

### DSBM and v3 Feature Controls

| Config Key | BMI Name | Description |
|-----------|----------|-------------|
| `control_ET_simulate_Priestley_Taylor` | `Priestley_Taylor_alpha` | P-T alpha coefficient (0 = disabled) |
| `control_soil_simulate_discrete_soil_moisture_true_false` | — | `TRUE` to enable DSBM |
| `control_soil_use_lookup_table_num_points` | — | LUT points (0 = analytic) |
| `control_ET_deepest_root_zone_discretization` | — | Deepest root zone layer (1-4) |
| `soil_ice_content_impervious_threshold` | `soil_ice_imperv_threshold` | Ice fraction threshold |

## Infiltration Excess Partitioning

1. **Schaake** (`partitioning_scheme_name=SCHAAKE`) — no additional parameters
2. **Xinanjiang** (`partitioning_scheme_name=XINANJIANG`) — requires the three Xinanjiang parameters above
