# CFE Configuration Files

CFE v3 supports two config file formats. Example configurations are provided
in this directory.

## Config Formats

### Legacy v2 format (`.cf2` / `.txt`)

Simple `key=value[units]` format, one parameter per line. Lines starting with
`#` are comments. Example: `legacy_cfe_config_cat87.cf2`.

### v3 keyword format (`.cf3`)

Extended `keyword=arg(s) [units]` format with `#` and `//` comments, richer
keyword names, and support for v3 features (DSBM, Priestley-Taylor PET,
output configuration). Example: `bmi_config_cat87_v3.cf3`.

See `clean_config.cf3` for a fully annotated template of all v3 options.

## Example Configs

| File | Format | Partitioning | Routing | Soil | Notes |
|------|--------|-------------|---------|------|-------|
| `cfe_config_cat_87_pass.txt` | v2 | Schaake | Nash Cascade | Linear reservoir | Used by unit tests |
| `legacy_cfe_config_cat87.cf2` | v2 | Xinanjiang | Nash Cascade | Linear reservoir | Used by integration test |
| `bmi_config_cat87_v3.cf3` | v3 | Schaake | GIUH | DSBM (4-layer) | Used by integration test |
| `cfe_config_cat87_v3.cf3` | v3 | Schaake | GIUH | DSBM (4-layer) | Standalone driver config |
| `clean_config.cf3` | v3 | — | — | — | Annotated template |

## Parameter Reference

Parameters marked with <sup>\*</sup> are calibratable via BMI `set_value` /
`get_value_ptr`. The v3 BMI canonical name and v2 config key are shown;
see `CHANGELOG.md` for the full v2→v3 name mapping.

### Soil Parameters

| v2 Config Key | v3 Config Key | BMI Name | Units | Description |
|--------------|--------------|----------|-------|-------------|
| `soil_params.depth` | `soil_depth_m` | `param_soil_depth_m` | m | Soil column depth |
| `soil_params.b`<sup>\*</sup> | `soil_Clapp_Hornberger_exponent_b` | `soil_Clapp_Hornberger_b` | - | Clapp-Hornberger exponent |
| `soil_params.satdk`<sup>\*</sup> | `soil_sat_hydraulic_conductivity_cm_per_h` | `soil_saturated_hydraulic_conductivity` | m s-1 (internal) | Saturated hydraulic conductivity |
| `soil_params.satpsi`<sup>\*</sup> | `soil_sat_capillary_head_cm` | `soil_saturated_capillary_head` | m (internal) | Saturated capillary head |
| `soil_params.smcmax`<sup>\*</sup> | `soil_effective_porosity` | `soil_effective_porosity` | - | Effective porosity |
| `soil_params.wltsmc`<sup>\*</sup> | `soil_wilting_point_moisture_content` | `soil_wilting_point` | - | Wilting point |
| `soil_params.slop`<sup>\*</sup> | `soil_to_gw_percolation_rate_limiter_0_to_1` | `soil_percolation_rate_limiter` | - | Percolation rate limiter (0-1) |
| `alpha_fc`<sup>\*</sup> | `soil_field_capacity_Pcap_over_Patm_0_1` | `soil_field_capacity_fraction` | - | Field capacity (Pcap/Patm) |
| `K_lf`<sup>\*</sup> | `soil_reservoir_rate_const_to_subsurface_lateral_flow` | `soil_lateral_flow_K` | h-1 | Lateral flow rate constant |
| `soil_storage` | `state_soil_reservoir_init_storage_m` | — | m | Initial soil storage |

### Groundwater Parameters

| v2 Config Key | v3 Config Key | BMI Name | Units | Description |
|--------------|--------------|----------|-------|-------------|
| `max_gw_storage`<sup>\*</sup> | `gw_reservoir_max_storage_m` | `gw_max_storage_m` | m | Maximum GW storage |
| `Cgw`<sup>\*</sup> | `gw_discharge_coeff_m_per_timestep` | `gw_discharge_coefficient` | m s-1 (internal) | GW discharge coefficient |
| `expon`<sup>\*</sup> | `gw_discharge_exponent` | `gw_discharge_exponent` | - | GW discharge exponent |
| `gw_storage` | `state_gw_reservoir_init_storage_m` | — | m | Initial GW storage |

### GIUH Surface Routing Parameters

| v3 Config Key | Units | Description |
|--------------|-------|-------------|
| `surface_routing_num_giuh_ordinates` | — | Number of GIUH ordinates (required) |
| `surface_routing_giuh_ordinates` | - | Comma-separated ordinates summing to 1.0 |
| `state_surface_routing_init_giuh_convolution_queue_m` | m | Initial convolution queue (one per ordinate) |

> **Deprecated:** `surface_runoff_scheme`, `N_nash_surface`, `K_nash_surface`,
> `Kinf_nash_surface`, `retention_depth_nash_surface` — accepted in legacy
> configs with a warning, but Nash Cascade surface routing is no longer available.

### Subsurface Routing Parameters

| v2 Config Key | v3 Config Key | BMI Name | Units | Description |
|--------------|--------------|----------|-------|-------------|
| `K_nash_subsurface`<sup>\*</sup> | `subsurface_routing_nash_reservoir_time_constant_k` | `subsurface_nash_K` | h-1 | Subsurface Nash time constant |

### Xinanjiang Parameters (when `surface_water_partitioning_scheme=Xinanjiang`)

| v2 Config Key | v3 Config Key | BMI Name | Units |
|--------------|--------------|----------|-------|
| `a_Xinanjiang_inflection_point_parameter`<sup>\*</sup> | `partitioning_Xinanjiang_tension_water_inflection_point` | `Xinanjiang_inflection_a` | - |
| `b_Xinanjiang_shape_parameter`<sup>\*</sup> | `partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent` | `Xinanjiang_shape_b` | - |
| `x_Xinanjiang_shape_parameter`<sup>\*</sup> | `partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent` | `Xinanjiang_shape_x` | - |

### v3-Only Parameters

| v3 Config Key | BMI Name | Units | Description |
|--------------|----------|-------|-------------|
| `control_ET_simulate_Priestley_Taylor` | `Priestley_Taylor_alpha` | - | P-T alpha coefficient (0 = disabled) |
| `control_soil_simulate_discrete_soil_moisture_true_false` | — | — | `TRUE` to enable DSBM |
| `control_soil_use_lookup_table_num_points` | — | — | LUT points (0 = analytic) |
| `control_ET_deepest_root_zone_discretization` | — | — | Deepest root zone layer (1-4) |
| `soil_ice_content_impervious_threshold` | `soil_ice_imperv_threshold` | - | Ice fraction threshold |

## Infiltration Excess Runoff Options

1. **Schaake** (`surface_water_partitioning_scheme=Schaake` or `partitioning_scheme_name=SCHAAKE`)
2. **Xinanjiang** (`surface_water_partitioning_scheme=Xinanjiang` or `partitioning_scheme_name=XINANJIANG`) — requires the three Xinanjiang parameters listed above.

## Surface Runoff Routing

CFE v3 uses **GIUH** (Geomorphological Instantaneous Unit Hydrograph) as the
only surface routing method. Requires `surface_routing_num_giuh_ordinates` and
`surface_routing_giuh_ordinates` in v3 configs.

Legacy v2 configs that specified `surface_runoff_scheme=NASH_CASCADE` are
accepted with a deprecation warning. If no GIUH ordinates are present, a
default unit impulse (1 ordinate = 1.0) is used.

> **Note:** Nash Cascade surface routing, retention depth, and runon infiltration
> were removed in CFE v3 as they did not produce added model skill.
