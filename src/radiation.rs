use crate::atmosphere::{LAYERS, RADIATION_SLABS, SLABS, TROPOSPHERE, slab_air};
use crate::material::Material;
use crate::scale::{KELVIN_AT_ZERO_CELSIUS, STEFAN_BOLTZMANN_W_PER_M2_PER_K4};

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Longwave {
    pub to_space_w_per_m2: f64,
    pub down_to_ground_w_per_m2: f64,
    pub slab_heating_w_per_m2: [f64; SLABS],
}

impl Longwave {
    pub fn air_cooling_w_per_m2(&self) -> f64 {
        -self.slab_heating_w_per_m2.iter().sum::<f64>()
    }

    pub fn layer_heating_w_per_m2(&self) -> [f64; LAYERS] {
        let mut layers = [0.0; LAYERS];
        for (slab, heating) in RADIATION_SLABS.iter().zip(&self.slab_heating_w_per_m2) {
            layers[slab.layer] += heating;
        }
        layers
    }
}

pub fn glow(kelvin: f64) -> f64 {
    STEFAN_BOLTZMANN_W_PER_M2_PER_K4 * kelvin * kelvin * kelvin * kelvin
}

pub fn capacity_celsius(air_celsius: f64, layer_celsius: &[f64; LAYERS]) -> [f64; LAYERS] {
    core::array::from_fn(|layer| {
        if layer == 0 {
            air_celsius
        } else {
            layer_celsius[layer]
        }
    })
}

pub fn slab_vapor_kg_per_m2(
    water: &Material,
    air_celsius: f64,
    slab_saturation: &[f64; SLABS],
    layer_humidity: &[f64; LAYERS],
) -> [f64; SLABS] {
    let weight = |slab: usize| {
        let at = RADIATION_SLABS[slab];
        if at.layer == 0 {
            at.air_kg_per_m2()
        } else {
            slab_saturation[slab] * at.thickness_metres()
        }
    };
    let mut layer_weight = [0.0; LAYERS];
    for slab in 0..SLABS {
        layer_weight[RADIATION_SLABS[slab].layer] += weight(slab);
    }
    let held: [f64; LAYERS] = core::array::from_fn(|layer| {
        let relative = layer_humidity[layer].max(0.0);
        if layer == 0 {
            relative * water.saturation(air_celsius) * TROPOSPHERE[0].thickness_metres()
        } else {
            relative * layer_weight[layer]
        }
    });
    core::array::from_fn(|slab| {
        let layer = RADIATION_SLABS[slab].layer;
        held[layer] * weight(slab) / layer_weight[layer]
    })
}

pub fn longwave(
    ground_celsius: f64,
    slab_celsius: &[f64; SLABS],
    slab_vapor_kg_per_m2: &[f64; SLABS],
    air: &Material,
    water: &Material,
) -> Longwave {
    let absorbing = 1.0 - air.longwave_window_share;
    let emissivity: [f64; SLABS] = core::array::from_fn(|slab| {
        let depth = air.longwave_absorption_m2_per_kg * RADIATION_SLABS[slab].air_kg_per_m2()
            + water.longwave_absorption_m2_per_kg * slab_vapor_kg_per_m2[slab];
        1.0 - libm::exp(-depth)
    });
    let emitted: [f64; SLABS] =
        core::array::from_fn(|slab| glow(slab_celsius[slab] + KELVIN_AT_ZERO_CELSIUS));
    let ground_glow = glow(ground_celsius + KELVIN_AT_ZERO_CELSIUS);
    let mut up_into = [0.0; SLABS];
    let mut up = ground_glow;
    for slab in 0..SLABS {
        up_into[slab] = up;
        up = up * (1.0 - emissivity[slab]) + emissivity[slab] * emitted[slab];
    }
    let mut down_into = [0.0; SLABS];
    let mut down = 0.0;
    for slab in (0..SLABS).rev() {
        down_into[slab] = down;
        down = down * (1.0 - emissivity[slab]) + emissivity[slab] * emitted[slab];
    }
    Longwave {
        to_space_w_per_m2: absorbing * up + air.longwave_window_share * ground_glow,
        down_to_ground_w_per_m2: absorbing * down,
        slab_heating_w_per_m2: core::array::from_fn(|slab| {
            absorbing * emissivity[slab] * (up_into[slab] + down_into[slab] - 2.0 * emitted[slab])
        }),
    }
}

pub struct Column {
    pub slab_celsius: [f64; SLABS],
    pub slab_saturation: [f64; SLABS],
    pub parcels: [f64; SLABS],
    pub longwave: Longwave,
}

pub fn column(
    ground_celsius: f64,
    air_celsius: f64,
    layer_humidity: &[f64; LAYERS],
    air: &Material,
    water: &Material,
) -> Column {
    column_from(
        None,
        ground_celsius,
        air_celsius,
        layer_humidity,
        air,
        water,
    )
}

pub fn column_from(
    parcel_guesses: Option<&[f64; SLABS]>,
    ground_celsius: f64,
    air_celsius: f64,
    layer_humidity: &[f64; LAYERS],
    air: &Material,
    water: &Material,
) -> Column {
    let slabs = slab_air(parcel_guesses, air_celsius, layer_humidity[0], air, water);
    let vapor = slab_vapor_kg_per_m2(water, air_celsius, &slabs.saturation, layer_humidity);
    Column {
        longwave: longwave(ground_celsius, &slabs.celsius, &vapor, air, water),
        slab_celsius: slabs.celsius,
        slab_saturation: slabs.saturation,
        parcels: slabs.celsius,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{air, water};

    const EARTH_CELSIUS: f64 = 288.0 - KELVIN_AT_ZERO_CELSIUS;
    const EARTH_TO_SPACE_W_PER_M2: f64 = 239.0;
    const EARTH_DOWN_TO_GROUND_W_PER_M2: f64 = 342.0;
    const EARTH_WINDOW_W_PER_M2: f64 = 40.0;
    pub(crate) const EARTH_LAYER_HUMIDITY: [f64; LAYERS] = [0.8, 0.5, 0.4];

    fn earth(ground_celsius: f64, air_celsius: f64, layer_humidity: [f64; LAYERS]) -> Longwave {
        column(
            ground_celsius,
            air_celsius,
            &layer_humidity,
            &air(),
            &water(),
        )
        .longwave
    }

    #[test]
    fn clear_air_lets_the_ground_glow_straight_to_space() {
        let mut clear = air();
        clear.longwave_absorption_m2_per_kg = 0.0;
        let longwave = longwave(
            EARTH_CELSIUS,
            &[5.0; SLABS],
            &[0.0; SLABS],
            &clear,
            &water(),
        );
        let ground = glow(EARTH_CELSIUS + KELVIN_AT_ZERO_CELSIUS);
        assert!((longwave.to_space_w_per_m2 - ground).abs() < 1e-9 * ground);
        assert_eq!(longwave.down_to_ground_w_per_m2, 0.0);
        assert_eq!(longwave.air_cooling_w_per_m2(), 0.0);
    }

    #[test]
    fn the_window_band_carries_the_ground_glow_past_even_opaque_air() {
        let mut opaque = air();
        opaque.longwave_absorption_m2_per_kg = 1.0;
        let cold = longwave(
            EARTH_CELSIUS,
            &[-60.0; SLABS],
            &[0.0; SLABS],
            &opaque,
            &water(),
        );
        let window = opaque.longwave_window_share * glow(EARTH_CELSIUS + KELVIN_AT_ZERO_CELSIUS);
        assert!(cold.to_space_w_per_m2 > window);
        assert!(
            (window - EARTH_WINDOW_W_PER_M2).abs() < 0.1 * EARTH_WINDOW_W_PER_M2,
            "{window}"
        );
    }

    #[test]
    fn the_air_cools_by_what_leaves_the_top_minus_what_the_ground_loses() {
        let longwave = earth(20.0, 14.0, [0.8; LAYERS]);
        let ground_loses = glow(20.0 + KELVIN_AT_ZERO_CELSIUS) - longwave.down_to_ground_w_per_m2;
        let expected = longwave.to_space_w_per_m2 - ground_loses;
        assert!(
            (longwave.air_cooling_w_per_m2() - expected).abs() < 1e-9,
            "{} vs {expected}",
            longwave.air_cooling_w_per_m2()
        );
    }

    #[test]
    fn earthlike_air_sends_earths_glow_to_space_and_back_to_the_ground() {
        let longwave = earth(EARTH_CELSIUS, EARTH_CELSIUS, EARTH_LAYER_HUMIDITY);
        let near = |value: f64, earth: f64| (value - earth).abs() < 0.03 * earth;
        assert!(
            near(longwave.to_space_w_per_m2, EARTH_TO_SPACE_W_PER_M2),
            "{}",
            longwave.to_space_w_per_m2
        );
        assert!(
            near(
                longwave.down_to_ground_w_per_m2,
                EARTH_DOWN_TO_GROUND_W_PER_M2
            ),
            "{}",
            longwave.down_to_ground_w_per_m2
        );
    }

    #[test]
    fn earthlike_air_cools_most_above_the_boundary_layer() {
        let cooling = earth(16.0, 14.2, EARTH_LAYER_HUMIDITY).layer_heating_w_per_m2();
        let boundary = -cooling[0];
        let above = -(cooling[1] + cooling[2]);
        assert!(boundary > 0.0 && above > 2.0 * boundary, "{cooling:?}");
    }

    #[test]
    fn moister_air_holds_more_glow_in_and_sends_more_back_down() {
        let (dry, moist) = (
            earth(EARTH_CELSIUS, EARTH_CELSIUS, [0.3; LAYERS]),
            earth(EARTH_CELSIUS, EARTH_CELSIUS, [0.9; LAYERS]),
        );
        assert!(moist.to_space_w_per_m2 < dry.to_space_w_per_m2);
        assert!(moist.down_to_ground_w_per_m2 > dry.down_to_ground_w_per_m2);
    }

    #[test]
    fn every_slab_stays_thin_in_the_moistest_column() {
        let water = water();
        let air = air();
        let tropical = 30.0;
        let lapse = crate::scale::LAPSE_RATE_CELSIUS_PER_METRE;
        let fading_per_metre =
            water.saturation_per_kelvin(tropical) / water.saturation(tropical) * lapse;
        let vapor_below =
            |metres: f64| water.saturation(tropical - lapse * metres) / fading_per_metre;
        for slab in RADIATION_SLABS {
            let vapor = vapor_below(slab.base_metres) - vapor_below(slab.top_metres);
            let depth = air.longwave_absorption_m2_per_kg * slab.air_kg_per_m2()
                + water.longwave_absorption_m2_per_kg * vapor;
            assert!(depth <= 1.0, "{slab:?}: {depth}");
        }
    }

    #[test]
    fn slabs_tile_the_water_layers() {
        for (layer, water) in TROPOSPHERE.iter().enumerate() {
            let mut inside = RADIATION_SLABS.iter().filter(|slab| slab.layer == layer);
            let first = inside.next().map(|slab| slab.base_metres);
            let last = inside.next_back().map(|slab| slab.top_metres);
            assert_eq!(first, Some(water.base_metres));
            assert_eq!(last, Some(water.top_metres));
        }
        for pair in RADIATION_SLABS.windows(2) {
            assert_eq!(pair[0].top_metres, pair[1].base_metres);
            assert_eq!(pair[0].top_hpa, pair[1].base_hpa);
        }
    }
}
