// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::string::String;
use alloc::vec::Vec;
use core::f64::consts::{PI, TAU};

use serde::Deserialize;

pub const GRAVITATIONAL_CONSTANT: f64 = 6.674_30e-11;

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Star {
    pub name: String,
    pub luminosity_watts: f64,
    pub mass_kg: f64,
}

impl Star {
    pub fn problems(&self) -> Vec<String> {
        let mut problems = Vec::new();
        if self.luminosity_watts <= 0.0 {
            problems.push(alloc::format!(
                "star `{}`: luminosity_watts must be positive",
                self.name
            ));
        }
        if self.mass_kg <= 0.0 {
            problems.push(alloc::format!(
                "star `{}`: mass_kg must be positive",
                self.name
            ));
        }
        problems
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    pub semi_major_axis_metres: f64,
    pub eccentricity: f64,
    pub axial_tilt_degrees: f64,
    pub sun_longitude_at_perihelion_degrees: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sunlight {
    pub flux_w_per_m2: f64,
    pub declination: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Companion {
    pub name: String,
    pub luminosity_watts: f64,
    pub distance_metres: f64,
    pub ecliptic_latitude_degrees: f64,
    pub ecliptic_longitude_degrees: f64,
    #[serde(default)]
    pub aether_watts: f64,
}

impl Companion {
    pub fn problems(&self) -> Vec<String> {
        let mut problems = Vec::new();
        if self.aether_watts < 0.0 {
            problems.push(alloc::format!(
                "companion `{}`: aether_watts must not be negative",
                self.name
            ));
        }
        for (value, field) in [
            (self.luminosity_watts, "luminosity_watts"),
            (self.distance_metres, "distance_metres"),
        ] {
            if value <= 0.0 {
                problems.push(alloc::format!(
                    "companion `{}`: {field} must be positive",
                    self.name
                ));
            }
        }
        if self.ecliptic_latitude_degrees.abs() > 90.0 {
            problems.push(alloc::format!(
                "companion `{}`: ecliptic_latitude_degrees must lie within ±90",
                self.name
            ));
        }
        problems
    }
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Heavens {
    pub star: Star,
    #[serde(default)]
    pub companions: Vec<Companion>,
}

impl Heavens {
    pub fn problems(&self) -> Vec<String> {
        let mut problems = self.star.problems();
        for companion in &self.companions {
            problems.extend(companion.problems());
        }
        problems
    }

    pub fn sky(&self, orbit: Orbit) -> Sky<'_> {
        Sky {
            star: &self.star,
            orbit,
            companions: &self.companions,
        }
    }
}

pub struct Sky<'a> {
    pub star: &'a Star,
    pub orbit: Orbit,
    pub companions: &'a [Companion],
}

struct Position {
    distance_metres: f64,
    sun_longitude: f64,
}

impl Sky<'_> {
    pub fn year_seconds(&self) -> f64 {
        let a = self.orbit.semi_major_axis_metres;
        TAU * libm::sqrt(a * a * a / (GRAVITATIONAL_CONSTANT * self.star.mass_kg))
    }

    fn position(&self, seconds: f64) -> Position {
        let orbit = &self.orbit;
        let e = orbit.eccentricity;
        let turns = seconds / self.year_seconds();
        let mean_anomaly = TAU * (turns - libm::floor(turns));
        let eccentric = eccentric_anomaly(mean_anomaly, e);
        let true_anomaly = 2.0
            * libm::atan2(
                libm::sqrt(1.0 + e) * libm::sin(eccentric / 2.0),
                libm::sqrt(1.0 - e) * libm::cos(eccentric / 2.0),
            );
        Position {
            distance_metres: orbit.semi_major_axis_metres * (1.0 - e * libm::cos(eccentric)),
            sun_longitude: true_anomaly + orbit.sun_longitude_at_perihelion_degrees.to_radians(),
        }
    }

    pub fn at(&self, seconds: f64) -> Sunlight {
        let Position {
            distance_metres: distance,
            sun_longitude,
        } = self.position(seconds);
        Sunlight {
            flux_w_per_m2: self.star.luminosity_watts / (2.0 * TAU * distance * distance),
            declination: libm::asin(
                libm::sin(self.orbit.axial_tilt_degrees.to_radians()) * libm::sin(sun_longitude),
            ),
        }
    }

    pub fn companion_light(&self, companion: &Companion, seconds: f64) -> Sunlight {
        let here = self.position(seconds);
        let planet_longitude = here.sun_longitude + PI;
        let planet = [
            here.distance_metres * libm::cos(planet_longitude),
            here.distance_metres * libm::sin(planet_longitude),
            0.0,
        ];
        let (latitude, longitude) = (
            companion.ecliptic_latitude_degrees.to_radians(),
            companion.ecliptic_longitude_degrees.to_radians(),
        );
        let far = companion.distance_metres;
        let toward = [
            far * libm::cos(latitude) * libm::cos(longitude) - planet[0],
            far * libm::cos(latitude) * libm::sin(longitude) - planet[1],
            far * libm::sin(latitude) - planet[2],
        ];
        let distance =
            libm::sqrt(toward[0] * toward[0] + toward[1] * toward[1] + toward[2] * toward[2]);
        let tilt = self.orbit.axial_tilt_degrees.to_radians();
        let along_axis = (toward[2] * libm::cos(tilt) + toward[1] * libm::sin(tilt)) / distance;
        Sunlight {
            flux_w_per_m2: companion.luminosity_watts / (2.0 * TAU * distance * distance),
            declination: libm::asin(along_axis.clamp(-1.0, 1.0)),
        }
    }

    pub fn daily_aether(&self, seconds: f64, latitude_north: f64) -> f64 {
        self.companions
            .iter()
            .filter(|c| c.aether_watts > 0.0)
            .map(|c| {
                let light = self.companion_light(c, seconds);
                let aether = Sunlight {
                    flux_w_per_m2: light.flux_w_per_m2 * c.aether_watts / c.luminosity_watts,
                    ..light
                };
                daily_insolation(aether, latitude_north)
            })
            .sum()
    }

    pub fn daily_insolation(&self, seconds: f64, latitude_north: f64) -> f64 {
        daily_insolation(self.at(seconds), latitude_north)
            + self
                .companions
                .iter()
                .map(|c| daily_insolation(self.companion_light(c, seconds), latitude_north))
                .sum::<f64>()
    }
}

fn eccentric_anomaly(mean_anomaly: f64, eccentricity: f64) -> f64 {
    let mut eccentric = mean_anomaly;
    let mut last_step = f64::INFINITY;
    loop {
        let step = (eccentric - eccentricity * libm::sin(eccentric) - mean_anomaly)
            / (1.0 - eccentricity * libm::cos(eccentric));
        if step == 0.0 || step.abs() >= last_step {
            return eccentric;
        }
        eccentric -= step;
        last_step = step.abs();
    }
}

pub fn daily_insolation(sunlight: Sunlight, latitude_north: f64) -> f64 {
    let (sin_lat, cos_lat) = (libm::sin(latitude_north), libm::cos(latitude_north));
    let (sin_dec, cos_dec) = (
        libm::sin(sunlight.declination),
        libm::cos(sunlight.declination),
    );
    let sunset_hour_angle =
        libm::acos((-(sin_lat * sin_dec) / (cos_lat * cos_dec)).clamp(-1.0, 1.0));
    sunlight.flux_w_per_m2 / PI
        * (sunset_hour_angle * sin_lat * sin_dec + cos_lat * cos_dec * libm::sin(sunset_hour_angle))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::earth::{self, sun};

    const DAY: f64 = 86_400.0;

    fn earth(star: &Star) -> Sky<'_> {
        Sky {
            star,
            orbit: earth::ORBIT,
            companions: &[],
        }
    }

    fn companion(latitude_degrees: f64, distance_metres: f64, luminosity_watts: f64) -> Companion {
        Companion {
            name: "ember".into(),
            luminosity_watts,
            distance_metres,
            ecliptic_latitude_degrees: latitude_degrees,
            ecliptic_longitude_degrees: 0.0,
            aether_watts: 0.0,
        }
    }

    #[test]
    fn aether_light_reaches_only_what_absorbs_it() {
        let star = sun();
        let plain = [companion(-60.0, 400.0 * AU, 1.0e3 * star.luminosity_watts)];
        let glowing = [Companion {
            aether_watts: 1.0e4 * star.luminosity_watts,
            ..plain[0].clone()
        }];
        let sky = |companions| Sky {
            companions,
            ..earth(&star)
        };
        let (plain, glowing) = (sky(&plain), sky(&glowing));
        let south = -70.0_f64.to_radians();
        for day in [0, 91, 182, 273] {
            let seconds = f64::from(day) * DAY;
            assert_eq!(
                plain.daily_insolation(seconds, south),
                glowing.daily_insolation(seconds, south)
            );
            assert_eq!(plain.daily_aether(seconds, south), 0.0);
            let ratio = glowing.daily_aether(seconds, south)
                / (glowing.daily_insolation(seconds, south)
                    - earth(&star).daily_insolation(seconds, south));
            assert!((ratio - 10.0).abs() < 1e-9, "day {day}: {ratio}");
        }
    }

    #[test]
    fn a_companion_where_the_sun_stands_shines_exactly_like_the_sun() {
        let star = sun();
        let sky = earth(&star);
        let twin = companion(0.0, 0.0, star.luminosity_watts);
        for day in 0..365 {
            let seconds = f64::from(day) * DAY;
            let (sun, twin) = (sky.at(seconds), sky.companion_light(&twin, seconds));
            assert!(
                close(sun.flux_w_per_m2, twin.flux_w_per_m2, 1e-9),
                "day {day}"
            );
            assert!(close(sun.declination, twin.declination, 1e-12), "day {day}");
        }
    }

    #[test]
    fn a_companion_below_the_orbit_lights_the_south_pole_all_year_and_never_the_north() {
        let star = sun();
        let ember = [companion(-90.0, 400.0 * AU, 1.0e3 * star.luminosity_watts)];
        let sky = Sky {
            companions: &ember,
            ..earth(&star)
        };
        let north = 80.0_f64.to_radians();
        for day in 0..365 {
            let light = sky.companion_light(&ember[0], f64::from(day) * DAY);
            assert!(daily_insolation(light, -north) > 0.0, "day {day}");
            assert_eq!(daily_insolation(light, north), 0.0, "day {day}");
        }
    }

    #[test]
    fn heavens_reject_a_companion_without_light_or_distance() {
        let heavens = Heavens {
            star: sun(),
            companions: alloc::vec![companion(-60.0, 0.0, 0.0)],
        };
        assert_eq!(heavens.problems().len(), 2);
    }

    const AU: f64 = 1.495_978_707e11;

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance
    }

    #[test]
    fn the_year_follows_from_keplers_third_law() {
        let sun = sun();
        let days = earth(&sun).year_seconds() / DAY;
        assert!(close(days, 365.256, 0.01), "{days}");
    }

    #[test]
    fn sunlight_at_earth_matches_the_measured_solar_constant() {
        let sun = sun();
        let sky = earth(&sun);
        let year = sky.year_seconds();
        let mean = (0..365)
            .map(|d| sky.at(f64::from(d) / 365.0 * year).flux_w_per_m2)
            .sum::<f64>()
            / 365.0;
        assert!(close(mean, 1361.0, 3.0), "{mean}");
        let perihelion = sky.at(0.0).flux_w_per_m2;
        let aphelion = sky.at(year / 2.0).flux_w_per_m2;
        assert!(
            close(perihelion / aphelion, 1.069, 0.002),
            "{perihelion} {aphelion}"
        );
    }

    #[test]
    fn the_sun_reaches_the_tropics_at_the_solstices() {
        let sun = sun();
        let sky = earth(&sun);
        let year = sky.year_seconds();
        let declinations: Vec<f64> = (0..3650)
            .map(|d| {
                sky.at(f64::from(d) / 3650.0 * year)
                    .declination
                    .to_degrees()
            })
            .collect();
        let highest = declinations.iter().copied().fold(f64::MIN, f64::max);
        let lowest = declinations.iter().copied().fold(f64::MAX, f64::min);
        assert!(close(highest, 23.44, 0.01) && close(lowest, -23.44, 0.01));
    }

    #[test]
    fn polar_night_gets_no_sun_and_polar_summer_gets_more_than_the_equator() {
        let solstice = Sunlight {
            flux_w_per_m2: 1361.0,
            declination: 23.44_f64.to_radians(),
        };
        let north_pole = 89.9_f64.to_radians();
        assert_eq!(daily_insolation(solstice, -north_pole), 0.0);
        let polar_day = daily_insolation(solstice, north_pole);
        let equator = daily_insolation(solstice, 0.0);
        assert!(polar_day > equator, "pole {polar_day} equator {equator}");
        assert!(close(polar_day, 1361.0 * 0.3978, 2.0), "{polar_day}");
    }

    #[test]
    fn an_equinox_lights_the_equator_by_the_solar_constant_over_pi() {
        let equinox = Sunlight {
            flux_w_per_m2: 1361.0,
            declination: 0.0,
        };
        assert!(close(daily_insolation(equinox, 0.0), 1361.0 / PI, 1e-9));
    }
}
