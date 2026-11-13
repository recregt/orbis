// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use core::f64::consts::TAU;

use crate::Seed;
use crate::hash::{Purpose, lattice, unit};

pub struct Field {
    pub seed: Seed,
    pub purpose: Purpose,
    pub base_frequency: f64,
    pub roughness: f64,
    pub max_layers: u64,
}

impl Field {
    pub fn sample_cylinder(&self, x: i64, y: u32, width: u32) -> f64 {
        let column = x.rem_euclid(i64::from(width)) as f64;
        let angle = column / f64::from(width) * TAU;
        let (sin, cos) = (libm::sin(angle), libm::cos(angle));
        let row = f64::from(y) + 0.5;

        let mut frequency = self.base_frequency;
        let mut amplitude = 1.0;
        let mut sum = 0.0;
        let mut norm = 0.0;
        let mut layer = 0;
        while f64::from(width) / frequency >= 2.0 && layer < self.max_layers {
            let radius = frequency / TAU;
            let z = row / f64::from(width) * frequency;
            sum += amplitude * self.value(layer, radius * cos, radius * sin, z);
            norm += amplitude;
            amplitude *= self.roughness;
            frequency *= 2.0;
            layer += 1;
        }
        if norm == 0.0 { 0.0 } else { sum / norm }
    }

    fn value(&self, layer: u64, x: f64, y: f64, z: f64) -> f64 {
        let (x0, y0, z0) = (libm::floor(x), libm::floor(y), libm::floor(z));
        let (fx, fy, fz) = (fade(x - x0), fade(y - y0), fade(z - z0));
        let (ix, iy, iz) = (x0 as i64, y0 as i64, z0 as i64);
        let corner = |dx: i64, dy: i64, dz: i64| {
            unit(lattice(
                self.seed,
                self.purpose,
                layer,
                ix + dx,
                iy + dy,
                iz + dz,
            )) * 2.0
                - 1.0
        };
        let face = |dz: i64| {
            let bottom = lerp(corner(0, 0, dz), corner(1, 0, dz), fx);
            let top = lerp(corner(0, 1, dz), corner(1, 1, dz), fx);
            lerp(bottom, top, fy)
        };
        lerp(face(0), face(1), fz)
    }
}

fn fade(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
