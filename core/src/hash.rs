// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use crate::Seed;

#[derive(Clone, Copy)]
pub enum Purpose {
    Elevation,
    Weather,
    Look,
    Relief,
    Sway,
    Wild,
    Fauna,
    Rain,
    Luck,
}

impl Purpose {
    fn id(self) -> u64 {
        match self {
            Self::Elevation => 1,
            Self::Weather => 2,
            Self::Look => 3,
            Self::Relief => 4,
            Self::Sway => 5,
            Self::Wild => 6,
            Self::Fauna => 7,
            Self::Rain => 8,
            Self::Luck => 9,
        }
    }
}

pub fn lattice(seed: Seed, purpose: Purpose, layer: u64, x: i64, y: i64, z: i64) -> u64 {
    let mut h = mix(seed.0 ^ mix(purpose.id()));
    for value in [layer, x as u64, y as u64, z as u64] {
        h = mix(h ^ value);
    }
    h
}

pub fn unit(h: u64) -> f64 {
    (h >> 11) as f64 / (1u64 << 53) as f64
}

fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}
