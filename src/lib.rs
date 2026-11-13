// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![cfg_attr(not(test), no_std)]

extern crate alloc;

pub mod atlas;
pub mod atmosphere;
pub mod climate;
pub mod cycle;
pub mod energy;
pub mod erosion;
pub mod hash;
pub mod hydrology;
pub mod life;
pub mod material;
pub mod noise;
pub mod ocean;
pub mod parallel;
pub mod radiation;
pub mod relief;
pub mod scale;
pub mod sky;
pub mod wind;
pub mod world;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Seed(pub u64);
