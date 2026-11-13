// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

use alloc::vec;
use alloc::vec::Vec;

use crate::cycle::Grid;

pub struct Drainage {
    pub inflow_mm_per_year: Vec<f64>,
    pub basin_salt_per_year: Vec<f64>,
    pub sea_discharge_mm_per_year: f64,
    pub sea_salt_per_year: f64,
}

pub fn drain(grid: &Grid, runoff_mm_per_year: &[f64], solubility: f64) -> Drainage {
    let cells = grid.len();
    let mut land: Vec<usize> = (0..cells).filter(|&i| !grid.ocean[i]).collect();
    land.sort_by(|&a, &b| {
        grid.elevation[b]
            .total_cmp(&grid.elevation[a])
            .then(a.cmp(&b))
    });

    let mut inflow = vec![0.0; cells];
    let mut salt_in = vec![0.0; cells];
    let mut basin_salt = vec![0.0; cells];
    let mut sea_discharge = 0.0;
    let mut sea_salt = 0.0;
    for i in land {
        let passing = inflow[i] + runoff_mm_per_year[i];
        let salt = salt_in[i] + solubility * runoff_mm_per_year[i];
        match grid.downhill(i) {
            Some(j) if grid.ocean[j] => {
                sea_discharge += passing * grid.area(i);
                sea_salt += salt * grid.area(i);
            }
            Some(j) => {
                inflow[j] += passing * grid.area_ratio(i, j);
                salt_in[j] += salt * grid.area_ratio(i, j);
            }
            None => basin_salt[i] += salt,
        }
    }
    Drainage {
        inflow_mm_per_year: inflow,
        basin_salt_per_year: basin_salt,
        sea_discharge_mm_per_year: sea_discharge,
        sea_salt_per_year: sea_salt,
    }
}

pub fn steady_salinity_grams_per_kg(
    sea_salt_mm_per_year: f64,
    residence_years: f64,
    total_sea_depth_metres: f64,
) -> f64 {
    if total_sea_depth_metres <= 0.0 {
        return 0.0;
    }
    sea_salt_mm_per_year * residence_years / total_sea_depth_metres
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Strip {
        elevation: Vec<f64>,
        celsius: Vec<f64>,
        ocean: Vec<bool>,
    }

    const WALL: f64 = 100.0;

    impl Strip {
        fn new(mut elevation: Vec<f64>) -> Self {
            elevation.push(WALL);
            Self {
                celsius: vec![10.0; elevation.len()],
                ocean: elevation.iter().map(|&e| e < 0.0).collect(),
                elevation,
            }
        }

        fn grid(&self) -> Grid<'_> {
            Grid {
                width: self.elevation.len() as u32,
                height: 1,
                elevation: &self.elevation,
                celsius: &self.celsius,
                condensation_celsius: &self.celsius,
                ocean: &self.ocean,
                sun_latitude: 0.0,
            }
        }
    }

    #[test]
    fn runoff_gathers_downstream_and_reaches_the_sea() {
        let strip = Strip::new(vec![-1.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        let drainage = drain(&strip.grid(), &[0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0], 0.0);
        assert_eq!(drainage.inflow_mm_per_year[1], 4.0);
        assert_eq!(drainage.sea_discharge_mm_per_year, 5.0);
    }

    #[test]
    fn runoff_dissolves_rock_once_and_rivers_deliver_it_to_the_sea() {
        let strip = Strip::new(vec![-1.0, 1.0, 2.0, 3.0]);
        let drainage = drain(&strip.grid(), &[0.0, 2.0, 2.0, 2.0, 0.0], 0.5);
        assert_eq!(drainage.sea_salt_per_year, 0.5 * (2.0 + 2.0 + 2.0));
        assert_eq!(
            drainage.sea_salt_per_year,
            0.5 * drainage.sea_discharge_mm_per_year
        );
    }

    #[test]
    fn closed_basins_keep_their_salt() {
        let strip = Strip::new(vec![-1.0, 3.0, 1.0, 3.0, 5.0]);
        let drainage = drain(&strip.grid(), &[0.0, 0.0, 1.0, 1.0, 1.0, 0.0], 0.5);
        assert_eq!(drainage.sea_salt_per_year, 0.0);
        assert_eq!(drainage.basin_salt_per_year[2], 0.5 * (1.0 + 1.0 + 1.0));
    }

    #[test]
    fn salinity_settles_where_input_balances_removal() {
        assert_eq!(steady_salinity_grams_per_kg(2.0, 1000.0, 500.0), 4.0);
    }
}
