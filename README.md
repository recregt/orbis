# Orbis

A deterministic planet simulation. From a seed and a recipe, Orbis derives a planet's
relief, climate, oceans, rivers and vegetation, using the same equations everywhere. A seed
gives the same world bit for bit, natively and in WebAssembly.

**Live demo:** [recregt.com/tools/orbis](https://recregt.com/tools/orbis)

## What it computes

| Module | Model |
| --- | --- |
| `sky` | Keplerian orbit, axial tilt, daily-mean insolation |
| `energy` | Seasonal energy balance with ice–albedo feedback and an ocean mixed layer |
| `radiation`, `atmosphere` | Layered longwave transfer, moist-adiabatic lapse rates |
| `wind` | Hadley, Ferrel and polar cells; Hadley width from the planet's spin; surface drag |
| `ocean` | Wind-driven Sverdrup gyres with western boundary currents |
| `cycle`, `hydrology` | Evaporation, transport and rain of water and ipse; drainage; sea salinity |
| `erosion`, `relief` | Stream-power incision and hillslope creep, from planet scale to regional detail |
| `life` | Plants that grow where temperature, water and ipse allow |

Ipse is the simulation's one invented substance: it evaporates, drifts and condenses like
water, and absorbs aether, a light given off by companion bodies in the sky.

## Use from JavaScript

```sh
npm install @recregt/orbis
```

```js
import init, { Layer, Shape, World } from '@recregt/orbis'

await init()

// The ground is ready at once.
const shape = new Shape(1n, 64, 32)
shape.elevation() // Float32Array, metres per cell, row 0 at the north pole
shape.relief() // eroded height field, shape.reliefWidth × shape.reliefHeight

// The climate takes seconds to settle, so run it in a Web Worker.
const world = new World(1n, 64, 32)
world.layer(Layer.Temperature) // °C per cell
world.settled // whether the climate reached a steady state
```

`Shape` and `World` are independent, so they can run in two workers at the same time.

## Build

Orbis needs Rust 1.99 or newer. The core crate is `no_std` and depends only on `libm` and `serde`.

```sh
cargo test --workspace
cd wasm && npm run build # needs the wasm32-unknown-unknown target and wasm-bindgen-cli
```

## License

[MPL-2.0](LICENSE)
