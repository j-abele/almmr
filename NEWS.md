# almmr 1.1.0

## New features

- `create_cost_surface()` can now incorporate water routing via two optional
  arguments:
  - `rivers` (SpatVector lines): edges between two river cells are re-weighted
    by travel time on water, faster downstream (`downstream_speed_kmh`) than
    upstream (`upstream_speed_kmh`). Flow direction is derived from the river
    geometry combined with the overall elevation trend, so DEM noise between
    adjacent cells does not flip it.
  - `waterbodies` (SpatVector polygons): movement between two waterbody cells is
    isotropic, using `waterbodies_speed_kmh`.
- Water routing works in both eager and lazy mode, so all analysis functions
  (`perform_tpla()`, `lcsc_territory()`, `sbr_network()`) pick it up
  automatically.


\# almmr 1.0.0



First stable release. This version completes the transition from the

gdistance-based prototype (0.1.0) to a native implementation and commits to

stable functions from here on.



\## Breaking changes (migration from 0.1.0)



\* Cost surfaces are now `almmr\_cs` objects returned by `create\_cost\_surface()`,

&#x20; replacing gdistance `TransitionLayer` objects. All analysis functions expect

&#x20; an `almmr\_cs` (or, where supported, build one internally from the DEM).

\* `sbr\_network()`: the `conductance` argument is renamed to `cost\_surface` and

&#x20; now requires an eager `almmr\_cs`; the function returns an `sf` route network

&#x20; instead of a `SpatialLinesDataFrame`.

\* `create\_cost\_surface()` and `lcsc\_territory()`: the `epsg` argument was

&#x20; removed. Coordinate reference systems are now handled automatically.

\* `lcsc\_territory()` now returns a `SpatVector` of territory polygons instead of

&#x20; an `sp` object.

\* `perform\_tpla()` and `perform\_tpla\_line\_based()` now take an `almmr\_cs` cost

&#x20; surface and return `sf`/`SpatRaster` objects rather than `sp` geometries.



\## Major changes



\* Removed the dependencies on \*\*gdistance\*\*, \*\*raster\*\* and \*\*sp\*\*. Movement

&#x20; modelling is now implemented natively with \*\*terra\*\* and graph operations

&#x20; (\*\*igraph\*\*).

\* New exported function `compute\_lcp()`: least-cost paths with optional

&#x20; hierarchical resolution refinement and bidirectional (anisotropic) paths.

\* `create\_cost\_surface()` gained a `lazy` mode: edge weights are computed on

&#x20; demand for the required extent only, keeping memory bounded on large DEMs.

\* `create\_cost\_surface()` extended `numberOfDirections` to support 32 and 48

&#x20; neighbours, and adds automatic CRS handling and clearer input validation.



\## Performance



\* `perform\_tpla()` and `perform\_tpla\_line\_based()` now build a single graph over

&#x20; the analysis region and query it one-to-many, instead of rebuilding a graph

&#x20; for every path pair — typically one to two orders of magnitude faster. Setting

&#x20; `resolutions` (hierarchical refinement) uses the slower per-pair path.



\## Documentation and packaging



\* All README examples were updated to the new functions and now render with

&#x20; figures.

\* Added a citation file (`inst/CITATION`, `CITATION.cff`) and a

&#x20; provenance/AI-assistance note.

\* Added a `testthat` test suite covering all exported functions.

\* Relicensed from MIT to \*\*GPL-3\*\*.

