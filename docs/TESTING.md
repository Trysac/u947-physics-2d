# Testing

## Running

```sh
cmake -P scripts/check.cmake -- tests   # Catch2 tests under llvm-cov
cmake -P scripts/check.cmake            # every gate, as CI and the factory run it
```

## Coverage

- `tests/physics_test.cpp`: Covers collision logic for circles and polygons (Voronoi regions), elastic collision conservation (energy/momentum), free-fall integration, resting/sleeping behavior, distance joint maintenance, and validation errors (invalid polygons, body counts, etc.).
- `tests/text_test.cpp`: Tests basic string normalization logic (trimming, case folding) and error handling for blank inputs.
- `tests/text_properties_test.cpp`: Property-based tests using Catch2 generators to verify idempotence and whitespace stripping.

## Policy

Tests exercise the public API (`trysac::physics_2d`). They assert observable effects: return values, thrown `PhysicsError` exceptions, and exact code strings. Tests are deterministic and avoid wall-clock time or external network dependency.
