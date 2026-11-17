# Architecture

`trysac::physics_2d` is a deterministic 2D rigid-body physics engine written in standard C++20 using IEEE-754 double precision.

## Module Layout

- `include/trysac/physics_2d/physics.hpp`: Public API declarations (`Vec2`, `Circle`, `Polygon`, `Shape`, `Transform`, `Manifold`, `collide`, `World`, `PhysicsError`).
- `src/physics.cpp`: Implementation of SAT geometry routines, mass property integration, sort-and-sweep broadphase, sequential impulse velocity and position solvers, sleeping, and FNV-1a state hashing.

## Fixed Simulation Step Order

Every call to `World::step(double dt)` executes the following deterministic sequence:

1. **Velocity Integration**: Awake dynamic bodies advance linear velocity by `gravity * dt`.
2. **Broadphase (Sort-and-Sweep)**: Body AABBs are sorted along `min_x` with ties broken by `BodyId`. Overlapping AABB pairs are emitted in canonical `(min_id, max_id)` order.
3. **Narrowphase (SAT)**:
   - **Polygon–Polygon**: Axis of minimum penetration with a `1e-9` relative tolerance favoring the first body as the reference face, followed by Sutherland–Hodgman segment clipping.
   - **Circle–Polygon**: Closest Voronoi feature (vertex or face) in the polygon's local frame.
   - **Circle–Circle**: Analytical distance test.
4. **Sequential Impulse Velocity Solver**: Warm starting is disabled. Runs `velocity_iterations` passes over contacts (in canonical body-pair order) and then distance joints (in creation order). Restitution bias uses the pre-iteration relative normal velocity whenever `restitution = max(e_a, e_b) > 0`. Tangential impulse is clamped to Coulomb friction `sqrt(f_a * f_b) * normal_impulse`.
5. **Position Integration**: Awake dynamic bodies advance `position += velocity * dt` and `angle += angular_velocity * dt`.
6. **Position Projection**: Runs `position_iterations` passes of Box2D-style non-linear projection (`slop = 0.005`, `factor = 0.2`, `max_correction = 0.2`), modifying positions and angles only and never velocities.
7. **Sleeping**: Without island construction, each dynamic body tracks how long it has stayed below `sleep_linear` and `sleep_angular`. A body transitions to sleep when both it and every body currently touching it have met the calm threshold for at least `sleep_time`.
