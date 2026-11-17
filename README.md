# physics-2d

Deterministic 2D rigid-body physics library in C++20 (`trysac::physics_2d`) featuring Separating Axis Theorem (SAT) collision detection for circles and convex polygons, sort-and-sweep broadphase, sequential impulse contact and distance-joint solver, non-linear position projection, and body sleeping.

## Usage

```cpp
#include "trysac/physics_2d/physics.hpp"

int main() {
  trysac::physics_2d::World world;
  trysac::physics_2d::BodyDef def{.type = trysac::physics_2d::BodyType::dynamic_body,
                                  .position = {0.0, 10.0},
                                  .shape = trysac::physics_2d::make_box(0.5, 0.5)};
  const auto body = world.create_body(def);
  world.step(1.0 / 60.0);
  return world.state(body).position.y < 10.0 ? 0 : 1;
}
```

## API Summary

- **Shapes**: `Circle` {radius}, `Polygon` {std::vector<Vec2>}, `make_box(half_width, half_height)`. Validations: 3-8 vertices, strictly convex, non-degenerate.
- **Collisions**: `collide(Shape A, Transform A, Shape B, Transform B)` returns `std::optional<Manifold>` (normal, penetration, points).
- **World**: `create_body(BodyDef)`, `destroy_body(BodyId)`, `add_distance_joint(...)`, `apply_impulse(...)`, `step(double dt)`.
- **Constraints/Dynamics**: Friction `sqrt(f_a * f_b)`, Restitution `max(e_a, e_b)`, distance joints, semi-implicit Euler integration.
- **Limits**: `dt` in (0, 0.1], `max_bodies` (default 10,000), requires finite `double` inputs.
- **Errors**: `PhysicsError` (code: 'invalid_mass', 'too_many_vertices', 'not_convex', 'invalid_step', 'too_many_bodies', 'invalid_value', 'unknown_body', 'invalid_joint').

## Why not Box2D?

Box2D is a mature industry standard but is large and complex; `physics-2d` is a compact, single-namespace C++20 implementation designed for study and precise, bit-level deterministic results across replays.

## Security

See [docs/SECURITY.md](docs/SECURITY.md) for the threat model.
