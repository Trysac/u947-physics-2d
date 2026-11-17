# Project record

## Objective

Simulate 2D rigid bodies deterministically with SAT collision for convex polygons and circles, sequential impulses with restitution and friction, distance joints and sleeping.

## Scope

- In scope: C++20 physics library implementing `trysac::physics_2d` namespace, `World`, `BodyDef`, `Shape`, `Vec2` types, SAT-based `collide` function, and deterministic sequential impulse solver with Sleeping/Waking support.
- Out of scope: continuous collision detection, joints other than distance, polygon radius (rounded polygons), concave shapes, broadphase trees (sort-and-sweep on X axis used), multithreading, and external rendering.

## Effort metrics

| Metric | Value |
|---|---|
| Tokens (input / output) | 12424419 / 74965 |
| Hours (human) | pending measurement |
| Cost (USD) | 2.3533 |

## Factory run u947 (measured from API usage)

| Date (UTC) | Phase | Model | Batch | Input | Output | Cache write | Cache read | Cost (USD) |
|---|---|---|---|---|---|---|---|---|
| 2026-11-17 | generation | google/gemini-4-argon | yes | 6895 | 26622 | 0 | 0 | 0.5600 |
| 2026-11-17 | fix-1 | google/gemini-3.8-flash | no | 191768 | 3563 | 0 | 1957704 | 0.3040 |
| 2026-11-17 | fix-2 | google/gemini-3.8-flash | no | 223747 | 1534 | 0 | 2541750 | 0.3642 |
| 2026-11-17 | fix-3 | google/gemini-3.8-flash | no | 221802 | 28802 | 0 | 3118344 | 0.5082 |
| 2026-11-17 | fix-4 | google/gemini-3.8-flash | no | 150970 | 11772 | 0 | 2460994 | 0.3419 |
| 2026-11-17 | audit | google/gemini-3.8-flash | no | 217403 | 1320 | 0 | 1305315 | 0.2659 |
| 2026-11-17 | docs | google/gemini-3.1-flash-lite | no | 27727 | 1352 | 0 | 0 | 0.0090 |

Total: $2.3533. Prices: vertex line, as measured by its adapter. Providers serve the current version of each model id, so the date tells which version answered; versions, launch stages and price sources per line are in the factory's model register (docs/MODELS.md).
