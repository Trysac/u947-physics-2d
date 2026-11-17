# Deterministic 2D Rigid-Body Step Pipeline and SAT Solver

Status: accepted (generated, pending human review)

## Context

Simulations, replays, and lockstep networked games require 2D rigid-body physics that produce identical bit-for-bit state across runs while conserving momentum and elastic collision energy without complex island or tree state.

## Decision

Implement a single-header/single-source C++20 engine in double precision using a fixed 7-stage step pipeline: semi-implicit Euler velocity integration, x-axis sort-and-sweep broadphase tied by BodyId, SAT narrowphase with 1e-9 reference-face tie-breaking, non-warm-started sequential impulses with pre-iteration restitution bias, position integration, non-linear position projection, and neighbor-checked calm-timer sleeping.

## Consequences

Simulations are bitwise deterministic and verifiable via 64-bit FNV-1a state hashing, elastic collisions conserve both momentum and kinetic energy within 1e-9, and resting stacks settle into sleep without velocity drift.
