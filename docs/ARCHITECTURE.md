# Architecture

## Purpose

`physics-2d` is a small, dependency-free C++20 library: Simulate 2D rigid bodies deterministically with SAT collision for convex polygons and circles, sequential impulses with restitution and friction, distance joints and sleeping

## Shape

```text
caller code --> include/trysac/physics_2d/*.hpp (namespace trysac::physics_2d)
                    |
                    +-- src/*.cpp: pure functions; no I/O, no global mutable state
                    +-- anonymous-namespace helpers and private headers under src/ (optional)
```

- **Public API**: every declaration of the public headers under `include/`. Code in `src/` and in
  anonymous namespaces may change without a major version.
- **Errors**: invalid input throws a `std::invalid_argument` subclass of the library (such as
  `BlankInputError`) documented with `@throws`. The library never prints or exits.
- **Side effects**: none. If a feature needs I/O, accept a `std::istream&`, a `std::span` or a small
  interface instead of opening resources, and record the choice in an ADR.

## Quality contract

`scripts/check.cmake` is the single definition of the gates; `fabrica.json` tells the factory and
CI how to run it. See [TESTING.md](TESTING.md).

## Decisions

Architecture decisions are recorded in [decisions/](decisions/) as numbered ADRs.
