# physics-2d

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![CI](https://github.com/Trysac/u947-physics-2d/actions/workflows/ci.yml/badge.svg)](https://github.com/Trysac/u947-physics-2d/actions/workflows/ci.yml)

Simulate 2D rigid bodies deterministically with SAT collision for convex polygons and circles, sequential impulses with restitution and friction, distance joints and sleeping

## Problem

Describe in two or three sentences the concrete problem this utility solves, who hits it, and what
they do today without it. Replace this paragraph when the first real feature lands.

## Install

A C++20 library built with CMake. Until it is published to a vcpkg registry, add this repository
at a release tag as a subdirectory (or through `FetchContent`) and link its target:

```cmake
add_subdirectory(u947-physics-2d)
target_link_libraries(app PRIVATE trysac::physics_2d)
```

Requires a C++20 compiler and CMake 3.29 or newer. No dependencies outside the C++ standard
library.

## Usage

```cpp
#include <trysac/physics_2d/text.hpp>

// normalize trims surrounding whitespace and lower-cases the text.
auto name = trysac::physics_2d::normalize("  Hello World ");  // "hello world"

// Blank input is rejected explicitly instead of returning an empty string.
trysac::physics_2d::normalize("   ");  // throws BlankInputError ("input is blank")
```

## Documentation

- [Architecture](docs/ARCHITECTURE.md) and [decisions](docs/decisions/)
- [Testing](docs/TESTING.md) and [security](docs/SECURITY.md)
- [Project record](docs/PROJECT.md) and [stack](docs/STACK.md)
- [Changelog](CHANGELOG.md) and [contributing guide](CONTRIBUTING.md)

## License

[MIT](LICENSE) (c) 2026 Miguel Gomez
