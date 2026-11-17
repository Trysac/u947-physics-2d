# Security notes

## Threat model

`physics-2d` is a pure C++ library designed for games and simulations. It performs no I/O operations.

- **Defended**: The library mitigates malformed inputs through explicit `PhysicsError` exceptions. All mathematical operations are checked for finiteness (NaN/Inf) using `std::isfinite`. It enforces strict constraints on polygon geometry (vertex count, convexity, non-degeneracy) and body limits (max bodies, time steps).
- **Not defended**: This is not a hardening library for arbitrary malicious execution of user-supplied scripts or untrusted binary code. It does not implement sandboxing for the host process. Memory usage is proportional to the number of bodies and joints added.
