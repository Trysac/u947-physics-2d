/// @file
/// @brief Deterministic 2D rigid body physics engine with SAT collisions, sequential impulses, and sleeping.
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

/// @brief The public API of the trysac::physics_2d library.
namespace trysac::physics_2d {

/// @brief Base exception class for physics errors with a stable machine-readable code.
class PhysicsError : public std::invalid_argument {
 public:
  /// @brief Constructs a PhysicsError with the given machine-readable error code and message.
  /// @param error_code The stable machine-readable error identifier.
  /// @param message Human-readable error description.
  PhysicsError(std::string_view error_code, const std::string& message);

  /// @brief Returns the stable machine-readable error code.
  /// @return The error code view.
  [[nodiscard]] std::string_view code() const noexcept;

 private:
  std::string code_;
};

/// @brief A 2D vector with double precision.
struct Vec2 {
  /// @brief The X coordinate.
  double x = 0.0;
  /// @brief The Y coordinate.
  double y = 0.0;

  /// @brief Default-constructs a zero vector.
  constexpr Vec2() noexcept = default;

  /// @brief Constructs a vector with coordinates (x_val, y_val).
  /// @param x_val The X component.
  /// @param y_val The Y component.
  constexpr Vec2(double x_val, double y_val) noexcept : x(x_val), y(y_val) {}

  /// @brief Unary negation.
  /// @return The negated vector.
  [[nodiscard]] constexpr Vec2 operator-() const noexcept { return {-x, -y}; }

  /// @brief Vector addition.
  /// @param rhs The vector to add.
  /// @return The sum vector.
  [[nodiscard]] constexpr Vec2 operator+(const Vec2& rhs) const noexcept { return {x + rhs.x, y + rhs.y}; }

  /// @brief Vector subtraction.
  /// @param rhs The vector to subtract.
  /// @return The difference vector.
  [[nodiscard]] constexpr Vec2 operator-(const Vec2& rhs) const noexcept { return {x - rhs.x, y - rhs.y}; }

  /// @brief Scalar multiplication.
  /// @param s The scalar factor.
  /// @return The scaled vector.
  [[nodiscard]] constexpr Vec2 operator*(double s) const noexcept { return {x * s, y * s}; }

  /// @brief Scalar division.
  /// @param s The scalar divisor.
  /// @return The divided vector.
  [[nodiscard]] constexpr Vec2 operator/(double s) const noexcept { return {x / s, y / s}; }

  /// @brief In-place vector addition.
  /// @param rhs The vector to add.
  /// @return Reference to this vector.
  constexpr Vec2& operator+=(const Vec2& rhs) noexcept {
    x += rhs.x;
    y += rhs.y;
    return *this;
  }

  /// @brief In-place vector subtraction.
  /// @param rhs The vector to subtract.
  /// @return Reference to this vector.
  constexpr Vec2& operator-=(const Vec2& rhs) noexcept {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
  }

  /// @brief In-place scalar multiplication.
  /// @param s The scalar factor.
  /// @return Reference to this vector.
  constexpr Vec2& operator*=(double s) noexcept {
    x *= s;
    y *= s;
    return *this;
  }

  /// @brief Equality comparison.
  /// @param rhs The vector to compare with.
  /// @return True if both components are equal.
  [[nodiscard]] constexpr bool operator==(const Vec2& rhs) const noexcept = default;

  /// @brief Returns the squared Euclidean length of the vector.
  /// @return Squared length.
  [[nodiscard]] constexpr double length_squared() const noexcept { return (x * x) + (y * y); }

  /// @brief Returns the Euclidean length of the vector.
  /// @return Length.
  [[nodiscard]] double length() const noexcept;

  /// @brief Returns a unit-length vector in the same direction, or (0, 0) if near zero.
  /// @return Normalized vector.
  [[nodiscard]] Vec2 normalized() const noexcept;
};

/// @brief Left-hand scalar multiplication.
/// @param s The scalar factor.
/// @param v The vector.
/// @return The scaled vector.
[[nodiscard]] constexpr Vec2 operator*(double s, const Vec2& v) noexcept { return {s * v.x, s * v.y}; }

/// @brief Computes the 2D dot product of two vectors.
/// @param a First vector.
/// @param b Second vector.
/// @return Dot product a.x * b.x + a.y * b.y.
[[nodiscard]] constexpr double dot(const Vec2& a, const Vec2& b) noexcept { return (a.x * b.x) + (a.y * b.y); }

/// @brief Computes the 2D scalar cross product (z-component of 3D cross product).
/// @param a First vector.
/// @param b Second vector.
/// @return Cross product a.x * b.y - a.y * b.x.
[[nodiscard]] constexpr double cross(const Vec2& a, const Vec2& b) noexcept { return (a.x * b.y) - (a.y * b.x); }

/// @brief Computes the cross product of a vector with a scalar (z-axis vector).
/// @param a The 2D vector.
/// @param s The scalar z-component.
/// @return Vector (s * a.y, -s * a.x).
[[nodiscard]] constexpr Vec2 cross(const Vec2& a, double s) noexcept { return {s * a.y, -s * a.x}; }

/// @brief Computes the cross product of a scalar (z-axis vector) with a 2D vector.
/// @param s The scalar z-component.
/// @param a The 2D vector.
/// @return Vector (-s * a.y, s * a.x).
[[nodiscard]] constexpr Vec2 cross(double s, const Vec2& a) noexcept { return {-s * a.y, s * a.x}; }

/// @brief A circle shape centered at the body's local origin.
struct Circle {
  /// @brief Radius of the circle (must be positive and finite).
  double radius = 1.0;
};

/// @brief A convex polygon shape with 3 to 8 vertices in local coordinates.
struct Polygon {
  /// @brief Vertices of the polygon (normalized to counter-clockwise order on validation).
  std::vector<Vec2> vertices;

  /// @brief Default-constructs an empty polygon.
  Polygon() = default;

  /// @brief Constructs and validates a convex polygon from vertices.
  /// @param verts Between 3 and 8 vertices; clockwise input is reversed to counter-clockwise.
  /// @throws PhysicsError with code `not_convex`, `too_many_vertices`, or `invalid_value`.
  explicit Polygon(std::vector<Vec2> verts);
};

/// @brief A geometric collision shape: either a Circle or a convex Polygon.
using Shape = std::variant<Circle, Polygon>;

/// @brief Creates an axis-aligned box polygon centered at the local origin.
/// @param half_width Half the width of the box (must be > 0 and finite).
/// @param half_height Half the height of the box (must be > 0 and finite).
/// @return A validated counter-clockwise 4-vertex Polygon.
/// @throws PhysicsError with code `invalid_value` if dimensions are not positive finite.
[[nodiscard]] Polygon make_box(double half_width, double half_height);

/// @brief A 2D rigid transform (translation and rotation).
struct Transform {
  /// @brief World position of the local origin.
  Vec2 position{0.0, 0.0};
  /// @brief Rotation angle in radians counter-clockwise.
  double angle = 0.0;
};

/// @brief Contact manifold describing the overlap between two shapes.
struct Manifold {
  /// @brief Unit contact normal pointing from shape A to shape B.
  Vec2 normal{0.0, 0.0};
  /// @brief Overlap depth along the contact normal.
  double penetration = 0.0;
  /// @brief World-space contact points (1 or 2 points).
  std::vector<Vec2> points;
};

/// @brief Computes the collision manifold between two shapes under given transforms.
/// @param a The first shape.
/// @param ta World transform of the first shape.
/// @param b The second shape.
/// @param tb World transform of the second shape.
/// @return The collision Manifold if the shapes overlap, or std::nullopt if separated.
/// @throws PhysicsError if any shape or transform is invalid.
[[nodiscard]] std::optional<Manifold> collide(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb);

/// @brief Classification of a rigid body's motion behavior.
enum class BodyType : std::uint8_t {
  /// @brief Zero velocity, infinite mass, immovable by forces or collisions.
  static_body,
  /// @brief Positive mass, simulated under forces, impulses, and constraints.
  dynamic_body,
};

/// @brief Definition parameters used to create a body in a World.
struct BodyDef {
  /// @brief Static or dynamic body type.
  BodyType type = BodyType::dynamic_body;
  /// @brief Initial world position.
  Vec2 position{0.0, 0.0};
  /// @brief Initial rotation angle in radians.
  double angle = 0.0;
  /// @brief Initial linear velocity.
  Vec2 velocity{0.0, 0.0};
  /// @brief Initial angular velocity in radians/second.
  double angular_velocity = 0.0;
  /// @brief Collision shape of the body.
  Shape shape = Circle{1.0};
  /// @brief Area mass density (must be > 0 for dynamic bodies).
  double density = 1.0;
  /// @brief Coefficient of restitution in [0, 1].
  double restitution = 0.0;
  /// @brief Coefficient of Coulomb friction (>= 0).
  double friction = 0.3;
};

/// @brief Simulation configuration settings for a World.
struct WorldSettings {
  /// @brief Gravitational acceleration applied to awake dynamic bodies.
  Vec2 gravity{0.0, -10.0};
  /// @brief Number of sequential impulse velocity solver iterations per step.
  int velocity_iterations = 10;
  /// @brief Number of non-linear position projection iterations per step.
  int position_iterations = 4;
  /// @brief Time in seconds a body and its contacts must remain calm before sleeping.
  double sleep_time = 0.5;
  /// @brief Linear speed threshold for sleeping.
  double sleep_linear = 0.01;
  /// @brief Angular speed threshold for sleeping.
  double sleep_angular = 0.01;
  /// @brief Maximum number of bodies allowed simultaneously in the world.
  std::size_t max_bodies = 10'000;
};

/// @brief Identifier for a rigid body inside a World.
using BodyId = std::uint64_t;

/// @brief Identifier for a constraint joint inside a World.
using JointId = std::uint64_t;

/// @brief Snapshot of the kinematic state of a body.
struct BodyState {
  /// @brief Current world position.
  Vec2 position{0.0, 0.0};
  /// @brief Current rotation angle in radians.
  double angle = 0.0;
  /// @brief Current linear velocity.
  Vec2 velocity{0.0, 0.0};
  /// @brief Current angular velocity in radians/second.
  double angular_velocity = 0.0;
};

/// @brief Information about an active contact between two bodies.
struct ContactInfo {
  /// @brief Identifier of the first body.
  BodyId body_a = 0;
  /// @brief Identifier of the second body.
  BodyId body_b = 0;
  /// @brief Contact manifold in world coordinates.
  Manifold manifold;
};

/// @brief Deterministic 2D physics simulation world.
class World {
 public:
  /// @brief Constructs a physics world with the given settings.
  /// @param s World configuration settings.
  /// @throws PhysicsError with code `invalid_value` if settings contain invalid values.
  explicit World(WorldSettings s = {});

  /// @brief Creates a rigid body from a definition.
  /// @param def The body definition.
  /// @return The unique BodyId assigned to the new body.
  /// @throws PhysicsError on invalid values, invalid mass, invalid shape, or exceeding max_bodies.
  BodyId create_body(const BodyDef& def);

  /// @brief Destroys an existing body and removes any attached joints and contacts.
  /// @param id The identifier of the body to destroy.
  /// @throws PhysicsError with code `unknown_body` if the body does not exist.
  void destroy_body(BodyId id);

  /// @brief Adds a distance constraint joint between two bodies.
  /// @param a First body identifier.
  /// @param b Second body identifier.
  /// @param local_anchor_a Anchor point in body A's local coordinates.
  /// @param local_anchor_b Anchor point in body B's local coordinates.
  /// @param length Target rest distance between the world anchors (>= 0).
  /// @return The unique JointId assigned to the new joint.
  /// @throws PhysicsError with code `unknown_body`, `invalid_joint`, or `invalid_value`.
  JointId add_distance_joint(BodyId a, BodyId b, Vec2 local_anchor_a, Vec2 local_anchor_b, double length);

  /// @brief Advances the physics simulation by time step `dt`.
  /// @param dt Time step in seconds, must be in (0, 0.1].
  /// @throws PhysicsError with code `invalid_step` or `invalid_value` if `dt` is out of bounds or non-finite.
  void step(double dt);

  /// @brief Applies an instantaneous linear and angular impulse to a dynamic body at a world point.
  ///        Wakes the body and any sleeping bodies currently touching it.
  /// @param id The target body identifier.
  /// @param impulse The impulse vector in N*s.
  /// @param world_point The application point in world coordinates.
  /// @throws PhysicsError with code `unknown_body` or `invalid_value`.
  void apply_impulse(BodyId id, Vec2 impulse, Vec2 world_point);

  /// @brief Queries the current kinematic state of a body.
  /// @param id The body identifier.
  /// @return The body's position, angle, velocity, and angular velocity.
  /// @throws PhysicsError with code `unknown_body` if the body does not exist.
  [[nodiscard]] BodyState state(BodyId id) const;

  /// @brief Checks whether a body is currently awake.
  /// @param id The body identifier.
  /// @return True if awake, false if sleeping (or static).
  /// @throws PhysicsError with code `unknown_body` if the body does not exist.
  [[nodiscard]] bool is_awake(BodyId id) const;

  /// @brief Returns the active contact manifolds from the most recent step.
  /// @return Vector of ContactInfo ordered deterministically.
  [[nodiscard]] std::vector<ContactInfo> contacts() const;

  /// @brief Computes a 64-bit FNV-1a hash of the IEEE-754 bit patterns of all bodies' kinematic states.
  /// @return The deterministic 64-bit state hash.
  [[nodiscard]] std::uint64_t state_hash() const;

 private:
  struct BodyRecord {
    BodyId id = 0;
    BodyType type = BodyType::dynamic_body;
    Vec2 position{0.0, 0.0};
    double angle = 0.0;
    Vec2 velocity{0.0, 0.0};
    double angular_velocity = 0.0;
    Shape shape = Circle{1.0};
    Vec2 local_centroid{0.0, 0.0};
    double mass = 0.0;
    double inv_mass = 0.0;
    double inertia = 0.0;
    double inv_inertia = 0.0;
    double restitution = 0.0;
    double friction = 0.3;
    bool awake = true;
    double calm_time = 0.0;
  };

  struct DistanceJointRecord {
    JointId id = 0;
    BodyId body_a = 0;
    BodyId body_b = 0;
    Vec2 local_anchor_a{0.0, 0.0};
    Vec2 local_anchor_b{0.0, 0.0};
    double length = 0.0;
  };

  struct ContactPointSolver {
    Vec2 point{0.0, 0.0};
    Vec2 local_a{0.0, 0.0};
    Vec2 local_b{0.0, 0.0};
    Vec2 r_a{0.0, 0.0};
    Vec2 r_b{0.0, 0.0};
    double normal_mass = 0.0;
    double tangent_mass = 0.0;
    double normal_impulse = 0.0;
    double tangent_impulse = 0.0;
    double velocity_bias = 0.0;
  };

  struct ContactConstraint {
    std::size_t index_a = 0;
    std::size_t index_b = 0;
    BodyId id_a = 0;
    BodyId id_b = 0;
    Vec2 normal{0.0, 0.0};
    Vec2 tangent{0.0, 0.0};
    double penetration = 0.0;
    double friction = 0.0;
    double restitution = 0.0;
    std::vector<ContactPointSolver> points;
    Manifold raw_manifold;
  };

  [[nodiscard]] std::size_t find_body_index(BodyId id) const;
  void wake_connected(std::size_t start_idx);

  WorldSettings settings_;
  std::vector<BodyRecord> bodies_;
  std::vector<DistanceJointRecord> joints_;
  std::vector<ContactInfo> last_contacts_;
  BodyId next_body_id_ = 1;
  JointId next_joint_id_ = 1;
};

}  // namespace trysac::physics_2d
