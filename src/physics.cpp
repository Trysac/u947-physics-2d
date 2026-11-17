#include "trysac/physics_2d/physics.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace trysac::physics_2d {

namespace {

constexpr std::size_t kMinPolygonVertices = 3;
constexpr std::size_t kMaxPolygonVertices = 8;
constexpr double kMaxTimeStep = 0.1;
constexpr double kLinearSlop = 0.005;
constexpr double kBaumgarteFactor = 0.2;
constexpr double kMaxLinearCorrection = 0.2;
constexpr double kSatRelativeTol = 1e-9;
constexpr double kCollinearEps = 1e-12;

bool is_finite(double v) noexcept { return std::isfinite(v); }

bool is_finite(const Vec2& v) noexcept { return std::isfinite(v.x) && std::isfinite(v.y); }

void require_finite(double v, const char* msg) {
  if (!is_finite(v)) {
    throw PhysicsError("invalid_value", msg);
  }
}

void require_finite(const Vec2& v, const char* msg) {
  if (!is_finite(v)) {
    throw PhysicsError("invalid_value", msg);
  }
}

struct Rotation {
  double c = 1.0;
  double s = 0.0;

  explicit Rotation(double angle) noexcept : c(std::cos(angle)), s(std::sin(angle)) {}
};

Vec2 rotate(const Rotation& q, const Vec2& v) noexcept {
  return {(q.c * v.x) - (q.s * v.y), (q.s * v.x) + (q.c * v.y)};
}

Vec2 inv_rotate(const Rotation& q, const Vec2& v) noexcept {
  return {(q.c * v.x) + (q.s * v.y), (-q.s * v.x) + (q.c * v.y)};
}

Vec2 transform_point(Vec2 pos, const Rotation& q, const Vec2& local) noexcept { return pos + rotate(q, local); }

Vec2 inv_transform_point(Vec2 pos, const Rotation& q, const Vec2& world) noexcept { return inv_rotate(q, world - pos); }

std::vector<Vec2> validate_and_orient_polygon(std::vector<Vec2> verts) {
  if (verts.size() > kMaxPolygonVertices) {
    throw PhysicsError("too_many_vertices", "polygon has more than 8 vertices");
  }
  if (verts.size() < kMinPolygonVertices) {
    throw PhysicsError("not_convex", "polygon must have at least 3 vertices");
  }
  for (const auto& v : verts) {
    require_finite(v, "polygon vertex must be finite");
  }

  double signed_area2 = 0.0;
  const std::size_t n = verts.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& p0 = verts.at(i);
    const Vec2& p1 = verts.at((i + 1) % n);
    signed_area2 += cross(p0, p1);
  }
  if (std::abs(signed_area2) <= kCollinearEps) {
    throw PhysicsError("not_convex", "polygon has zero area or degenerate vertices");
  }
  if (signed_area2 < 0.0) {
    std::ranges::reverse(verts);
  }

  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = verts.at(i);
    const Vec2& b = verts.at((i + 1) % n);
    const Vec2& c = verts.at((i + 2) % n);
    const Vec2 ab = b - a;
    const Vec2 bc = c - b;
    const double ab_len2 = ab.length_squared();
    const double bc_len2 = bc.length_squared();
    if (ab_len2 <= kCollinearEps || bc_len2 <= kCollinearEps) {
      throw PhysicsError("not_convex", "polygon has duplicate vertices");
    }
    const double cr = cross(ab, bc);
    if (cr <= kCollinearEps * std::sqrt(ab_len2 * bc_len2)) {
      throw PhysicsError("not_convex", "polygon is not strictly convex");
    }
  }

  // Check that every vertex is strictly on the interior side of every non-adjacent edge
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = verts.at(i);
    const Vec2& b = verts.at((i + 1) % n);
    const Vec2 edge = b - a;
    for (std::size_t j = 0; j < n; ++j) {
      if (j == i || j == (i + 1) % n) {
        continue;
      }
      if (cross(edge, verts.at(j) - a) <= kCollinearEps) {
        throw PhysicsError("not_convex", "polygon is self-intersecting or non-convex");
      }
    }
  }

  return verts;
}

struct MassData {
  double mass = 0.0;
  double inertia = 0.0;
  Vec2 centroid{0.0, 0.0};
};

MassData compute_shape_mass(const Shape& shape, double density) {
  if (const auto* circle = std::get_if<Circle>(&shape)) {
    require_finite(circle->radius, "circle radius must be finite");
    if (circle->radius <= 0.0) {
      throw PhysicsError("invalid_value", "circle radius must be positive");
    }
    const double r2 = circle->radius * circle->radius;
    const double area = std::numbers::pi * r2;
    const double mass = density * area;
    const double inertia = 0.5 * mass * r2;
    return {.mass = mass, .inertia = inertia, .centroid = {0.0, 0.0}};
  }

  const auto& poly = std::get<Polygon>(shape);
  const auto verts = validate_and_orient_polygon(poly.vertices);
  const std::size_t n = verts.size();

  double area = 0.0;
  Vec2 centroid{0.0, 0.0};
  double second_moment_origin = 0.0;

  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& p0 = verts.at(i);
    const Vec2& p1 = verts.at((i + 1) % n);
    const double d = cross(p0, p1);
    const double tri_area = 0.5 * d;
    area += tri_area;
    centroid += (p0 + p1) * (tri_area / 3.0);
    const double ix = (p0.x * p0.x) + (p0.x * p1.x) + (p1.x * p1.x);
    const double iy = (p0.y * p0.y) + (p0.y * p1.y) + (p1.y * p1.y);
    second_moment_origin += (d / 12.0) * (ix + iy);
  }

  centroid = centroid / area;
  const double mass = density * area;
  const double inertia_origin = density * second_moment_origin;
  // Shift inertia to centroid if needed, though local origin is body position
  const double inertia_cm = inertia_origin - (mass * centroid.length_squared());
  return {.mass = mass, .inertia = std::max(inertia_cm, 1e-12 * mass), .centroid = centroid};
}

struct Aabb {
  double min_x = 0.0;
  double min_y = 0.0;
  double max_x = 0.0;
  double max_y = 0.0;
};

Aabb compute_aabb(const Shape& shape, const Transform& t) {
  const Rotation q(t.angle);
  if (const auto* circle = std::get_if<Circle>(&shape)) {
    return {
        .min_x = t.position.x - circle->radius,
        .min_y = t.position.y - circle->radius,
        .max_x = t.position.x + circle->radius,
        .max_y = t.position.y + circle->radius,
    };
  }
  const auto& poly = std::get<Polygon>(shape);
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  for (const auto& v : poly.vertices) {
    const Vec2 w = transform_point(t.position, q, v);
    min_x = std::min(min_x, w.x);
    min_y = std::min(min_y, w.y);
    max_x = std::max(max_x, w.x);
    max_y = std::max(max_y, w.y);
  }
  return {.min_x = min_x, .min_y = min_y, .max_x = max_x, .max_y = max_y};
}

struct WorldPolygon {
  std::vector<Vec2> verts;
  std::vector<Vec2> normals;
};

WorldPolygon build_world_polygon(const Polygon& poly, const Transform& t) {
  const auto ccw_verts = validate_and_orient_polygon(poly.vertices);
  const Rotation q(t.angle);
  const std::size_t n = ccw_verts.size();
  WorldPolygon wp;
  wp.verts.resize(n);
  wp.normals.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    wp.verts.at(i) = transform_point(t.position, q, ccw_verts.at(i));
  }
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2 edge = wp.verts.at((i + 1) % n) - wp.verts.at(i);
    // Outward normal for CCW polygon is (edge.y, -edge.x)
    wp.normals.at(i) = Vec2{edge.y, -edge.x}.normalized();
  }
  return wp;
}

struct EdgeSeparation {
  double max_separation = -std::numeric_limits<double>::infinity();
  std::size_t best_edge = 0;
};

EdgeSeparation find_max_separation(const WorldPolygon& poly_a, const WorldPolygon& poly_b) {
  EdgeSeparation result;
  for (std::size_t i = 0; i < poly_a.verts.size(); ++i) {
    const Vec2& n = poly_a.normals.at(i);
    const Vec2& v = poly_a.verts.at(i);
    double deepest = std::numeric_limits<double>::infinity();
    for (const auto& vb : poly_b.verts) {
      const double sep = dot(n, vb - v);
      deepest = std::min(deepest, sep);
    }
    if (deepest > result.max_separation) {
      result.max_separation = deepest;
      result.best_edge = i;
    }
  }
  return result;
}

std::vector<Vec2> clip_segment_to_line(const std::vector<Vec2>& in_pts, const Vec2& normal, double offset) {
  std::vector<Vec2> out;
  if (in_pts.size() < 2) {
    return in_pts;
  }
  const double d0 = dot(normal, in_pts.at(0)) - offset;
  const double d1 = dot(normal, in_pts.at(1)) - offset;
  if (d0 <= 0.0) {
    out.push_back(in_pts.at(0));
  }
  if (d1 <= 0.0) {
    out.push_back(in_pts.at(1));
  }
  if (d0 * d1 < 0.0) {
    const double alpha = d0 / (d0 - d1);
    out.push_back(in_pts.at(0) + (in_pts.at(1) - in_pts.at(0)) * alpha);
  }
  return out;
}

std::optional<Manifold> collide_circles(const Circle& a, const Transform& ta, const Circle& b, const Transform& tb) {
  require_finite(a.radius, "circle radius must be finite");
  require_finite(b.radius, "circle radius must be finite");
  if (a.radius <= 0.0 || b.radius <= 0.0) {
    throw PhysicsError("invalid_value", "circle radius must be positive");
  }
  const Vec2 d = tb.position - ta.position;
  const double dist2 = d.length_squared();
  const double r_sum = a.radius + b.radius;
  if (dist2 > r_sum * r_sum) {
    return std::nullopt;
  }
  const double dist = std::sqrt(dist2);
  Manifold m;
  if (dist <= 1e-12) {
    m.normal = {1.0, 0.0};
    m.penetration = r_sum;
    m.points.push_back(ta.position);
  } else {
    m.normal = d / dist;
    m.penetration = r_sum - dist;
    m.points.push_back(ta.position + (m.normal * (a.radius - (0.5 * m.penetration))));
  }
  return m;
}

std::optional<Manifold> collide_polygons(const Polygon& a, const Transform& ta, const Polygon& b, const Transform& tb) {
  const WorldPolygon wa = build_world_polygon(a, ta);
  const WorldPolygon wb = build_world_polygon(b, tb);

  const EdgeSeparation sep_a = find_max_separation(wa, wb);
  if (sep_a.max_separation > 0.0) {
    return std::nullopt;
  }
  const EdgeSeparation sep_b = find_max_separation(wb, wa);
  if (sep_b.max_separation > 0.0) {
    return std::nullopt;
  }

  const double pen_a = -sep_a.max_separation;
  const double pen_b = -sep_b.max_separation;

  // Relative tolerance 1e-9 favouring the first body
  const bool flip = (pen_b < pen_a * (1.0 - kSatRelativeTol));
  const WorldPolygon& ref_poly = flip ? wb : wa;
  const WorldPolygon& inc_poly = flip ? wa : wb;
  const std::size_t ref_edge = flip ? sep_b.best_edge : sep_a.best_edge;

  const Vec2 ref_v1 = ref_poly.verts.at(ref_edge);
  const Vec2 ref_v2 = ref_poly.verts.at((ref_edge + 1) % ref_poly.verts.size());
  const Vec2 ref_normal = ref_poly.normals.at(ref_edge);
  const Vec2 side_tangent = (ref_v2 - ref_v1).normalized();

  // Find incident edge on inc_poly most anti-parallel to ref_normal
  std::size_t inc_edge = 0;
  double min_dot = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < inc_poly.normals.size(); ++i) {
    const double d = dot(ref_normal, inc_poly.normals.at(i));
    if (d < min_dot) {
      min_dot = d;
      inc_edge = i;
    }
  }

  std::vector<Vec2> clip_pts = {
      inc_poly.verts.at(inc_edge),
      inc_poly.verts.at((inc_edge + 1) % inc_poly.verts.size()),
  };

  const double side_offset1 = -dot(side_tangent, ref_v1);
  const double side_offset2 = dot(side_tangent, ref_v2);

  clip_pts = clip_segment_to_line(clip_pts, -side_tangent, side_offset1);
  if (clip_pts.size() < 2) {
    return std::nullopt;
  }
  clip_pts = clip_segment_to_line(clip_pts, side_tangent, side_offset2);
  if (clip_pts.empty()) {
    return std::nullopt;
  }

  const double front_offset = dot(ref_normal, ref_v1);
  std::vector<Vec2> contact_points;
  double max_pen = 0.0;
  for (const auto& pt : clip_pts) {
    const double separation = dot(ref_normal, pt) - front_offset;
    if (separation <= 1e-12) {
      const double pen = -separation;
      max_pen = std::max(max_pen, pen);
      contact_points.push_back(pt);
    }
  }

  if (contact_points.empty()) {
    return std::nullopt;
  }

  Manifold m;
  m.normal = flip ? -ref_normal : ref_normal;
  m.penetration = flip ? pen_b : pen_a;
  m.points = std::move(contact_points);
  return m;
}

std::optional<Manifold> collide_polygon_circle(const Polygon& poly_in, const Transform& tp, const Circle& circle,
                                               const Transform& tc) {
  require_finite(circle.radius, "circle radius must be finite");
  if (circle.radius <= 0.0) {
    throw PhysicsError("invalid_value", "circle radius must be positive");
  }
  const auto verts = validate_and_orient_polygon(poly_in.vertices);
  const std::size_t n = verts.size();
  if (n < kMinPolygonVertices) {
    return std::nullopt;
  }
  const Rotation qp(tp.angle);
  const Vec2 c_local = inv_transform_point(tp.position, qp, tc.position);

  double max_sep = -std::numeric_limits<double>::infinity();
  std::size_t best_edge = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2 v1 = verts.at(i);
    const Vec2 v2 = verts.at((i + 1) % n);
    const Vec2 edge = v2 - v1;
    const Vec2 normal = Vec2{edge.y, -edge.x}.normalized();
    const double s = dot(normal, c_local - v1);
    if (s > circle.radius) {
      return std::nullopt;
    }
    if (s > max_sep) {
      max_sep = s;
      best_edge = i;
    }
  }

  const Vec2 v1 = verts.at(best_edge);
  const Vec2 v2 = verts.at((best_edge + 1) % n);
  const Vec2 edge = v2 - v1;
  const Vec2 face_normal = Vec2{edge.y, -edge.x}.normalized();

  // Center is inside polygon
  if (max_sep <= 1e-12) {
    Manifold m;
    m.normal = rotate(qp, face_normal);
    m.penetration = circle.radius - max_sep;
    m.points.push_back(tc.position - m.normal * circle.radius);
    return m;
  }

  const double u1 = dot(c_local - v1, v2 - v1);
  const double u2 = dot(c_local - v2, v1 - v2);

  if (u1 <= 0.0) {
    const Vec2 d = c_local - v1;
    const double dist2 = d.length_squared();
    if (dist2 > circle.radius * circle.radius) {
      return std::nullopt;
    }
    const double dist = std::sqrt(dist2);
    Manifold m;
    const Vec2 local_n = (dist > 1e-12) ? (d / dist) : face_normal;
    m.normal = rotate(qp, local_n);
    m.penetration = circle.radius - dist;
    m.points.push_back(transform_point(tp.position, qp, v1));
    return m;
  }

  if (u2 <= 0.0) {
    const Vec2 d = c_local - v2;
    const double dist2 = d.length_squared();
    if (dist2 > circle.radius * circle.radius) {
      return std::nullopt;
    }
    const double dist = std::sqrt(dist2);
    Manifold m;
    const Vec2 local_n = (dist > 1e-12) ? (d / dist) : face_normal;
    m.normal = rotate(qp, local_n);
    m.penetration = circle.radius - dist;
    m.points.push_back(transform_point(tp.position, qp, v2));
    return m;
  }

  // Face region
  const double face_dist = dot(c_local - v1, face_normal);
  if (face_dist > circle.radius) {
    return std::nullopt;
  }
  Manifold m;
  m.normal = rotate(qp, face_normal);
  m.penetration = circle.radius - face_dist;
  m.points.push_back(tc.position - m.normal * circle.radius);
  return m;
}

void fnv1a_append_double(std::uint64_t& hash, double value) noexcept {
  // Canonicalize -0.0 to +0.0 for deterministic IEEE-754 hashing
  const double canonical = (value == 0.0) ? 0.0 : value;
  const auto bits = std::bit_cast<std::uint64_t>(canonical);
  for (unsigned int i = 0; i < 8; ++i) {
    const auto byte = static_cast<std::uint8_t>((bits >> (i * 8U)) & 0xFFU);
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
}

}  // namespace

PhysicsError::PhysicsError(std::string_view error_code, const std::string& message)
    : std::invalid_argument(message), code_(error_code) {}

std::string_view PhysicsError::code() const noexcept { return code_; }

double Vec2::length() const noexcept { return std::sqrt(length_squared()); }

Vec2 Vec2::normalized() const noexcept {
  const double len = length();
  if (len <= 1e-15) {
    return {0.0, 0.0};
  }
  return {x / len, y / len};
}

Polygon::Polygon(std::vector<Vec2> verts) : vertices(validate_and_orient_polygon(std::move(verts))) {}

Polygon make_box(double half_width, double half_height) {
  require_finite(half_width, "box half_width must be finite");
  require_finite(half_height, "box half_height must be finite");
  if (half_width <= 0.0 || half_height <= 0.0) {
    throw PhysicsError("invalid_value", "box half extents must be positive");
  }
  return Polygon({
      {-half_width, -half_height},
      {half_width, -half_height},
      {half_width, half_height},
      {-half_width, half_height},
  });
}

std::optional<Manifold> collide(const Shape& a, const Transform& ta, const Shape& b, const Transform& tb) {
  require_finite(ta.position, "transform A position must be finite");
  require_finite(ta.angle, "transform A angle must be finite");
  require_finite(tb.position, "transform B position must be finite");
  require_finite(tb.angle, "transform B angle must be finite");

  if (const auto* ca = std::get_if<Circle>(&a)) {
    if (const auto* cb = std::get_if<Circle>(&b)) {
      return collide_circles(*ca, ta, *cb, tb);
    }
    const auto& pb = std::get<Polygon>(b);
    auto m = collide_polygon_circle(pb, tb, *ca, ta);
    if (m.has_value()) {
      m->normal = -m->normal;
    }
    return m;
  }

  const auto& pa = std::get<Polygon>(a);
  if (const auto* cb = std::get_if<Circle>(&b)) {
    return collide_polygon_circle(pa, ta, *cb, tb);
  }
  const auto& pb = std::get<Polygon>(b);
  return collide_polygons(pa, ta, pb, tb);
}

World::World(WorldSettings s) : settings_(s) {
  require_finite(settings_.gravity, "gravity must be finite");
  require_finite(settings_.sleep_time, "sleep_time must be finite");
  require_finite(settings_.sleep_linear, "sleep_linear must be finite");
  require_finite(settings_.sleep_angular, "sleep_angular must be finite");
  if (settings_.velocity_iterations < 1 || settings_.position_iterations < 0 || settings_.sleep_time < 0.0 ||
      settings_.sleep_linear < 0.0 || settings_.sleep_angular < 0.0 || settings_.max_bodies == 0) {
    throw PhysicsError("invalid_value", "invalid world settings");
  }
}

BodyId World::create_body(const BodyDef& def) {
  if (bodies_.size() >= settings_.max_bodies) {
    throw PhysicsError("too_many_bodies", "maximum number of bodies reached");
  }
  require_finite(def.position, "body position must be finite");
  require_finite(def.angle, "body angle must be finite");
  require_finite(def.velocity, "body velocity must be finite");
  require_finite(def.angular_velocity, "body angular_velocity must be finite");
  require_finite(def.density, "body density must be finite");
  require_finite(def.restitution, "body restitution must be finite");
  require_finite(def.friction, "body friction must be finite");

  if (def.restitution < 0.0 || def.restitution > 1.0 || def.friction < 0.0) {
    throw PhysicsError("invalid_value", "body restitution must be in [0, 1] and friction >= 0");
  }

  if (def.type == BodyType::dynamic_body && def.density <= 0.0) {
    throw PhysicsError("invalid_mass", "dynamic body must have positive density");
  }

  Shape validated_shape = def.shape;
  if (auto* poly = std::get_if<Polygon>(&validated_shape)) {
    poly->vertices = validate_and_orient_polygon(poly->vertices);
  }

  const double effective_density = (def.type == BodyType::dynamic_body) ? def.density : 1.0;
  const MassData md = compute_shape_mass(validated_shape, effective_density);
  if (def.type == BodyType::dynamic_body && (!is_finite(md.mass) || md.mass <= 0.0 || !is_finite(md.inertia))) {
    throw PhysicsError("invalid_mass", "computed mass or inertia is invalid");
  }

  BodyRecord rec;
  rec.id = next_body_id_++;
  rec.type = def.type;
  rec.position = def.position;
  rec.angle = def.angle;
  rec.velocity = (def.type == BodyType::dynamic_body) ? def.velocity : Vec2{0.0, 0.0};
  rec.angular_velocity = (def.type == BodyType::dynamic_body) ? def.angular_velocity : 0.0;
  rec.shape = std::move(validated_shape);
  rec.local_centroid = md.centroid;
  rec.restitution = def.restitution;
  rec.friction = def.friction;
  rec.awake = (def.type == BodyType::dynamic_body);
  rec.calm_time = 0.0;

  if (def.type == BodyType::dynamic_body) {
    rec.mass = md.mass;
    rec.inv_mass = 1.0 / md.mass;
    rec.inertia = md.inertia;
    rec.inv_inertia = 1.0 / md.inertia;
  }

  bodies_.push_back(std::move(rec));
  return bodies_.back().id;
}

void World::destroy_body(BodyId id) {
  const std::size_t idx = find_body_index(id);
  bodies_.erase(bodies_.begin() + static_cast<std::ptrdiff_t>(idx));
  std::erase_if(joints_, [id](const DistanceJointRecord& j) { return j.body_a == id || j.body_b == id; });
  std::erase_if(last_contacts_, [id](const ContactInfo& c) { return c.body_a == id || c.body_b == id; });
}

JointId World::add_distance_joint(BodyId a, BodyId b, Vec2 local_anchor_a, Vec2 local_anchor_b, double length) {
  const std::size_t idx_a = find_body_index(a);
  const std::size_t idx_b = find_body_index(b);
  require_finite(local_anchor_a, "local_anchor_a must be finite");
  require_finite(local_anchor_b, "local_anchor_b must be finite");
  require_finite(length, "joint length must be finite");
  if (a == b || length < 0.0 ||
      (bodies_.at(idx_a).type == BodyType::static_body && bodies_.at(idx_b).type == BodyType::static_body)) {
    throw PhysicsError("invalid_joint",
                       "distance joint requires two distinct bodies, at least one dynamic, length >= 0");
  }

  DistanceJointRecord j;
  j.id = next_joint_id_++;
  j.body_a = a;
  j.body_b = b;
  j.local_anchor_a = local_anchor_a;
  j.local_anchor_b = local_anchor_b;
  j.length = length;
  joints_.push_back(j);
  return j.id;
}

void World::apply_impulse(BodyId id, Vec2 impulse, Vec2 world_point) {
  const std::size_t idx = find_body_index(id);
  require_finite(impulse, "impulse must be finite");
  require_finite(world_point, "world_point must be finite");

  BodyRecord& b = bodies_.at(idx);
  if (b.type != BodyType::dynamic_body) {
    return;
  }
  b.awake = true;
  b.calm_time = 0.0;
  const Rotation q(b.angle);
  const Vec2 world_cm = transform_point(b.position, q, b.local_centroid);
  const Vec2 r = world_point - world_cm;
  b.velocity += impulse * b.inv_mass;
  b.angular_velocity += cross(r, impulse) * b.inv_inertia;

  wake_connected(idx);
}

void World::wake_connected(std::size_t start_idx) {
  std::vector<std::size_t> queue = {start_idx};
  while (!queue.empty()) {
    const std::size_t curr = queue.back();
    queue.pop_back();
    const BodyRecord& ba = bodies_.at(curr);
    const Transform ta{.position = ba.position, .angle = ba.angle};
    for (std::size_t other = 0; other < bodies_.size(); ++other) {
      BodyRecord& bb = bodies_.at(other);
      if (other == curr || bb.type != BodyType::dynamic_body || bb.awake) {
        continue;
      }
      const Transform tb{.position = bb.position, .angle = bb.angle};
      if (collide(ba.shape, ta, bb.shape, tb).has_value()) {
        bb.awake = true;
        bb.calm_time = 0.0;
        queue.push_back(other);
      }
    }
  }
}

BodyState World::state(BodyId id) const {
  const std::size_t idx = find_body_index(id);
  const BodyRecord& b = bodies_.at(idx);
  return {
      .position = b.position,
      .angle = b.angle,
      .velocity = b.velocity,
      .angular_velocity = b.angular_velocity,
  };
}

bool World::is_awake(BodyId id) const {
  const std::size_t idx = find_body_index(id);
  return bodies_.at(idx).awake;
}

std::vector<ContactInfo> World::contacts() const { return last_contacts_; }

std::uint64_t World::state_hash() const {
  std::uint64_t hash = 14695981039346656037ULL;
  for (const auto& b : bodies_) {
    fnv1a_append_double(hash, b.position.x);
    fnv1a_append_double(hash, b.position.y);
    fnv1a_append_double(hash, b.angle);
    fnv1a_append_double(hash, b.velocity.x);
    fnv1a_append_double(hash, b.velocity.y);
    fnv1a_append_double(hash, b.angular_velocity);
  }
  return hash;
}

std::size_t World::find_body_index(BodyId id) const {
  const auto it = std::ranges::lower_bound(bodies_, id, {}, &BodyRecord::id);
  if (it == bodies_.end() || it->id != id) {
    throw PhysicsError("unknown_body", "body id does not exist");
  }
  return static_cast<std::size_t>(std::distance(bodies_.begin(), it));
}

void World::step(double dt) {
  if (!is_finite(dt) || dt <= 0.0 || dt > kMaxTimeStep) {
    throw PhysicsError("invalid_step", "step dt must be finite and in (0, 0.1]");
  }

  const std::size_t n_bodies = bodies_.size();

  // 1. Integrate velocities with gravity for awake dynamic bodies
  for (auto& b : bodies_) {
    if (b.type == BodyType::dynamic_body && b.awake) {
      b.velocity += settings_.gravity * dt;
    }
  }

  // 2. Broadphase sort-and-sweep by min_x with ties broken by body id
  struct Proxy {
    std::size_t body_idx = 0;
    BodyId id = 0;
    Aabb aabb;
  };
  std::vector<Proxy> proxies(n_bodies);
  for (std::size_t i = 0; i < n_bodies; ++i) {
    proxies.at(i) = {
        .body_idx = i,
        .id = bodies_.at(i).id,
        .aabb = compute_aabb(bodies_.at(i).shape, {.position = bodies_.at(i).position, .angle = bodies_.at(i).angle}),
    };
  }
  std::ranges::sort(proxies, [](const Proxy& a, const Proxy& b) {
    if (a.aabb.min_x != b.aabb.min_x) {
      return a.aabb.min_x < b.aabb.min_x;
    }
    return a.id < b.id;
  });

  struct Pair {
    std::size_t idx_a = 0;
    std::size_t idx_b = 0;
    BodyId id_a = 0;
    BodyId id_b = 0;
  };
  std::vector<Pair> candidate_pairs;
  for (std::size_t i = 0; i < n_bodies; ++i) {
    const Proxy& pa = proxies.at(i);
    for (std::size_t j = i + 1; j < n_bodies; ++j) {
      const Proxy& pb = proxies.at(j);
      if (pb.aabb.min_x > pa.aabb.max_x) {
        break;
      }
      if (pa.aabb.max_y < pb.aabb.min_y || pb.aabb.max_y < pa.aabb.min_y) {
        continue;
      }
      const BodyRecord& ba = bodies_.at(pa.body_idx);
      const BodyRecord& bb = bodies_.at(pb.body_idx);
      if (ba.type == BodyType::static_body && bb.type == BodyType::static_body) {
        continue;
      }
      if (pa.id < pb.id) {
        candidate_pairs.push_back({.idx_a = pa.body_idx, .idx_b = pb.body_idx, .id_a = pa.id, .id_b = pb.id});
      } else {
        candidate_pairs.push_back({.idx_a = pb.body_idx, .idx_b = pa.body_idx, .id_a = pb.id, .id_b = pa.id});
      }
    }
  }

  std::ranges::sort(candidate_pairs, [](const Pair& a, const Pair& b) {
    if (a.id_a != b.id_a) {
      return a.id_a < b.id_a;
    }
    return a.id_b < b.id_b;
  });

  // 3. Narrowphase SAT
  std::vector<ContactConstraint> constraints;
  last_contacts_.clear();
  std::vector<std::vector<std::size_t>> adjacency(n_bodies);

  for (const auto& pair : candidate_pairs) {
    const BodyRecord& ba = bodies_.at(pair.idx_a);
    const BodyRecord& bb = bodies_.at(pair.idx_b);
    const Transform ta{.position = ba.position, .angle = ba.angle};
    const Transform tb{.position = bb.position, .angle = bb.angle};
    auto manifold_opt = collide(ba.shape, ta, bb.shape, tb);
    if (!manifold_opt.has_value()) {
      continue;
    }

    adjacency.at(pair.idx_a).push_back(pair.idx_b);
    adjacency.at(pair.idx_b).push_back(pair.idx_a);

    last_contacts_.push_back({.body_a = pair.id_a, .body_b = pair.id_b, .manifold = *manifold_opt});

    const Rotation qa(ba.angle);
    const Rotation qb(bb.angle);
    const Vec2 cm_a = transform_point(ta.position, qa, ba.local_centroid);
    const Vec2 cm_b = transform_point(tb.position, qb, bb.local_centroid);

    ContactConstraint cc;
    cc.index_a = pair.idx_a;
    cc.index_b = pair.idx_b;
    cc.id_a = pair.id_a;
    cc.id_b = pair.id_b;
    cc.normal = manifold_opt->normal;
    cc.tangent = {-cc.normal.y, cc.normal.x};
    cc.penetration = manifold_opt->penetration;
    cc.friction = std::sqrt(ba.friction * bb.friction);
    cc.restitution = std::max(ba.restitution, bb.restitution);
    cc.raw_manifold = *manifold_opt;

    for (const auto& pt : manifold_opt->points) {
      ContactPointSolver cps;
      cps.point = pt;
      cps.local_a = inv_transform_point(ta.position, qa, pt);
      cps.local_b = inv_transform_point(tb.position, qb, pt);
      cps.r_a = pt - cm_a;
      cps.r_b = pt - cm_b;

      const double rn_a = cross(cps.r_a, cc.normal);
      const double rn_b = cross(cps.r_b, cc.normal);
      const double k_normal =
          ba.inv_mass + bb.inv_mass + (ba.inv_inertia * rn_a * rn_a) + (bb.inv_inertia * rn_b * rn_b);
      cps.normal_mass = (k_normal > 0.0) ? (1.0 / k_normal) : 0.0;

      const double rt_a = cross(cps.r_a, cc.tangent);
      const double rt_b = cross(cps.r_b, cc.tangent);
      const double k_tangent =
          ba.inv_mass + bb.inv_mass + (ba.inv_inertia * rt_a * rt_a) + (bb.inv_inertia * rt_b * rt_b);
      cps.tangent_mass = (k_tangent > 0.0) ? (1.0 / k_tangent) : 0.0;

      const Vec2 v_a = ba.velocity + cross(ba.angular_velocity, cps.r_a);
      const Vec2 v_b = bb.velocity + cross(bb.angular_velocity, cps.r_b);
      const double v_rel_n = dot(v_b - v_a, cc.normal);
      if (cc.restitution > 0.0 && v_rel_n < 0.0) {
        cps.velocity_bias = -cc.restitution * v_rel_n;
      } else {
        cps.velocity_bias = 0.0;
      }

      cc.points.push_back(cps);
    }
    constraints.push_back(std::move(cc));
  }

  // Wake sleeping dynamic bodies touching an awake dynamic body
  bool woke_any = true;
  while (woke_any) {
    woke_any = false;
    for (const auto& cc : constraints) {
      BodyRecord& ba = bodies_.at(cc.index_a);
      BodyRecord& bb = bodies_.at(cc.index_b);
      if (ba.type == BodyType::dynamic_body && ba.awake && bb.type == BodyType::dynamic_body && !bb.awake) {
        bb.awake = true;
        bb.calm_time = 0.0;
        woke_any = true;
      } else if (bb.type == BodyType::dynamic_body && bb.awake && ba.type == BodyType::dynamic_body && !ba.awake) {
        ba.awake = true;
        ba.calm_time = 0.0;
        woke_any = true;
      }
    }
  }

  // Precompute distance joint solver data
  struct JointSolverData {
    std::size_t idx_a = 0;
    std::size_t idx_b = 0;
    Vec2 r_a{0.0, 0.0};
    Vec2 r_b{0.0, 0.0};
    Vec2 axis{0.0, 0.0};
    double effective_mass = 0.0;
    double bias = 0.0;
  };
  std::vector<JointSolverData> joint_solvers;
  joint_solvers.reserve(joints_.size());
  for (const auto& j : joints_) {
    const std::size_t idx_a = find_body_index(j.body_a);
    const std::size_t idx_b = find_body_index(j.body_b);
    BodyRecord& ba = bodies_.at(idx_a);
    BodyRecord& bb = bodies_.at(idx_b);
    if (ba.type == BodyType::dynamic_body && ba.awake && bb.type == BodyType::dynamic_body && !bb.awake) {
      bb.awake = true;
      bb.calm_time = 0.0;
    } else if (bb.type == BodyType::dynamic_body && bb.awake && ba.type == BodyType::dynamic_body && !ba.awake) {
      ba.awake = true;
      ba.calm_time = 0.0;
    }

    const Rotation qa(ba.angle);
    const Rotation qb(bb.angle);
    const Vec2 world_a = transform_point(ba.position, qa, j.local_anchor_a);
    const Vec2 world_b = transform_point(bb.position, qb, j.local_anchor_b);
    const Vec2 cm_a = transform_point(ba.position, qa, ba.local_centroid);
    const Vec2 cm_b = transform_point(bb.position, qb, bb.local_centroid);

    JointSolverData jsd;
    jsd.idx_a = idx_a;
    jsd.idx_b = idx_b;
    jsd.r_a = world_a - cm_a;
    jsd.r_b = world_b - cm_b;
    const Vec2 d = world_b - world_a;
    const double current_len = d.length();
    jsd.axis = (current_len > 1e-12) ? (d / current_len) : Vec2{1.0, 0.0};
    const double cr_a = cross(jsd.r_a, jsd.axis);
    const double cr_b = cross(jsd.r_b, jsd.axis);
    const double k = ba.inv_mass + bb.inv_mass + (ba.inv_inertia * cr_a * cr_a) + (bb.inv_inertia * cr_b * cr_b);
    jsd.effective_mass = (k > 0.0) ? (1.0 / k) : 0.0;
    jsd.bias = (kBaumgarteFactor / dt) * (current_len - j.length);
    joint_solvers.push_back(jsd);
  }

  // 4. Velocity iterations: sequential impulses over contacts then joints in creation order
  for (int iter = 0; iter < settings_.velocity_iterations; ++iter) {
    for (auto& cc : constraints) {
      BodyRecord& ba = bodies_.at(cc.index_a);
      BodyRecord& bb = bodies_.at(cc.index_b);
      if (!ba.awake && !bb.awake) {
        continue;
      }

      for (auto& cps : cc.points) {
        // Normal impulse
        const Vec2 v_a = ba.velocity + cross(ba.angular_velocity, cps.r_a);
        const Vec2 v_b = bb.velocity + cross(bb.angular_velocity, cps.r_b);
        const double vn = dot(v_b - v_a, cc.normal);
        double d_lambda_n = -cps.normal_mass * (vn - cps.velocity_bias);
        const double old_normal_impulse = cps.normal_impulse;
        cps.normal_impulse = std::max(old_normal_impulse + d_lambda_n, 0.0);
        d_lambda_n = cps.normal_impulse - old_normal_impulse;

        const Vec2 impulse_n = cc.normal * d_lambda_n;
        ba.velocity -= impulse_n * ba.inv_mass;
        ba.angular_velocity -= cross(cps.r_a, impulse_n) * ba.inv_inertia;
        bb.velocity += impulse_n * bb.inv_mass;
        bb.angular_velocity += cross(cps.r_b, impulse_n) * bb.inv_inertia;

        // Tangent (Coulomb friction) impulse
        const Vec2 v_a_t = ba.velocity + cross(ba.angular_velocity, cps.r_a);
        const Vec2 v_b_t = bb.velocity + cross(bb.angular_velocity, cps.r_b);
        const double vt = dot(v_b_t - v_a_t, cc.tangent);
        double d_lambda_t = -cps.tangent_mass * vt;
        const double max_friction = cc.friction * cps.normal_impulse;
        const double old_tangent_impulse = cps.tangent_impulse;
        cps.tangent_impulse = std::clamp(old_tangent_impulse + d_lambda_t, -max_friction, max_friction);
        d_lambda_t = cps.tangent_impulse - old_tangent_impulse;

        const Vec2 impulse_t = cc.tangent * d_lambda_t;
        ba.velocity -= impulse_t * ba.inv_mass;
        ba.angular_velocity -= cross(cps.r_a, impulse_t) * ba.inv_inertia;
        bb.velocity += impulse_t * bb.inv_mass;
        bb.angular_velocity += cross(cps.r_b, impulse_t) * bb.inv_inertia;
      }
    }

    for (const auto& jsd : joint_solvers) {
      BodyRecord& ba = bodies_.at(jsd.idx_a);
      BodyRecord& bb = bodies_.at(jsd.idx_b);
      if (!ba.awake && !bb.awake) {
        continue;
      }
      const Vec2 v_a = ba.velocity + cross(ba.angular_velocity, jsd.r_a);
      const Vec2 v_b = bb.velocity + cross(bb.angular_velocity, jsd.r_b);
      const double c_dot = dot(v_b - v_a, jsd.axis);
      const double lambda = -jsd.effective_mass * (c_dot + jsd.bias);
      const Vec2 impulse = jsd.axis * lambda;
      ba.velocity -= impulse * ba.inv_mass;
      ba.angular_velocity -= cross(jsd.r_a, impulse) * ba.inv_inertia;
      bb.velocity += impulse * bb.inv_mass;
      bb.angular_velocity += cross(jsd.r_b, impulse) * bb.inv_inertia;
    }
  }

  // 5. Integrate positions for awake dynamic bodies
  for (auto& b : bodies_) {
    if (b.type == BodyType::dynamic_body && b.awake) {
      b.position += b.velocity * dt;
      b.angle += b.angular_velocity * dt;
    }
  }

  // 6. Position iterations: non-linear position projection (changes positions only, never velocities)
  for (int iter = 0; iter < settings_.position_iterations; ++iter) {
    for (const auto& cc : constraints) {
      BodyRecord& ba = bodies_.at(cc.index_a);
      BodyRecord& bb = bodies_.at(cc.index_b);
      if (!ba.awake && !bb.awake) {
        continue;
      }
      const Transform ta{.position = ba.position, .angle = ba.angle};
      const Transform tb{.position = bb.position, .angle = bb.angle};
      auto m_opt = collide(ba.shape, ta, bb.shape, tb);
      if (!m_opt.has_value() || m_opt->penetration <= kLinearSlop) {
        continue;
      }
      const Vec2 normal = m_opt->normal;
      const double correction = std::min(kBaumgarteFactor * (m_opt->penetration - kLinearSlop), kMaxLinearCorrection);
      const Rotation qa(ba.angle);
      const Rotation qb(bb.angle);
      const Vec2 cm_a = transform_point(ta.position, qa, ba.local_centroid);
      const Vec2 cm_b = transform_point(tb.position, qb, bb.local_centroid);

      const double inv_m_sum = ba.inv_mass + bb.inv_mass;
      if (inv_m_sum <= 0.0) {
        continue;
      }
      for (const auto& contact_pt : m_opt->points) {
        const Vec2 r_a = contact_pt - cm_a;
        const Vec2 r_b = contact_pt - cm_b;
        const double rn_a = cross(r_a, normal);
        const double rn_b = cross(r_b, normal);
        const double k = ba.inv_mass + bb.inv_mass + (ba.inv_inertia * rn_a * rn_a) + (bb.inv_inertia * rn_b * rn_b);
        if (k <= 0.0) {
          continue;
        }
        const double impulse_mag = (correction / static_cast<double>(m_opt->points.size())) / k;
        const Vec2 p = normal * impulse_mag;
        ba.position -= p * ba.inv_mass;
        ba.angle -= cross(r_a, p) * ba.inv_inertia;
        bb.position += p * bb.inv_mass;
        bb.angle += cross(r_b, p) * bb.inv_inertia;
      }
    }

    for (const auto& j : joints_) {
      const std::size_t idx_a = find_body_index(j.body_a);
      const std::size_t idx_b = find_body_index(j.body_b);
      BodyRecord& ba = bodies_.at(idx_a);
      BodyRecord& bb = bodies_.at(idx_b);
      if (!ba.awake && !bb.awake) {
        continue;
      }
      const Rotation qa(ba.angle);
      const Rotation qb(bb.angle);
      const Vec2 world_a = transform_point(ba.position, qa, j.local_anchor_a);
      const Vec2 world_b = transform_point(bb.position, qb, j.local_anchor_b);
      const Vec2 d = world_b - world_a;
      const double current_len = d.length();
      const double c_err = current_len - j.length;
      if (std::abs(c_err) <= 1e-9) {
        continue;
      }
      const Vec2 axis = (current_len > 1e-12) ? (d / current_len) : Vec2{1.0, 0.0};
      const Vec2 cm_a = transform_point(ba.position, qa, ba.local_centroid);
      const Vec2 cm_b = transform_point(bb.position, qb, bb.local_centroid);
      const Vec2 r_a = world_a - cm_a;
      const Vec2 r_b = world_b - cm_b;
      const double cr_a = cross(r_a, axis);
      const double cr_b = cross(r_b, axis);
      const double k = ba.inv_mass + bb.inv_mass + (ba.inv_inertia * cr_a * cr_a) + (bb.inv_inertia * cr_b * cr_b);
      if (k <= 0.0) {
        continue;
      }
      const double clamped_err = std::clamp(c_err, -kMaxLinearCorrection, kMaxLinearCorrection);
      const double lambda = -clamped_err / k;
      const Vec2 p = axis * lambda;
      ba.position -= p * ba.inv_mass;
      ba.angle -= cross(r_a, p) * ba.inv_inertia;
      bb.position += p * bb.inv_mass;
      bb.angle += cross(r_b, p) * bb.inv_inertia;
    }
  }

  // 7. Sleeping: islands are not built; a body sleeps when it and every body it touches have been below both speed
  // limits for sleep_time
  if (settings_.sleep_time > 0.0) {
    const double lin_tol2 = settings_.sleep_linear * settings_.sleep_linear;
    const double ang_tol = settings_.sleep_angular;
    for (auto& b : bodies_) {
      if (b.type != BodyType::dynamic_body) {
        continue;
      }
      if (!b.awake) {
        b.calm_time += dt;
        continue;
      }
      if (b.velocity.length_squared() <= lin_tol2 && std::abs(b.angular_velocity) <= ang_tol) {
        b.calm_time += dt;
      } else {
        b.calm_time = 0.0;
      }
    }

    std::vector<bool> can_sleep(n_bodies, false);
    for (std::size_t i = 0; i < n_bodies; ++i) {
      const BodyRecord& b = bodies_.at(i);
      if (b.type != BodyType::dynamic_body || !b.awake) {
        continue;
      }
      if (b.calm_time < settings_.sleep_time) {
        continue;
      }
      bool all_neighbors_calm = true;
      for (const std::size_t nbr_idx : adjacency.at(i)) {
        const BodyRecord& nbr = bodies_.at(nbr_idx);
        if (nbr.type == BodyType::dynamic_body && nbr.calm_time < settings_.sleep_time) {
          all_neighbors_calm = false;
          break;
        }
      }
      can_sleep.at(i) = all_neighbors_calm;
    }

    for (std::size_t i = 0; i < n_bodies; ++i) {
      if (can_sleep.at(i)) {
        bodies_.at(i).awake = false;
        bodies_.at(i).velocity = {0.0, 0.0};
        bodies_.at(i).angular_velocity = 0.0;
      }
    }
  }
}

}  // namespace trysac::physics_2d
