#include "trysac/physics_2d/physics.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace phys = trysac::physics_2d;

TEST_CASE("collide on boxes and circles matches acceptance criteria", "[collide]") {
  const phys::Shape box = phys::make_box(0.5, 0.5);

  // Two boxes with half extents 0.5 centered at x = 0 and x = 0.9
  const auto m1 =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, box, {.position = {0.9, 0.0}, .angle = 0.0});
  if (m1) {
    CHECK(std::abs(m1->normal.x - 1.0) <= 1e-12);
    CHECK(std::abs(m1->normal.y - 0.0) <= 1e-12);
    CHECK(std::abs(m1->penetration - 0.1) <= 1e-12);
    CHECK(m1->points.size() == 2);
  } else {
    FAIL("expected collision between boxes");
  }

  // Separated at x = 1.1
  const auto m2 =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, box, {.position = {1.1, 0.0}, .angle = 0.0});
  CHECK_FALSE(m2.has_value());

  // Circles of radius 1 at (0, 0) and (1.5, 0)
  const phys::Shape c1 = phys::Circle{1.0};
  const auto m3 = phys::collide(c1, {.position = {0.0, 0.0}, .angle = 0.0}, c1, {.position = {1.5, 0.0}, .angle = 0.0});
  if (m3) {
    CHECK(std::abs(m3->normal.x - 1.0) <= 1e-12);
    CHECK(std::abs(m3->normal.y - 0.0) <= 1e-12);
    CHECK(std::abs(m3->penetration - 0.5) <= 1e-12);
    CHECK(m3->points.size() == 1);
  } else {
    FAIL("expected collision between circles");
  }

  // Separated circles
  CHECK_FALSE(phys::collide(c1, {.position = {0.0, 0.0}, .angle = 0.0}, c1, {.position = {2.5, 0.0}, .angle = 0.0})
                  .has_value());
  // Coincident circles
  const auto m_same =
      phys::collide(c1, {.position = {0.0, 0.0}, .angle = 0.0}, c1, {.position = {0.0, 0.0}, .angle = 0.0});
  if (m_same) {
    CHECK(std::abs(m_same->penetration - 2.0) <= 1e-12);
  } else {
    FAIL("expected collision between coincident circles");
  }
}

TEST_CASE("collide handles circle-polygon in all Voronoi regions", "[collide]") {
  const phys::Shape box = phys::make_box(1.0, 1.0);
  const phys::Shape circle = phys::Circle{0.5};

  // Face contact (poly A, circle B)
  const auto m_face =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {1.25, 0.0}, .angle = 0.0});
  if (m_face) {
    CHECK(std::abs(m_face->normal.x - 1.0) <= 1e-12);
    CHECK(std::abs(m_face->penetration - 0.25) <= 1e-12);
  } else {
    FAIL("expected collision in face contact");
  }

  // Reversed order (circle A, poly B)
  const auto m_rev =
      phys::collide(circle, {.position = {1.25, 0.0}, .angle = 0.0}, box, {.position = {0.0, 0.0}, .angle = 0.0});
  if (m_rev) {
    CHECK(std::abs(m_rev->normal.x + 1.0) <= 1e-12);
    CHECK(std::abs(m_rev->penetration - 0.25) <= 1e-12);
  } else {
    FAIL("expected collision in reversed face contact");
  }

  // Vertex v1 region
  const auto m_v1 =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {1.2, 1.2}, .angle = 0.0});
  if (m_v1) {
    CHECK(m_v1->penetration > 0.0);
  } else {
    FAIL("expected collision in v1 region");
  }

  // Vertex v2 region
  const auto m_v2 =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {1.2, -1.2}, .angle = 0.0});
  if (m_v2) {
    CHECK(m_v2->penetration > 0.0);
  } else {
    FAIL("expected collision in v2 region");
  }

  // Circle center deeply inside box
  const auto m_inside =
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {0.1, 0.0}, .angle = 0.0});
  if (m_inside) {
    CHECK(m_inside->penetration > 0.5);
  } else {
    FAIL("expected collision inside box");
  }

  // Separated in vertex regions
  CHECK_FALSE(
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {1.45, 1.45}, .angle = 0.0})
          .has_value());
  CHECK_FALSE(
      phys::collide(box, {.position = {0.0, 0.0}, .angle = 0.0}, circle, {.position = {1.45, -1.45}, .angle = 0.0})
          .has_value());
}

TEST_CASE("Elastic head-on collision conserves kinetic energy and linear momentum", "[world][elastic]") {
  phys::WorldSettings ws;
  ws.gravity = {0.0, 0.0};
  phys::World world(ws);

  phys::BodyDef def_a;
  def_a.type = phys::BodyType::dynamic_body;
  def_a.position = {-1.0, 0.0};
  def_a.velocity = {2.0, 0.0};
  def_a.shape = phys::Circle{0.5};
  def_a.density = 1.0;
  def_a.restitution = 1.0;
  def_a.friction = 0.0;

  phys::BodyDef def_b = def_a;
  def_b.position = {1.0, 0.0};
  def_b.velocity = {-2.0, 0.0};

  const auto id_a = world.create_body(def_a);
  const auto id_b = world.create_body(def_b);

  const double mass = std::acos(-1.0) * 0.5 * 0.5 * 1.0;
  const double initial_ke = (0.5 * mass * 4.0) + (0.5 * mass * 4.0);
  const double dt = 1.0 / 60.0;

  for (int i = 0; i < 60; ++i) {
    world.step(dt);
    const auto sa = world.state(id_a);
    const auto sb = world.state(id_b);
    const phys::Vec2 total_p = (sa.velocity + sb.velocity) * mass;
    CHECK(std::abs(total_p.x) <= 1e-12);
    CHECK(std::abs(total_p.y) <= 1e-12);
  }

  const auto final_a = world.state(id_a);
  const auto final_b = world.state(id_b);
  CHECK(std::abs(final_a.velocity.x - (-2.0)) <= 1e-9);
  CHECK(std::abs(final_a.velocity.y - 0.0) <= 1e-9);
  CHECK(std::abs(final_b.velocity.x - 2.0) <= 1e-9);
  CHECK(std::abs(final_b.velocity.y - 0.0) <= 1e-9);

  const double final_ke =
      (0.5 * mass * final_a.velocity.length_squared()) + (0.5 * mass * final_b.velocity.length_squared());
  CHECK(std::abs(final_ke - initial_ke) <= 1e-9);
}

TEST_CASE("Free fall follows semi-implicit Euler integration", "[world][freefall]") {
  phys::WorldSettings ws;
  ws.gravity = {0.0, -10.0};
  phys::World world(ws);

  phys::BodyDef def;
  def.type = phys::BodyType::dynamic_body;
  def.position = {0.0, 10.0};
  def.velocity = {0.0, 0.0};
  def.shape = phys::Circle{0.5};
  const auto id = world.create_body(def);

  for (int i = 0; i < 100; ++i) {
    world.step(0.01);
  }

  const auto st = world.state(id);
  CHECK(std::abs(st.velocity.y - (-10.0)) <= 1e-9);
  CHECK(std::abs(st.position.y - 4.95) <= 1e-9);
}

TEST_CASE("Resting and sleeping behavior and waking via apply_impulse", "[world][sleep]") {
  phys::World world;

  phys::BodyDef ground_def;
  ground_def.type = phys::BodyType::static_body;
  ground_def.position = {0.0, -0.5};
  ground_def.shape = phys::make_box(10.0, 0.5);
  ground_def.restitution = 0.0;
  ground_def.friction = 0.5;
  world.create_body(ground_def);

  phys::BodyDef box_def;
  box_def.type = phys::BodyType::dynamic_body;
  box_def.position = {0.0, 0.6};
  box_def.shape = phys::make_box(0.5, 0.5);
  box_def.restitution = 0.0;
  box_def.friction = 0.5;
  const auto box_id = world.create_body(box_def);

  phys::BodyDef top_box_def = box_def;
  top_box_def.position = {0.0, 1.6};
  const auto top_box_id = world.create_body(top_box_def);

  const double dt = 1.0 / 60.0;
  for (int step = 1; step <= 120; ++step) {
    world.step(dt);
  }
  CHECK(std::abs(world.state(box_id).velocity.y) < 0.001);

  for (int step = 121; step <= 240; ++step) {
    world.step(dt);
  }
  CHECK_FALSE(world.is_awake(box_id));
  CHECK_FALSE(world.is_awake(top_box_id));
  CHECK_FALSE(world.contacts().empty());

  // apply_impulse wakes the box and the touching top box
  world.apply_impulse(box_id, {0.0, 5.0}, {0.0, 0.5});
  CHECK(world.is_awake(box_id));
  CHECK(world.is_awake(top_box_id));
}

TEST_CASE("Determinism across 50 bodies with seed 99 for 600 steps", "[world][determinism]") {
  const auto run_scene = [] {
    phys::World world;
    phys::BodyDef ground;
    ground.type = phys::BodyType::static_body;
    ground.position = {0.0, -1.0};
    ground.shape = phys::make_box(20.0, 1.0);
    world.create_body(ground);

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937_64 rng(99);
    std::uniform_real_distribution<double> x_dist(-5.0, 5.0);
    std::uniform_real_distribution<double> y_dist(1.0, 15.0);
    std::uniform_real_distribution<double> v_dist(-1.0, 1.0);

    for (int i = 0; i < 49; ++i) {
      phys::BodyDef d;
      d.type = phys::BodyType::dynamic_body;
      d.position = {x_dist(rng), y_dist(rng)};
      d.velocity = {v_dist(rng), v_dist(rng)};
      if (i % 2 == 0) {
        d.shape = phys::Circle{0.25};
      } else {
        d.shape = phys::make_box(0.25, 0.25);
      }
      d.restitution = 0.2;
      d.friction = 0.3;
      world.create_body(d);
    }

    std::vector<std::uint64_t> hashes;
    hashes.reserve(600);
    for (int s = 0; s < 600; ++s) {
      world.step(1.0 / 60.0);
      hashes.push_back(world.state_hash());
    }
    return hashes;
  };

  const auto hashes_1 = run_scene();
  const auto hashes_2 = run_scene();
  REQUIRE(hashes_1.size() == 600);
  CHECK(hashes_1 == hashes_2);
}

TEST_CASE("Distance joint maintains pendulum anchor distance within 0.001 for 300 steps", "[world][joint]") {
  phys::World world;

  phys::BodyDef anchor_def;
  anchor_def.type = phys::BodyType::static_body;
  anchor_def.position = {0.0, 5.0};
  anchor_def.shape = phys::Circle{0.1};
  const auto anchor_id = world.create_body(anchor_def);

  phys::BodyDef bob_def;
  bob_def.type = phys::BodyType::dynamic_body;
  bob_def.position = {2.0, 5.0};
  bob_def.shape = phys::Circle{0.2};
  const auto bob_id = world.create_body(bob_def);

  world.add_distance_joint(anchor_id, bob_id, {0.0, 0.0}, {0.0, 0.0}, 2.0);

  const double dt = 1.0 / 60.0;
  for (int i = 0; i < 300; ++i) {
    world.step(dt);
    const auto bob_state = world.state(bob_id);
    const double dist = (bob_state.position - phys::Vec2{0.0, 5.0}).length();
    CHECK(std::abs(dist - 2.0) <= 0.001);
  }
}

TEST_CASE("Error codes and validation checks", "[errors]") {
  const auto check_code = [](auto&& fn, std::string_view expected_code) {
    bool caught = false;
    try {
      fn();
    } catch (const phys::PhysicsError& err) {
      caught = true;
      CHECK(err.code() == expected_code);
    }
    CHECK(caught);
  };

  // Dynamic body density <= 0 -> invalid_mass
  check_code(
      [] {
        phys::World w;
        phys::BodyDef d;
        d.type = phys::BodyType::dynamic_body;
        d.density = 0.0;
        w.create_body(d);
      },
      "invalid_mass");

  check_code(
      [] {
        phys::World w;
        phys::BodyDef d;
        d.type = phys::BodyType::dynamic_body;
        d.density = -1.0;
        w.create_body(d);
      },
      "invalid_mass");

  // Polygon with 9 vertices -> too_many_vertices
  check_code(
      [] {
        std::vector<phys::Vec2> verts(9);
        for (int i = 0; i < 9; ++i) {
          const double a = 2.0 * std::numbers::pi * static_cast<double>(i) / 9.0;
          verts.at(static_cast<std::size_t>(i)) = {std::cos(a), std::sin(a)};
        }
        const phys::Polygon p(verts);
      },
      "too_many_vertices");

  // Non-convex polygon or collinear vertices -> not_convex
  check_code([] { const phys::Polygon p({{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {1.0, 1.0}}); }, "not_convex");

  check_code(
      [] {
        // Star / concave shape
        const phys::Polygon p({{0.0, 0.0}, {2.0, 0.0}, {1.0, 0.5}, {2.0, 2.0}, {0.0, 2.0}});
      },
      "not_convex");

  check_code([] { const phys::Polygon p({{0.0, 0.0}, {1.0, 0.0}}); }, "not_convex");

  // Clockwise convex polygon is automatically reversed to CCW
  const phys::Polygon cw_poly({{0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}});
  CHECK(cw_poly.vertices.size() == 4);

  // step(0) or step(0.2) -> invalid_step
  check_code(
      [] {
        phys::World w;
        w.step(0.0);
      },
      "invalid_step");
  check_code(
      [] {
        phys::World w;
        w.step(0.2);
      },
      "invalid_step");

  // 10_001st body -> too_many_bodies
  check_code(
      [] {
        phys::World w;
        phys::BodyDef d;
        d.type = phys::BodyType::static_body;
        for (std::size_t i = 0; i < 10'000; ++i) {
          w.create_body(d);
        }
        w.create_body(d);
      },
      "too_many_bodies");

  // NaN position -> invalid_value
  check_code(
      [] {
        phys::World w;
        phys::BodyDef d;
        d.position = {std::numeric_limits<double>::quiet_NaN(), 0.0};
        w.create_body(d);
      },
      "invalid_value");

  // Unknown body and invalid joint
  check_code(
      [] {
        const phys::World w;
        static_cast<void>(w.state(999));
      },
      "unknown_body");

  check_code(
      [] {
        phys::World w;
        const phys::BodyDef d;
        const auto id = w.create_body(d);
        w.add_distance_joint(id, id, {0.0, 0.0}, {0.0, 0.0}, 1.0);
      },
      "invalid_joint");

  // Destroy body removes joints and body
  {
    phys::World w;
    const phys::BodyDef d;
    const auto a = w.create_body(d);
    const auto b = w.create_body(d);
    w.add_distance_joint(a, b, {0.0, 0.0}, {0.0, 0.0}, 1.0);
    w.destroy_body(a);
    check_code([&] { static_cast<void>(w.state(a)); }, "unknown_body");
  }
}

TEST_CASE("Vec2 arithmetic and properties", "[vec2]") {
  const double x = GENERATE(take(10, random(-100.0, 100.0)));
  const double y = GENERATE(take(10, random(-100.0, 100.0)));
  phys::Vec2 v{x, y};
  const phys::Vec2 neg = -v;
  CHECK((v + neg) == phys::Vec2{0.0, 0.0});
  CHECK(phys::dot(v, {1.0, 0.0}) == x);
  CHECK(phys::cross(v, v) == 0.0);
  v += {1.0, 2.0};
  v -= {1.0, 2.0};
  v *= 2.0;
  CHECK(v == (2.0 * phys::Vec2{x, y}));
  CHECK((phys::Vec2{0.0, 0.0}).normalized() == phys::Vec2{0.0, 0.0});
}
