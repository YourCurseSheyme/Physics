//! Second task

use nalgebra::Vector2;

#[derive(Clone, Copy, Debug)]
pub struct Environment {
  pub wall_normal: Vector2<f64>,
}

impl Environment {
  pub fn new(wall_normal: Vector2<f64>) -> Self {
    let n = unit(wall_normal);
    Self { wall_normal: n }
  }
}

#[derive(Clone, Copy, Debug)]
pub struct Ball {
  pub m: f64,
  pub x: Vector2<f64>,
  pub v: Vector2<f64>,
}

impl Ball {
  pub fn new(m: f64, x: Vector2<f64>, v: Vector2<f64>) -> Self {
    Self { m, x, v }
  }
}

#[derive(Clone, Copy, Debug)]
pub enum ContactLaw {
  Hooke { k: f64 }, // F = -k * delta
  Hertz { k: f64 }, // F = -k * delta^(3/2)
}

pub struct BilliardsSimulator {
  env: Environment,
  dt: f64,
  max_steps: usize,
}

impl BilliardsSimulator {
  pub fn new(env: Environment) -> Self {
    Self { env, dt: 1e-5, max_steps: 2_000_000 }
  }

  pub fn with_time_step(mut self, dt: f64) -> Self {
    self.dt = dt;
    self
  }

  pub fn with_max_steps(mut self, max_steps: usize) -> Self {
    self.max_steps = max_steps;
    self
  }

  // Ball-Wall elastic
  pub fn collide_wall_elastic(&self, ball: &Ball) -> Vector2<f64> {
    let n = self.env.wall_normal;
    let v = ball.v;
    v - 2.0 * v.dot(&n) * n
  }

  /// Ball-Ball elastic
  pub fn collide_balls_elastic(&self, b1: &Ball, b2: &Ball) -> (Vector2<f64>, Vector2<f64>) {
    let n = unit(b2.x - b1.x);
    let t = Vector2::new(-n.y, n.x);

    let u1n = b1.v.dot(&n);
    let u1t = b1.v.dot(&t);
    let u2n = b2.v.dot(&n);
    let u2t = b2.v.dot(&t);

    let m1 = b1.m;
    let m2 = b2.m;
    let denom = m1 + m2;

    let v1n = ((m1 - m2) / denom) * u1n + (2.0 * m2 / denom) * u2n;
    let v2n = (2.0 * m1 / denom) * u1n + ((m2 - m1) / denom) * u2n;

    let v1 = v1n * n + u1t * t;
    let v2 = v2n * n + u2t * t;
    (v1, v2)
  }

  // Ball-Wall
  pub fn collide_wall_contact(
    &self,
    ball: &Ball,
    law: ContactLaw,
  ) -> (Vector2<f64>, f64, f64) {
    let n = self.env.wall_normal;
    let s = ball.v.dot(&n);
    if s >= 0.0 {
      return (ball.v, 0.0, 0.0);
    }

    let v_in = -s;
    let m_eff = ball.m;

    let (ddelta_out, t_c, delta_max) =
      integrate_contact(v_in, law, m_eff, self.dt, self.max_steps);

    let s_out = -ddelta_out;  // >= 0
    let j = m_eff * (s_out - s);

    let v_out = ball.v + (j / ball.m) * n;
    (v_out, t_c, delta_max)
  }


  /// Ball-Ball
  pub fn collide_balls_contact(
    &self,
    b1: &Ball,
    b2: &Ball,
    law: ContactLaw,
  ) -> (Vector2<f64>, Vector2<f64>, f64, f64) {
    let n = unit(b2.x - b1.x);
    let v_rel = b2.v - b1.v;
    let s = v_rel.dot(&n);
    if s >= 0.0 {
      return (b1.v, b2.v, 0.0, 0.0);
    }

    let v_in = -s;
    let m_eff = (b1.m * b2.m) / (b1.m + b2.m);

    let (ddelta_out, t_c, delta_max) =
      integrate_contact(v_in, law, m_eff, self.dt, self.max_steps);
    let s_out = -ddelta_out;

    let j = m_eff * (s - s_out);
    let v1_out = b1.v + ( j / b1.m) * n;
    let v2_out = b2.v - ( j / b2.m) * n;
    (v1_out, v2_out, t_c, delta_max)
  }

  pub fn check_energy_momentum_two_balls(
    &self,
    before: (&Ball, &Ball),
    after_v: (Vector2<f64>, Vector2<f64>)
  ) -> (f64, f64, f64) {
    let (b1, b2) = before;
    let (v1p, v2p) = after_v;

    let e0 = kinetic_energy(b1.m, b1.v) + kinetic_energy(b2.m, b2.v);
    let e1 = kinetic_energy(b1.m, v1p) + kinetic_energy(b2.m, v2p);
    let de_rel = rel(e1 - e0, e0);

    let p0 = b1.m * b1.v + b2.m * b2.v;
    let p1 = b1.m * v1p + b2.m * v2p;
    let dp_rel = rel((p1 - p0).norm(), p0.norm());

    (e0, de_rel, dp_rel)
  }

  pub fn aim_ghost_ball(
    &self,
    cue_center: Vector2<f64>,
    target_center: Vector2<f64>,
    pocket_pos: Vector2<f64>,
    r: f64,
  ) -> Option<(Vector2<f64>, Vector2<f64>)> {
    let u = pocket_pos - target_center;
    let d = u.norm();
    if d <= 1e-12 { return None; }
    let uhat = u / d;
    let ghost = target_center - 2.0 * r * uhat;
    let dir = ghost - cue_center;
    let nd = dir.norm();
    if nd <= 1e-12 { return None; }
    Some((ghost, dir / nd))
  }
}

#[inline]
fn unit(v: Vector2<f64>) -> Vector2<f64> {
  let n = v.norm();
  if n == 0.0 { Vector2::new(1.0, 0.0) } else { v / n }
}

#[inline]
fn kinetic_energy(m: f64, v: Vector2<f64>) -> f64 {
  0.5 * m * v.norm_squared()
}

#[inline]
fn rel(x: f64, scale: f64) -> f64 {
  let s = scale.abs().max(1e-12);
  (x / s).abs()
}

fn a_from_delta(law: ContactLaw, m_eff: f64, delta: f64) -> f64 {
  match law {
    ContactLaw::Hooke { k } => -(k / m_eff) * delta,
    ContactLaw::Hertz { k } => {
      if delta <= 0.0 { 0.0 } else { -(k / m_eff) * delta.powf(1.5) }
    }
  }
}

fn integrate_contact(
  v_n_in: f64,
  law: ContactLaw,
  m_eff: f64,
  dt: f64,
  max_steps: usize,
) -> (f64, f64, f64) {
  // δ(0)=0, δ̇(0)=v_n_in (>0)
  let mut delta = 0.0;
  let mut ddelta = v_n_in.max(0.0);
  let mut a = a_from_delta(law, m_eff, delta);

  let mut t = 0.0;
  let mut delta_max = 0.0;

  for _ in 0..max_steps {
    let mut delta_new = delta + ddelta * dt + 0.5 * a * dt * dt;
    if delta_new < 0.0 { delta_new = 0.0; }
    let a_new = a_from_delta(law, m_eff, delta_new);
    let ddelta_new = ddelta + 0.5 * (a + a_new) * dt;
    if delta_new > delta_max { delta_max = delta_new; }

    t += dt;
    if delta > 0.0 && delta_new == 0.0 && ddelta_new <= 0.0 {
      return (ddelta_new, t, delta_max);
    }

    delta = delta_new;
    ddelta = ddelta_new;
    a = a_new;
  }
  (ddelta, t, delta_max)
}

fn test_elastic_two_balls(sim: &BilliardsSimulator) {
  println!("Test #1: Elastic (analytic) two balls");
  let b1 = Ball::new(1.0, Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0));
  let b2 = Ball::new(2.0, Vector2::new(1.0, 0.0), Vector2::new(0.0, 0.0));

  let (v1p, v2p) = sim.collide_balls_elastic(&b1, &b2);
  let (_e0, de_rel, dp_rel) = sim.check_energy_momentum_two_balls((&b1, &b2), (v1p, v2p));

  println!("v1' = {:?}, v2' = {:?}", v1p, v2p);
  println!("Energy drift (should be 0): {:.3e}, Momentum drift: {:.3e}", de_rel, dp_rel);
}

fn test_hooke_contact_vs_elastic(sim: &BilliardsSimulator) {
  println!("Test #2: Hooke contact vs analytic");
  let b1 = Ball::new(1.0, Vector2::new(0.0, 0.0), Vector2::new(1.0, 0.0));
  let b2 = Ball::new(2.0, Vector2::new(1.0, 0.0), Vector2::new(0.0, 0.0));

  let (v1_a, v2_a) = sim.collide_balls_elastic(&b1, &b2);

  let law = ContactLaw::Hooke { k: 1e5 };
  let (v1_c, v2_c, t_c, delta_max) = sim.collide_balls_contact(&b1, &b2, law);

  let diff_v1 = (v1_c - v1_a).norm();
  let diff_v2 = (v2_c - v2_a).norm();

  let (_e0, de_rel, dp_rel) = sim.check_energy_momentum_two_balls((&b1, &b2), (v1_c, v2_c));

  println!("analytic v1'={:?}, v2'={:?}", v1_a, v2_a);
  println!("contact  v1'={:?}, v2'={:?}, t_contact={:.6} s, delta_max={:.6} m", v1_c, v2_c, t_c, delta_max);
  println!("|Δv1|={:.3e}, |Δv2|={:.3e}; energy drift={:.3e}, momentum drift={:.3e}",
           diff_v1, diff_v2, de_rel, dp_rel);
}

fn test_hertz_contact_vs_elastic(sim: &BilliardsSimulator) {
  println!("Test #3: Hertz contact vs analytic");
  let b1 = Ball::new(1.0, Vector2::new(0.0, 0.0), Vector2::new(1.5, 0.0));
  let b2 = Ball::new(1.0, Vector2::new(0.3, 0.0), Vector2::new(0.0, 0.0));

  let (v1_a, v2_a) = sim.collide_balls_elastic(&b1, &b2);

  let sim = BilliardsSimulator::new(sim.env).with_time_step(5e-6).with_max_steps(sim.max_steps);
  let law = ContactLaw::Hertz { k: 2e6 };
  let (v1_c, v2_c, t_c, delta_max) = sim.collide_balls_contact(&b1, &b2, law);

  let diff_v1 = (v1_c - v1_a).norm();
  let diff_v2 = (v2_c - v2_a).norm();

  let (_e0, de_rel, dp_rel) = sim.check_energy_momentum_two_balls((&b1, &b2), (v1_c, v2_c));

  println!("analytic v1'={:?}, v2'={:?}", v1_a, v2_a);
  println!("contact  v1'={:?}, v2'={:?}, t_contact={:.6} s, delta_max={:.6} m", v1_c, v2_c, t_c, delta_max);
  println!("|Δv1|={:.3e}, |Δv2|={:.3e}; energy drift={:.3e}, momentum drift={:.3e}",
           diff_v1, diff_v2, de_rel, dp_rel);
}

fn test_wall_elastic() {
  println!("Test #4: Wall (elastic)");
  let sim_wall = BilliardsSimulator::new(Environment::new(Vector2::new(1.0, 0.0)));

  let speed = 1.2345;
  let ball = Ball::new(1.0, Vector2::new(1.0, 0.0), Vector2::new(-speed, 0.3));

  let v_ref = sim_wall.collide_wall_elastic(&ball);

  let v0 = ball.v.norm();
  let v1 = v_ref.norm();
  let dv_abs = (v1 - v0).abs();

  println!("v_in  = {:?}", ball.v);
  println!("v_out = {:?}", v_ref);
  println!("|v_out| - |v_in| = {:.3e} (≈ 0)", dv_abs);
}

fn test_wall_contact_vs_elastic() {
  println!("Test #5: Wall contact (Hooke/Hertz) vs elastic");
  let sim_wall = BilliardsSimulator::new(Environment::new(Vector2::new(1.0, 0.0)))
      .with_time_step(2e-6)
      .with_max_steps(5_000_000);

  let ball = Ball::new(1.0, Vector2::new(0.5, 0.2), Vector2::new(-1.0, 0.0));

  let v_el = sim_wall.collide_wall_elastic(&ball);

  // Hooke
  let (v_hooke, t_h, dmax_h) = sim_wall.collide_wall_contact(&ball, ContactLaw::Hooke { k: 2e5 });
  let diff_hooke = (v_hooke - v_el).norm();

  // Hertz
  let sim_wall_hertz = BilliardsSimulator::new(sim_wall.env).with_time_step(1e-6).with_max_steps(sim_wall.max_steps);
  let (v_hertz, t_z, dmax_z) = sim_wall_hertz.collide_wall_contact(&ball, ContactLaw::Hertz { k: 5e6 });
  let diff_hertz = (v_hertz - v_el).norm();

  println!("elastic v'   = {:?}", v_el);
  println!("Hooke  v'    = {:?}  (Δ= {:.3e}),  t_contact={:.6} s, δ_max={:.6} m", v_hooke, diff_hooke, t_h, dmax_h);
  println!("Hertz  v'    = {:?}  (Δ= {:.3e}),  t_contact={:.6} s, δ_max={:.6} m", v_hertz, diff_hertz, t_z, dmax_z);
}

fn test_ghost_ball(sim: &BilliardsSimulator) {
  println!("Test #6: Ghost ball aiming");

  let r = 0.05715;
  let target = Vector2::new(0.5, 0.2);
  let pocket = Vector2::new(1.2, 0.2);
  let cue_start = Vector2::new(0.1, 0.0);

  let Some((ghost, aim_dir)) = sim.aim_ghost_ball(cue_start, target, pocket, r) else {
    println!("ghost_ball: geometry degenerate");
    return;
  };

  println!("ghost center = {:?}", ghost);
  println!("aim dir      = {:?}", aim_dir);

  let m = 1.0;
  let hit_dir = (target - ghost).normalize();
  let b1 = Ball::new(m, ghost, 1.0 * hit_dir);
  let b2 = Ball::new(m, target, Vector2::new(0.0, 0.0));

  let (v1_after, v2_after) = sim.collide_balls_elastic(&b1, &b2);

  let to_pocket = (pocket - target).normalize();
  let proj_mag = v2_after.dot(&to_pocket);
  let perp = (v2_after - proj_mag * to_pocket).norm();

  println!("v2_after      = {:?}", v2_after);
  println!("projection    = {:.6} m/s along pocket line,  perpendicular drift = {:.3e} m/s", proj_mag, perp);

  println!("v1_after norm = {:.6} m/s (should be ~0 for equal masses)", v1_after.norm());
}

pub fn run_billiard_tests() -> Result<(), Box<dyn std::error::Error>> {
  println!("> Running billiards tests...");
  let env = Environment::new(Vector2::new(0.0, 1.0));
  let sim = BilliardsSimulator::new(env)
    .with_time_step(1e-5)
    .with_max_steps(2_000_000);

  println!("+---------------------------");
  test_elastic_two_balls(&sim);
  println!("+---------------------------");
  test_hooke_contact_vs_elastic(&sim);
  println!("+---------------------------");
  test_hertz_contact_vs_elastic(&sim);
  println!("+---------------------------");
  test_wall_elastic();
  println!("+---------------------------");
  test_wall_contact_vs_elastic();
  println!("+---------------------------");
  test_ghost_ball(&sim);
  println!("+---------------------------");

  println!("> Done!");
  Ok(())
}
