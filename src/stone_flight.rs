//! First task — refactored

use ode_solvers::*;
use plotters::prelude::*;

const ASSETS_PATH: &str = "plots/m1.png";

#[derive(Debug, Clone, Copy)]
pub enum DragModel {
  NoDrag,
  Viscous(f64), // F = -k v
  Frontal(f64), // F = -k |v| v
}

#[derive(Debug, Clone, Copy)]
pub struct SimulationParameters {
  pub mass: f64,
  pub g: f64,
  pub drag_model: DragModel,
}

type R2 = Vector2<f64>;
type R4 = Vector4<f64>;

struct ProjectileSystem {
  params: SimulationParameters,
}

impl System<f64, R4> for ProjectileSystem {
  fn system(&self, _t: f64, y: &R4, dy: &mut R4) {
    let v = R2::new(y[2], y[3]);

    let gravity = R2::new(0.0, -self.params.mass * self.params.g);
    let drag = match self.params.drag_model {
      DragModel::NoDrag => R2::zeros(),
      DragModel::Viscous(k) => -k * v,
      DragModel::Frontal(k) => -k * v.norm() * v,
    };
    let a = (gravity + drag) / self.params.mass;

    dy[0] = v.x;
    dy[1] = v.y;
    dy[2] = a.x;
    dy[3] = a.y;
  }

  fn solout(&mut self, _t: f64, y: &R4, _dy: &R4) -> bool {
    y[1] < 0.0
  }
}

fn hit_ground_interpolate(p1: (f64, f64), p2: (f64, f64)) -> (f64, f64) {
  let (x1, y1) = p1;
  let (x2, y2) = p2;
  let denom = y1 - y2;
  if denom.abs() < f64::EPSILON {
    return (x1, 0.0);
  }
  let t = (0.0 - y1) / denom; // y1 + t*(y2-y1) = 0
  (x1 + t * (x2 - x1), 0.0)
}

fn clamp_to_ground(mut xy: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
  if xy.len() >= 2 {
    let n = xy.len();
    let y_last = xy[n - 1].1;
    let y_prev = xy[n - 2].1;
    if y_last < 0.0 && y_prev >= 0.0 {
      xy[n - 1] = hit_ground_interpolate(xy[n - 2], xy[n - 1]);
    } else if y_last < 0.0 {
      // если обе последние ниже нуля, ищем последнее пересечение
      for i in (1..n).rev() {
        if xy[i - 1].1 >= 0.0 && xy[i].1 < 0.0 {
          let h = hit_ground_interpolate(xy[i - 1], xy[i]);
          xy.truncate(i); // обрежем всё после пересечения
          xy.push(h);
          break;
        }
      }
    }
  }
  xy
}

pub struct StoneFlightSimulator {
  initial_velocity: f64,
  initial_angle_deg: f64,
  params: SimulationParameters,
  atol: f64,
  rtol: f64,
  h0: f64,
}

impl StoneFlightSimulator {
  pub fn new(initial_velocity: f64, initial_angle_deg: f64, params: SimulationParameters) -> Self {
    Self {
      initial_velocity,
      initial_angle_deg,
      params,
      atol: 1.0e-6,
      rtol: 1.0e-6,
      h0: 0.05,
    }
  }

  pub fn calculate_trajectory(&self) -> Result<Vec<(f64, f64)>, String> {
    let angle_rad = self.initial_angle_deg.to_radians();
    let v0x = self.initial_velocity * angle_rad.cos();
    let v0y = self.initial_velocity * angle_rad.sin();

    // оценим разумный t_end: без сопротивления T ≈ 2 v0 sinθ / g
    let t_est = if self.params.g > 0.0 {
      (2.0 * self.initial_velocity * angle_rad.sin() / self.params.g).max(1.0)
    } else {
      5.0
    };
    let t_end = t_est * 1.5;

    let y0 = R4::new(0.0, 0.0, v0x, v0y);
    let system = ProjectileSystem { params: self.params };
    let mut stepper =
      Dopri5::new(system, 0.0, t_end, self.h0, y0, self.atol, self.rtol);

    match stepper.integrate() {
      Ok(_stats) => {
        let states = stepper.y_out();
        if states.is_empty() {
          return Err("empty integration output".into());
        }
        let mut pts: Vec<(f64, f64)> = states.iter().map(|s| (s[0], s[1])).collect();
        pts = clamp_to_ground(pts);
        if pts.is_empty() {
          return Err("no trajectory points after clamping".into());
        }
        Ok(pts)
      }
      Err(e) => Err(format!("integrator error: {e}")),
    }
  }

  pub fn calculate_trajectory_by_curvature(&self) -> Result<Vec<(f64, f64)>, String> {
    let mut trajectory = Vec::new();
    let angle_rad = self.initial_angle_deg.to_radians();
    let mut pos = Vector2::new(0.0, 0.0);
    let mut vel = Vector2::new(
      self.initial_velocity * angle_rad.cos(),
      self.initial_velocity * angle_rad.sin(),
    );
    trajectory.push((pos.x, pos.y));

    let d_theta = 0.005;
    let max_steps = 5000;

    for _ in 0..max_steps {
      if pos.y < 0.0 { break; }

      let v_mag = vel.norm();
      if v_mag < 1e-6 { break; }

      let gravity_force = Vector2::new(0.0, -self.params.mass * self.params.g);
      let drag_force = match self.params.drag_model {
        DragModel::NoDrag => Vector2::zeros(),
        DragModel::Viscous(k) => -k * vel,
        DragModel::Frontal(k) => -k * v_mag * vel,
      };
      let total_force = gravity_force + drag_force;
      let acc = total_force / self.params.mass;

      let vel_unit = vel / v_mag;
      let acc_tangential_scalar = acc.dot(&vel_unit);
      let acc_tangential_vec = acc_tangential_scalar * vel_unit;
      let acc_normal_vec = acc - acc_tangential_vec;
      let acc_normal_mag = acc_normal_vec.norm();

      if acc_normal_mag < 1e-6 {
        let dt = 0.01;
        pos += vel * dt;
        vel += acc * dt;
        trajectory.push((pos.x, pos.y));
        continue;
      }

      let turn_sign = (vel.x * acc.y - vel.y * acc.x).signum();
      if turn_sign == 0.0 {
        let dt = 0.01;
        pos += vel * dt;
        vel += acc * dt;
        trajectory.push((pos.x, pos.y));
        continue;
      }
      let signed_d_theta = d_theta * turn_sign;

      let radius = v_mag.powi(2) / acc_normal_mag;
      let normal_unit = acc_normal_vec / acc_normal_mag;
      let center = pos + radius * normal_unit;
      let r_vec = pos - center;

      let cos_dt = signed_d_theta.cos();
      let sin_dt = signed_d_theta.sin();
      let new_r_vec = Vector2::new(
        r_vec.x * cos_dt - r_vec.y * sin_dt,
        r_vec.x * sin_dt + r_vec.y * cos_dt,
      );
      pos = center + new_r_vec;
      trajectory.push((pos.x, pos.y));

      vel = Vector2::new(
        vel.x * cos_dt - vel.y * sin_dt,
        vel.x * sin_dt + vel.y * cos_dt,
      );

      let dt = (radius * d_theta) / v_mag;
      let new_v_mag = v_mag + acc_tangential_scalar * dt;
      vel = vel.normalize() * new_v_mag;
    }

    Ok(trajectory)
  }
}

fn test_no_drag(initial_velocity: f64, initial_angle_deg: f64) -> Result<Vec<(f64, f64)>, String> {
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::NoDrag,
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  simulator.calculate_trajectory()
}

fn test_viscous_drag(
  initial_velocity: f64,
  initial_angle_deg: f64,
) -> Result<Vec<(f64, f64)>, String> {
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::Viscous(0.1),
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  simulator.calculate_trajectory()
}

fn test_frontal_drag(
  initial_velocity: f64,
  initial_angle_deg: f64,
) -> Result<Vec<(f64, f64)>, String> {
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::Frontal(0.01),
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  simulator.calculate_trajectory()
}

fn test_curvature_method(
  initial_velocity: f64,
  initial_angle_deg: f64,
) -> Result<Vec<(f64, f64)>, String> {
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::NoDrag,
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  simulator.calculate_trajectory_by_curvature()
}

pub fn run_stone_tests() -> Result<(), Box<dyn std::error::Error>> {
  let initial_velocity = 50.0;
  let initial_angle_deg = 45.0;
  println!("> Running stone flight tests...");
  println!("+---------------");

  let trajectory_no_drag = test_no_drag(initial_velocity, initial_angle_deg)?;
  println!("No drag: range ~ {:.2} m",
           trajectory_no_drag.last().map(|p| p.0).unwrap_or(0.0));
  println!("+---------------");

  let trajectory_viscous = test_viscous_drag(initial_velocity, initial_angle_deg)?;
  println!("Viscous: range ~ {:.2} m",
           trajectory_viscous.last().map(|p| p.0).unwrap_or(0.0));
  println!("+---------------");

  let trajectory_frontal = test_frontal_drag(initial_velocity, initial_angle_deg)?;
  println!("Frontal: range ~ {:.2} m",
           trajectory_frontal.last().map(|p| p.0).unwrap_or(0.0));
  println!("+---------------");

  let trajectory_curvature = test_curvature_method(initial_velocity, initial_angle_deg)?;
  println!("Curvature (no drag): range ~ {:.2} m",
           trajectory_curvature.last().map(|p| p.0).unwrap_or(0.0));
  println!("+---------------");

  let all_series: [&[(f64, f64)]; 4] = [
    &trajectory_no_drag,
    &trajectory_viscous,
    &trajectory_frontal,
    &trajectory_curvature,
  ];

  let mut max_x = 0.0f64;
  let mut max_y = 0.0f64;
  for ser in all_series.iter() {
    for (x, y) in (*ser).iter() {
      if *x > max_x {
        max_x = *x;
      }
      if *y > max_y {
        max_y = *y;
      }
    }
  }
  if max_x <= 0.0 || max_y <= 0.0 {
    return Err("nothing to plot".into());
  }
  max_x *= 1.1;
  max_y *= 1.1;

  let root_area = BitMapBackend::new(ASSETS_PATH, (800, 600)).into_drawing_area();
  root_area.fill(&WHITE)?;

  let mut chart = ChartBuilder::on(&root_area)
    .caption("Trajectory of stone flight", ("sans-serif", 30))
    .margin(10)
    .x_label_area_size(40)
    .y_label_area_size(40)
    .build_cartesian_2d(0.0..max_x, 0.0..max_y)?;

  chart
    .configure_mesh()
    .x_desc("x, m")
    .y_desc("y, m")
    .draw()?;

  chart
    .draw_series(LineSeries::new(trajectory_no_drag.clone(), &RED))?
    .label("No drag (ODE)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

  chart
    .draw_series(LineSeries::new(trajectory_viscous, &BLUE))?
    .label("Viscous drag")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

  chart
    .draw_series(LineSeries::new(trajectory_frontal, &GREEN))?
    .label("Frontal drag")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

  chart
    .draw_series(PointSeries::of_element(
      trajectory_curvature,
      2,
      &MAGENTA,
      &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
    ))?
    .label("Curvature method")
    .legend(|(x, y)| Circle::new((x + 10, y), 3, MAGENTA.filled()));

  chart
    .configure_series_labels()
    .background_style(&WHITE.mix(0.8))
    .border_style(&BLACK)
    .draw()?;

  root_area.present()?;
  println!("Plot saved to {}", ASSETS_PATH);
  println!("> Done!");

  Ok(())
}
