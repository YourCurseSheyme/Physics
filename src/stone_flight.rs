//! First task (M1)

use ode_solvers::*;
use plotters::prelude::*;
use ode_solvers::{Vector2, Vector4};

const ASSETS_PATH: &str = "/home/sheyme/workspace/physics/plots/stone_flight.png";

#[derive(Debug, Clone, Copy)]
pub enum DragModel {
  NoDrag,
  Viscous(f64),
  Frontal(f64),
}

#[derive(Debug, Clone, Copy)]
pub struct SimulationParameters {
  pub mass: f64,
  pub g: f64,
  pub drag_model: DragModel,
}

struct ProjectileSystem {
  params: SimulationParameters,
}

impl System<f64, Vector4<f64>> for ProjectileSystem {
  fn system(&self, _t: f64, y: &Vector4<f64>, dy: &mut Vector4<f64>) {
    let velocity = Vector2::new(y[2], y[3]);
    let gravity_force = Vector2::new(0.0, -self.params.mass * self.params.g);
    let drag_force = match self.params.drag_model {
      DragModel::NoDrag => Vector2::zeros(),
      DragModel::Viscous(k) => -k * velocity,
      DragModel::Frontal(k) => -k * velocity.norm() * velocity,
    };
    let total_force = gravity_force + drag_force;
    let acceleration = total_force / self.params.mass;

    dy[0] = velocity.x;
    dy[1] = velocity.y;
    dy[2] = acceleration.x;
    dy[3] = acceleration.y;
  }

  fn solout(&mut self, _x: f64, y: &Vector4<f64>, _dy: &Vector4<f64>) -> bool {
    y[1] < 0.0
  }
}

pub struct StoneFlightSimulator {
  initial_velocity: f64,
  initial_angle_deg: f64,
  params: SimulationParameters,
}

impl StoneFlightSimulator {
  pub fn new(
    initial_velocity: f64,
    initial_angle_deg: f64,
    params: SimulationParameters
  ) -> Self {
    Self {
      initial_velocity, initial_angle_deg, params,
    }
  }

  pub fn calculate_trajectory(&self) -> Vec<(f64, f64)> {
    let angle_rad = self.initial_angle_deg.to_radians();
    let v0_x = self.initial_velocity * angle_rad.cos();
    let v0_y = self.initial_velocity * angle_rad.sin();

    let y0 = Vector4::new(0.0, 0.0, v0_x, v0_y);
    let system = ProjectileSystem { params: self.params };
    let mut stepper = Dopri5::new(system, 0.0, 100.0, 0.1, y0, 1.0e-6, 1.0e-6);

    let res = stepper.integrate();
    match res {
      Ok(stats) => {
        println!("{}", stats);
        let states = stepper.y_out();
        states.iter().map(|s| (s[0], s[1])).collect()
      }
      Err(e) => {
        println!("Error while integrating: {}", e);
        Vec::new()
      }
    }
  }
}

fn test_no_drag(initial_velocity: f64, initial_angle_deg: f64) -> Vec<(f64, f64)> {
  println!("First test without drag");
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::NoDrag,
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  let trajectory = simulator.calculate_trajectory();
  let landing_point = trajectory.last().unwrap_or(&(0.0, 0.0));
  println!("Landing point: {:.2} m", landing_point.0);
  trajectory
}

fn test_viscous_drag(initial_velocity: f64, initial_angle_deg: f64) -> Vec<(f64, f64)> {
  println!("Second test viscous drag");
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::Viscous(0.1),
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  let trajectory = simulator.calculate_trajectory();
  let landing_point = trajectory.last().unwrap_or(&(0.0, 0.0));
  println!("Landing point: {:.2} m", landing_point.0);
  trajectory
}

fn test_frontal_drag(initial_velocity: f64, initial_angle_deg: f64) -> Vec<(f64, f64)> {
  println!("Third test frontal drag");
  let params = SimulationParameters {
    mass: 1.0,
    g: 9.81,
    drag_model: DragModel::Frontal(0.01),
  };
  let simulator = StoneFlightSimulator::new(initial_velocity, initial_angle_deg, params);
  let trajectory = simulator.calculate_trajectory();
  let landing_point = trajectory.last().unwrap_or(&(0.0, 0.0));
  println!("Landing point: {:.2} m", landing_point.0);
  trajectory
}

pub fn run_stone_tests() -> Result<(), Box<dyn std::error::Error>> {
  let initial_velocity = 50.0;
  let initial_angle_deg = 45.0;

  let trajectory_no_drag = test_no_drag(initial_velocity, initial_angle_deg);
  println!("+---------------");
  let trajectory_viscous = test_viscous_drag(initial_velocity, initial_angle_deg);
  println!("+---------------");
  let trajectory_frontal = test_frontal_drag(initial_velocity, initial_angle_deg);
  println!("+---------------");

  let root_area = BitMapBackend::new(ASSETS_PATH, (800, 600)).into_drawing_area();
  root_area.fill(&WHITE)?;

  let landing_point_no_drag = trajectory_no_drag.last().unwrap_or(&(0.0, 0.0));
  let max_x = landing_point_no_drag.0 * 1.1;
  let max_y = trajectory_no_drag.iter().map(|p| p.1).fold(0.0, f64::max) * 1.1;

  let mut chart = ChartBuilder::on(&root_area)
    .caption("Trajectory of stone flight", ("sans-serif", 30))
    .margin(10)
    .x_label_area_size(40)
    .y_label_area_size(40)
    .build_cartesian_2d(0.0..max_x, 0.0..max_y)?;

  chart.configure_mesh().draw()?;

  chart.draw_series(LineSeries::new(trajectory_no_drag, &RED))?
    .label("No drag").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

  chart.draw_series(LineSeries::new(trajectory_viscous, &BLUE))?
    .label("Viscous drag").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

  chart.draw_series(LineSeries::new(trajectory_frontal, &GREEN))?
    .label("Frontal drag").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

  chart.configure_series_labels()
    .background_style(&WHITE.mix(0.8))
    .border_style(&BLACK)
    .draw()?;

  root_area.present()?;
  println!("Plot have been saved to {}", ASSETS_PATH);
  println!("Tests have been done!");

  Ok(())
}
