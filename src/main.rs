//! Testing. One task - one module

mod stone_flight;  // M1
use stone_flight::run_stone_tests;

fn main() {
    run_stone_tests().expect("Error while plotting");
}
