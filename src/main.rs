//! Testing. One task - one module

mod stone_flight;
mod billiard;

use stone_flight::run_stone_tests;
use billiard::run_billiard_tests;

// mod fly_me_to_the_Mars;

fn main() {
    run_stone_tests().expect("Error while plotting");
    run_billiard_tests().expect("Error while simulating billiard");

}
