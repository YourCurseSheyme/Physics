//! Testing. One task - one module

mod stone_flight;
mod billiard;
// mod fly_me_to_the_mars;

use stone_flight::run_stone_tests;
use billiard::run_billiard_tests;
// use fly_me_to_the_mars::run_fly_to_the_mars_tests;

fn main() {
    run_stone_tests().expect("Error while plotting");
    run_billiard_tests().expect("Error while simulating billiard");
    // run_fly_to_the_mars_tests().expect("Error while simulating flight");
}
