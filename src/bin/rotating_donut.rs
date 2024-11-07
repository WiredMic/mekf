#![deny(unsafe_code)]
#![deny(warnings)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unused_imports)]
#![allow(dead_code)]
#![no_main]
#![no_std]

use cortex_m::delay;
use hal::timer::SysDelay;
use libm::ceilf;
//Setup serial
use rtt_target::{rprint, rprintln, rtt_init_print};

use nalgebra::Point;
//Define panic behaviour
use panic_probe as _;

//Hal library
use hal::{pac, prelude::*};
use stm32f4xx_hal as hal;

use cortex_m_rt::entry;

use blinky::geometric_algebra::geometric_algebra_bivector::GaBivector;
use blinky::geometric_algebra::geometric_algebra_multivector::GaMultivector;
use blinky::geometric_algebra::geometric_algebra_rotor::GaRotor;
use blinky::geometric_algebra::geometric_algebra_vector::GaVector;
use blinky::geometric_algebra::rotations::ga_rotation;
use core::f32::consts::PI;

#[entry] //Intro point of program
fn main() -> ! {
    rtt_init_print!();
    //Cortex core peripherals
    let cp = cortex_m::peripheral::Peripherals::take().unwrap();

    //Hal peripherals
    let dp = pac::Peripherals::take().unwrap(); //Extract peripherals

    let rcc = dp.RCC.constrain(); //Grab clocks
    let clocks = rcc.cfgr.use_hse(8.MHz()).sysclk(84.MHz()).freeze(); //Set sysclk compared to external 8Mhz oscillator:
    rprintln!("Sysclk running at: {}", clocks.sysclk().raw());

    let mut _delay = cp.SYST.delay(&clocks);

    rprint!("{}[2J", 27 as char);
    rprint!("{}[?25l", 27 as char);
    rprint!("{}[0:0H", 27 as char);
    let minor_radius: f32 = 1.0;
    let major_radius: f32 = 3.0;
    const NUM_POINT_IN_CIRC: usize = 12;
    const NUM_CIRC_IN_TORUS: usize = 12;

    let mut point_array = [GaVector::zero(); (NUM_POINT_IN_CIRC * NUM_CIRC_IN_TORUS)];
    let first_point = GaVector::new(minor_radius, 0.0, 0.0);
    point_array[0] = first_point;

    make_torus(
        &mut _delay,
        &mut point_array,
        major_radius,
        NUM_POINT_IN_CIRC,
        NUM_CIRC_IN_TORUS,
    );

    let angle3_i = PI / (256.0 * 3.0);
    let mut angle3 = angle3_i;
    let rotation_plane3 = GaBivector::new(1.0, 4.0, 2.0);
    let angle4_i = PI / (128.0 * 2.0);
    let mut angle4 = angle4_i;
    let rotation_plane4 = GaBivector::new(7.0, 4.0, -1.0);

    loop {
        rprint!("{}[2J", 27 as char);
        let mut torus = point_array;
        angle3 += angle3_i;
        if angle3 > 2.0 * PI {
            angle3 -= PI;
        }
        let rotor3 = GaRotor::new(angle3, rotation_plane3);

        angle4 += angle4_i;
        if angle4 > 2.0 * PI {
            angle4 -= PI;
        }
        let rotor4 = GaRotor::new(angle4, rotation_plane4);
        let rotor5 = rotor3 * rotor4;
        rotate_donut(&mut _delay, &mut point_array, rotor5);

        project_onto_screen(&mut _delay, &mut torus);

        // _delay.delay_ms(10);
    }
}

fn make_circle(
    delay: &mut SysDelay,
    point_array: &mut [GaVector],
    major_radius: f32,
    num_point_in_circ: usize,
) {
    // \[R_0^\dag\vec{r_0}R_0 + \vec{r}_1 \]
    // let mut sum = 0;

    // rotor to make circle
    let circ_rot_plane: GaBivector = GaBivector::new(1.0, 0.0, 0.0);
    // The angle of rotaion
    let angle1: f32 = PI * 2.0 / (num_point_in_circ as f32);
    let rotor1 = GaRotor::new(angle1, circ_rot_plane);

    for i in 1..num_point_in_circ {
        point_array[i] = ga_rotation(rotor1, point_array[i - 1]);
        // rprintln!(
        //     "the {} point in the {} circle is:, x: {}, y: {}, z: {}",
        //     i,
        //     "0",
        //     point_array[i].mvec[1],
        //     point_array[i].mvec[2],
        //     point_array[i].mvec[3],
        // );
    }

    let major_radius = GaVector::new(major_radius, 0.0, 0.0);
    for i in 0..num_point_in_circ {
        // sum += i;
        point_array[i] = point_array[i] + major_radius;
        // rprintln!(
        //     "the {} point in the {} circle is:, x: {}, y: {}, z: {}",
        //     i,
        //     "0",
        //     point_array[i].mvec[1],
        //     point_array[i].mvec[2],
        //     point_array[i].mvec[3],
        // );
        delay.delay_ms(1);
    }
}

fn make_torus(
    delay: &mut SysDelay,
    point_array: &mut [GaVector],
    major_radius: f32,
    num_point_in_circ: usize,
    num_circ_in_torus: usize,
) {
    make_circle(delay, point_array, major_radius, num_point_in_circ);

    // rotor to make torus
    let torus_rot_plane: GaBivector = GaBivector::new(0.0, 1.0, 0.0);
    // half the angle of rotation divided by the number of circles in the torus
    let angle2 = PI * 2.0 / (num_circ_in_torus as f32);
    let rotor2 = GaRotor::new(angle2, torus_rot_plane);

    let mut max = 0;

    for i in 1..num_circ_in_torus {
        // rotate circle
        // rprintln!("i: {}", i);
        delay.delay_ms(1);
        for j in 0..num_point_in_circ {
            // num_points to excluding num_points * num_circ
            let index = i * num_point_in_circ + j;
            delay.delay_ms(1);
            point_array[i * num_circ_in_torus + j] =
                ga_rotation(rotor2, point_array[(i - 1) * num_circ_in_torus + j]);

            // rprintln!(
            //     "athe {} point in the {} circle is:, x: {}, y: {}, z: {}",
            //     j,
            //     i,
            //     point_array[i * num_circ_in_torus + j].mvec[1],
            //     point_array[i * num_circ_in_torus + j].mvec[2],
            //     point_array[i * num_circ_in_torus + j].mvec[3],
            // );

            // delay.delay_ms(10);
            // rprintln!("  j: {}", j);
            // rprintln!("  index: {}", index);
            max = index;
        }
    }
    rprintln!("{}", max);
}

fn rotate_donut(delay: &mut SysDelay, point_array: &mut [GaVector], rotor5: GaRotor) {
    // rprintln!("The rotor is");
    // rprintln!("    Scalar: {}", rotor5[0]);
    // rprintln!("    e1e2: {}", rotor5[4]);
    // rprintln!("    e3e1: {}", rotor5[5]);
    // rprintln!("    e23: {}", rotor5[6]);

    for i in 0..point_array.len() {
        point_array[i] = ga_rotation(rotor5, point_array[i]);
        delay.delay_ms(1);

        // rprintln!(
        //     "index: {} of the rotated torus is:, x: {}, y: {}, z: {}",
        //     i,
        //     point_array[i].mvec[1],
        //     point_array[i].mvec[2],
        //     point_array[i].mvec[3],
        // );
    }
}

fn project_onto_screen(delay: &mut SysDelay, point_array: &mut [GaVector]) {
    // \[\vec{v}_\parallel=(\vec{v}\cdot\vec{B})\vec{B}^{-1} \]
    use libm::floorf;
    let screen_plane = GaBivector::new(0.0, 1.0, 0.0);
    let screen_plane_inv = screen_plane.Inverse();

    for i in 0..point_array.len() {
        point_array[i] = blinky::geometric_algebra::geometric_algebra_vector::GaVector {
            mvec: (point_array[i] | screen_plane) * screen_plane_inv,
        };
        delay.delay_ms(1);
        // iterate over all element. Scale them up so that the 4 * 1 torus fits on a 100 * 100 grit
        // 25.0/4.0
        // x = x + 25
        // y = 25 - z
        let x = (ceilf(point_array[i].mvec[1] * 25.0 / 4.0) as i32) + 35;
        let y = 15 - (ceilf(point_array[i].mvec[3] * 15.0 / 4.0) as i32);

        // rprintln!("index: {}, x: {}, y: {}", i, x, y);
        rprint!("{}[{};{}H", 27 as char, y, x);
        rprint!("#");
        // rprintln!("{}", i);
    }
}
