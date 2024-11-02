#![deny(unsafe_code)]
#![deny(warnings)]
#![allow(unused_imports)]
#![allow(dead_code)]
#![no_main]
#![no_std]

//Setup serial
use defmt as _;
use defmt_brtt as _;

//Define panic behaviour
use panic_probe as _;

//Hal library
use hal::{pac, prelude::*};
use stm32f4xx_hal as hal;

use cortex_m_rt::entry;

//Matrix Library
use libm::powf;
use nalgebra as na;

use blinky::geometric_algebra::geometric_algebra_bivector::GaBivector;
use blinky::geometric_algebra::geometric_algebra_multivector::GaMultivector;
use blinky::geometric_algebra::geometric_algebra_rotor::GaRotor;
use blinky::geometric_algebra::geometric_algebra_vector::GaVector;
use blinky::geometric_algebra::rotations::ga_rotation;
use core::f32::consts::PI;

#[entry] //Intro point of program
fn main() -> ! {
    //Cortex core peripherals
    let cp = cortex_m::peripheral::Peripherals::take().unwrap();

    //Hal peripherals
    let dp = pac::Peripherals::take().unwrap(); //Extract peripherals

    let rcc = dp.RCC.constrain(); //Grab clocks
    let clocks = rcc.cfgr.use_hse(8.MHz()).sysclk(84.MHz()).freeze(); //Set sysclk compared to external 8Mhz oscillator:
    defmt::info!("Sysclk running at: {}", clocks.sysclk().raw());

    let mut _delay = cp.SYST.delay(&clocks);

    let gpioa = dp.GPIOA.split(); //Extract GPIOA
    let mut led = gpioa.pa5.into_push_pull_output(); //Set pin A5 as output

    // https://docs.rs/stm32f4xx-hal/latest/stm32f4xx_hal/i2c/struct.I2c.html
    // https://electropeak.com/learn/interfacing-gy-91-9-axis-mpu9250-bmp280-module-with-arduino/
    // https://docs.rs/mpu9250/latest/mpu9250/
    // i2c
    // hal::i2c::

    let mut vector = GaVector::new(3.0, 0.0, 0.0);

    let angle: f32 = PI / 2.0;
    let bivector = GaBivector::new(1.0, 0.0, 0.0);
    let rotor = GaRotor::new(angle, bivector);

    defmt::debug!("Starting blink sequence!");
    loop {
        // defmt::info!("Led on:");
        led.set_high();
        // _delay.delay_ms(250);

        // defmt::info!("Led off");
        // led.set_low();
        // _delay.delay_ms(250);

        defmt::info!("The input vector is {}", vector);
        let vector_rot = ga_rotation(rotor, vector);
        defmt::info!("The rotated vector is {}", vector_rot);
        vector = vector_rot;
    }
}
