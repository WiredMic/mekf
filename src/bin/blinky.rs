#![deny(unsafe_code)]
#![deny(warnings)]
#![allow(unused_imports)]
#![allow(dead_code)]
#![no_main]
#![no_std]

// use blinky::geometric_algebra::ga_rotation;
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

// use blinky::geometric_algebra;
// use blinky::quarterions;

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

    defmt::debug!("Starting blink sequence!");
    loop {
        defmt::info!("Led on:");
        led.set_high();
        _delay.delay_ms(250);

        defmt::info!("Led off");
        led.set_low();
        _delay.delay_ms(250);
    }

    // All of this notation is writen with matrix algebra
    //
    // Skew symmetric matrix
    // \[ [ \vec{x} \times ] = \begin{bmatrix} 0 & -x_3 & x_2\\ x_3 & 0 & -x_1 \\ -x_2 & x_1 & 0 \end{bmatrix} \]

    // Zero vector
    // \[ \vec{0}_{3\times 1} = \begin{bmatrix} 0\\ 0\\ 0\\ \end{bmatrix}\]

    fn update_step() {
        // \[\vec{a}_{k+1|k}=\vec{0}_{3\times 1}\]
        let v_a = na::Vector3::from_element(0.0);
        assert!(v_a.x == 0.0 && v_a.y == 0.0 && v_a.z == 0.0);

        // \[\hat{\omega}_{k+1|k}=\hat{\omega}_{m,k+1}\]

        // \[P_{a,k+1|k}=[\hat{\omega}_{k+1|k}\times]P_{a,k|k}[\hat{\omega}_{k+1|k}\times]^T+\Delta t Q_{k+1}\]
        // \[\boldsymbol{q}_{k+1|k}=\boldsymbol{q}_{k|k}+\Delta t \frac{1}{2}\boldsymbol{q}_{k|k}\otimes\hat{\omega}_{k+1|k}\]
        // \[\boldsymbol{q}_{k+1|k}=\frac{\boldsymbol{q}_{k+1|k}}{||\boldsymbol{q}_{k+1|k}||}\]
    }

    fn measurement_step() {
        // \[ A(\boldsymbol{q}_{k+1|k})=I_{3\times 3} \hat{q}_4^2 -2 \hat{q}_4 [\vec{q}\times] + [\vec{q}\times]^2 + \vec{q}\,\vec{q}^{\,T} \]
        // \[h_{k+1|k}=A(\boldsymbol{q}_{k+1|k})^i \vec{z}_{i,m,k+1}\]
        // \[ ^i \vec{y}_{k+1|k}=^i \vec{z}_{z,m,k+1}-h_{k+1|k}\]
        // \[ P_{y,k+1|k}=[h_{k+1|k}\times]P_{a,k+1|k}H^{T}_{a,k+1}P^{-1}_{y,k+1|k}\]
        // \[\vec{a}_{k+1|k+1}=K_{k+1}\vec{y}_{k+1|k}\]
    }

    fn reset_step() {
        // \[\delta \boldsymbol{q}_{k+1|k+1} = \left[  \frac{\frac{1}{2}\vec{a}_{k+1|k+1}}{\sqrt{1-\frac{1}{4}\left|\vec{a}_{k+1|k+1}\right|^2}} \right]\]
        // \[\boldsymbol{q}_{k+1|k+1} = \boldsymbol{q}_{k+1|k}\otimes \delta \boldsymbol{q}_{k+1|k+1}\]
        // \[\hat{\omega}_{k+1|k+1} = \hat{\omega}_{k+1|k} + K_{k+1}\vec{y}_{k+1|k}\]
    }
}
