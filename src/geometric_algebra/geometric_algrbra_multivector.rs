// Written by a generator written by enki.
// https://bivector.net/tools.html?p=3&q=0&r=0
// Modified by Rasmus Enevoldsen

#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
// #![feature(const_slice_len)]

use core::ops::{Add, BitAnd, BitOr, BitXor, Div, Index, IndexMut, Mul, Not, Sub};
use libm::{cosf, powf, sinf, sqrtf};

// use std::f32::consts::PI;
const PI: f32 = 3.14159265358979323846;

// the 8 bases of 3D vectorspace geometric algebra
const basis: &'static [&'static str] = &["1", "e1", "e2", "e3", "e12", "e13", "e23", "e123"];
const basis_count: usize = basis.len();

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct GaMultivector {
    mvec: [f32; basis_count],
}

impl GaMultivector {
    pub const fn zero() -> Self {
        Self {
            mvec: [0.0; basis_count],
        }
    }

    pub const fn new(f: f32, idx: usize) -> Self {
        let mut ret = Self::zero();
        ret.mvec[idx] = f;
        ret
    }
}

// basis vectors are available as global constants.
const e1: GaMultivector = GaMultivector::new(1.0, 1);
const e2: GaMultivector = GaMultivector::new(1.0, 2);
const e3: GaMultivector = GaMultivector::new(1.0, 3);
const e12: GaMultivector = GaMultivector::new(1.0, 4);
const e13: GaMultivector = GaMultivector::new(1.0, 5);
const e23: GaMultivector = GaMultivector::new(1.0, 6);
const e123: GaMultivector = GaMultivector::new(1.0, 7);

impl Index<usize> for GaMultivector {
    type Output = f32;

    fn index<'a>(&'a self, index: usize) -> &'a Self::Output {
        &self.mvec[index]
    }
}

impl IndexMut<usize> for GaMultivector {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Self::Output {
        &mut self.mvec[index]
    }
}

// Reverse
// Reverse the order of the basis blades.
impl GaMultivector {
    pub fn Reverse(self: Self) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0];
        res[1] = a[1];
        res[2] = a[2];
        res[3] = a[3];
        res[4] = -a[4];
        res[5] = -a[5];
        res[6] = -a[6];
        res[7] = -a[7];
        res
    }
}

// Dual
// Poincare duality operator.
impl GaMultivector {
    pub fn Dual(self: Self) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = -a[7];
        res[1] = -a[6];
        res[2] = a[5];
        res[3] = -a[4];
        res[4] = a[3];
        res[5] = -a[2];
        res[6] = a[1];
        res[7] = a[0];
        res
    }
}

impl Not for GaMultivector {
    type Output = GaMultivector;

    fn not(self: Self) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = -a[7];
        res[1] = -a[6];
        res[2] = a[5];
        res[3] = -a[4];
        res[4] = a[3];
        res[5] = -a[2];
        res[6] = a[1];
        res[7] = a[0];
        res
    }
}

// Conjugate
// Clifford Conjugation
impl GaMultivector {
    pub fn Conjugate(self: Self) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0];
        res[1] = -a[1];
        res[2] = -a[2];
        res[3] = -a[3];
        res[4] = -a[4];
        res[5] = -a[5];
        res[6] = -a[6];
        res[7] = a[7];
        res
    }
}

// Involute
// Main involution
impl GaMultivector {
    pub fn Involute(self: Self) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0];
        res[1] = -a[1];
        res[2] = -a[2];
        res[3] = -a[3];
        res[4] = a[4];
        res[5] = a[5];
        res[6] = a[6];
        res[7] = -a[7];
        res
    }
}

// Mul
// The geometric product.
impl Mul for GaMultivector {
    type Output = GaMultivector;

    fn mul(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = b[0] * a[0] + b[1] * a[1] + b[2] * a[2] + b[3] * a[3]
            - b[4] * a[4]
            - b[5] * a[5]
            - b[6] * a[6]
            - b[7] * a[7];
        res[1] = b[1] * a[0] + b[0] * a[1] - b[4] * a[2] - b[5] * a[3] + b[2] * a[4] + b[3] * a[5]
            - b[7] * a[6]
            - b[6] * a[7];
        res[2] = b[2] * a[0] + b[4] * a[1] + b[0] * a[2] - b[6] * a[3] - b[1] * a[4]
            + b[7] * a[5]
            + b[3] * a[6]
            + b[5] * a[7];
        res[3] = b[3] * a[0] + b[5] * a[1] + b[6] * a[2] + b[0] * a[3]
            - b[7] * a[4]
            - b[1] * a[5]
            - b[2] * a[6]
            - b[4] * a[7];
        res[4] = b[4] * a[0] + b[2] * a[1] - b[1] * a[2] + b[7] * a[3] + b[0] * a[4] - b[6] * a[5]
            + b[5] * a[6]
            + b[3] * a[7];
        res[5] = b[5] * a[0] + b[3] * a[1] - b[7] * a[2] - b[1] * a[3] + b[6] * a[4] + b[0] * a[5]
            - b[4] * a[6]
            - b[2] * a[7];
        res[6] = b[6] * a[0] + b[7] * a[1] + b[3] * a[2] - b[2] * a[3] - b[5] * a[4]
            + b[4] * a[5]
            + b[0] * a[6]
            + b[1] * a[7];
        res[7] = b[7] * a[0] + b[6] * a[1] - b[5] * a[2] + b[4] * a[3] + b[3] * a[4] - b[2] * a[5]
            + b[1] * a[6]
            + b[0] * a[7];
        res
    }
}

// the norm of a multivector |A|
// \[|A|^2=\left< A\^dag A \right>_0\]
impl GaMultivector {
    pub fn Norm(self: Self) -> f32 {
        let a = self;
        let res = sqrtf((a.Reverse() * a)[0]);
        res
    }
}

// the inverse of a multivector A^-1
// \[A^{-1}=\frac{A^\dag}{|A|^2}\]
impl GaMultivector {
    pub fn Inverse(self: Self) -> GaMultivector {
        // TODO test if self is the zero multivector
        let a = self;
        a.Reverse() * (1.0 / (a.Reverse() * a)[0])
    }
}

// Wedge
// The outer product. (MEET)
impl BitXor for GaMultivector {
    type Output = GaMultivector;

    fn bitxor(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = b[0] * a[0];
        res[1] = b[1] * a[0] + b[0] * a[1];
        res[2] = b[2] * a[0] + b[0] * a[2];
        res[3] = b[3] * a[0] + b[0] * a[3];
        res[4] = b[4] * a[0] + b[2] * a[1] - b[1] * a[2] + b[0] * a[4];
        res[5] = b[5] * a[0] + b[3] * a[1] - b[1] * a[3] + b[0] * a[5];
        res[6] = b[6] * a[0] + b[3] * a[2] - b[2] * a[3] + b[0] * a[6];
        res[7] = b[7] * a[0] + b[6] * a[1] - b[5] * a[2] + b[4] * a[3] + b[3] * a[4] - b[2] * a[5]
            + b[1] * a[6]
            + b[0] * a[7];
        res
    }
}

// Vee
// The regressive product. (JOIN)
impl BitAnd for GaMultivector {
    type Output = GaMultivector;

    fn bitand(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[7] = 1.0 * (a[7] * b[7]);
        res[6] = 1.0 * (a[6] * b[7] + a[7] * b[6]);
        res[5] = -1.0 * (a[5] * -1.0 * b[7] + a[7] * b[5] * -1.0);
        res[4] = 1.0 * (a[4] * b[7] + a[7] * b[4]);
        res[3] = 1.0 * (a[3] * b[7] + a[5] * -1.0 * b[6] - a[6] * b[5] * -1.0 + a[7] * b[3]);
        res[2] = -1.0 * (a[2] * -1.0 * b[7] + a[4] * b[6] - a[6] * b[4] + a[7] * b[2] * -1.0);
        res[1] = 1.0 * (a[1] * b[7] + a[4] * b[5] * -1.0 - a[5] * -1.0 * b[4] + a[7] * b[1]);
        res[0] = 1.0
            * (a[0] * b[7] + a[1] * b[6] - a[2] * -1.0 * b[5] * -1.0 + a[3] * b[4] + a[4] * b[3]
                - a[5] * -1.0 * b[2] * -1.0
                + a[6] * b[1]
                + a[7] * b[0]);
        res
    }
}

// Dot
// The inner product.
impl BitOr for GaMultivector {
    type Output = GaMultivector;

    fn bitor(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = b[0] * a[0] + b[1] * a[1] + b[2] * a[2] + b[3] * a[3]
            - b[4] * a[4]
            - b[5] * a[5]
            - b[6] * a[6]
            - b[7] * a[7];
        res[1] = b[1] * a[0] + b[0] * a[1] - b[4] * a[2] - b[5] * a[3] + b[2] * a[4] + b[3] * a[5]
            - b[7] * a[6]
            - b[6] * a[7];
        res[2] = b[2] * a[0] + b[4] * a[1] + b[0] * a[2] - b[6] * a[3] - b[1] * a[4]
            + b[7] * a[5]
            + b[3] * a[6]
            + b[5] * a[7];
        res[3] = b[3] * a[0] + b[5] * a[1] + b[6] * a[2] + b[0] * a[3]
            - b[7] * a[4]
            - b[1] * a[5]
            - b[2] * a[6]
            - b[4] * a[7];
        res[4] = b[4] * a[0] + b[7] * a[3] + b[0] * a[4] + b[3] * a[7];
        res[5] = b[5] * a[0] - b[7] * a[2] + b[0] * a[5] - b[2] * a[7];
        res[6] = b[6] * a[0] + b[7] * a[1] + b[0] * a[6] + b[1] * a[7];
        res[7] = b[7] * a[0] + b[0] * a[7];
        res
    }
}

// Add
// Multivector addition
impl Add for GaMultivector {
    type Output = GaMultivector;

    fn add(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] + b[0];
        res[1] = a[1] + b[1];
        res[2] = a[2] + b[2];
        res[3] = a[3] + b[3];
        res[4] = a[4] + b[4];
        res[5] = a[5] + b[5];
        res[6] = a[6] + b[6];
        res[7] = a[7] + b[7];
        res
    }
}

// Sub
// Multivector subtraction
impl Sub for GaMultivector {
    type Output = GaMultivector;

    fn sub(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] - b[0];
        res[1] = a[1] - b[1];
        res[2] = a[2] - b[2];
        res[3] = a[3] - b[3];
        res[4] = a[4] - b[4];
        res[5] = a[5] - b[5];
        res[6] = a[6] - b[6];
        res[7] = a[7] - b[7];
        res
    }
}

// smul
// scalar/multivector multiplication
impl Mul<GaMultivector> for f32 {
    type Output = GaMultivector;

    fn mul(self: f32, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a * b[0];
        res[1] = a * b[1];
        res[2] = a * b[2];
        res[3] = a * b[3];
        res[4] = a * b[4];
        res[5] = a * b[5];
        res[6] = a * b[6];
        res[7] = a * b[7];
        res
    }
}

// muls
// multivector/scalar multiplication
impl Mul<f32> for GaMultivector {
    type Output = GaMultivector;

    fn mul(self: GaMultivector, b: f32) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] * b;
        res[1] = a[1] * b;
        res[2] = a[2] * b;
        res[3] = a[3] * b;
        res[4] = a[4] * b;
        res[5] = a[5] * b;
        res[6] = a[6] * b;
        res[7] = a[7] * b;
        res
    }
}

// sadd
// scalar/multivector addition
impl Add<GaMultivector> for f32 {
    type Output = GaMultivector;

    fn add(self: f32, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a + b[0];
        res[1] = b[1];
        res[2] = b[2];
        res[3] = b[3];
        res[4] = b[4];
        res[5] = b[5];
        res[6] = b[6];
        res[7] = b[7];
        res
    }
}

// adds
// multivector/scalar addition
impl Add<f32> for GaMultivector {
    type Output = GaMultivector;

    fn add(self: GaMultivector, b: f32) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] + b;
        res[1] = a[1];
        res[2] = a[2];
        res[3] = a[3];
        res[4] = a[4];
        res[5] = a[5];
        res[6] = a[6];
        res[7] = a[7];
        res
    }
}

// new ga vector type

impl GaMultivector {
    pub fn new_vector(v1: f32, v2: f32, v3: f32) -> GaMultivector {
        let mut res = GaMultivector::zero();
        res[1] = v1;
        res[2] = v2;
        res[3] = v3;
        res
    }

    pub fn new_bivector(b1: f32, b2: f32, b3: f32) -> GaMultivector {
        let mut res = GaMultivector::zero();
        res[4] = b1;
        res[5] = b2;
        res[6] = b3;
        res
    }

    // angle in radians and the 3 bivector componets
    pub fn new_rotor(angle_radians: f32, b1: f32, b2: f32, b3: f32) -> GaMultivector {
        // nomilize the bivector
        let bivector_norm = sqrtf(powf(b1, 2.0) + powf(b2, 2.0) + powf(b3, 2.0));

        // init mvec
        let mut res = GaMultivector::zero();

        res[0] = cosf(angle_radians / 2.0);
        res[4] = sinf(angle_radians / 2.0) * b1 / bivector_norm;
        res[5] = sinf(angle_radians / 2.0) * b2 / bivector_norm;
        res[6] = sinf(angle_radians / 2.0) * b3 / bivector_norm;
        res
    }
}
