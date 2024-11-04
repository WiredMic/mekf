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
extern crate defmt;
use core::f32::consts::PI;
use defmt::Format;

// the 8 bases of 3D vectorspace geometric algebra
const basis: &'static [&'static str] = &["1", "e1", "e2", "e3", "e12", "e31", "e23", "e123"];
pub const basis_count: usize = basis.len();

#[derive(Default, Debug, Clone, Copy, PartialEq, Format)]
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

    pub const fn new_mvec(
        scalar: f32,
        e1: f32,
        e2: f32,
        e3: f32,
        e12: f32,
        e31: f32,
        e23: f32,
        e123: f32,
    ) -> Self {
        let mut ret = Self::zero();
        ret.mvec[0] = scalar;
        ret.mvec[1] = e1;
        ret.mvec[2] = e2;
        ret.mvec[3] = e3;
        ret.mvec[4] = e12;
        ret.mvec[5] = e31;
        ret.mvec[6] = e23;
        ret.mvec[7] = e123;
        ret
    }
}

// basis vectors are available as global constants.
// const e1: GaMultivector = GaMultivector::new(1.0, 1);
// const e2: GaMultivector = GaMultivector::new(1.0, 2);
// const e3: GaMultivector = GaMultivector::new(1.0, 3);
// const e12: GaMultivector = GaMultivector::new(1.0, 4);
// const e23: GaMultivector = GaMultivector::new(1.0, 5);
// const e31: GaMultivector = GaMultivector::new(1.0, 6);
// const e123: GaMultivector = GaMultivector::new(1.0, 7);

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
        res[0] = -a[7]; // e123e123 = -scalar
        res[1] = -a[6]; // e23e123 = -e1
        res[2] = -a[5]; // e31e123 = -e2
        res[3] = -a[4]; // e12e123 = -e3
        res[4] = a[3]; // e3e123 = e12
        res[5] = a[2]; // e2e123 = e31
        res[6] = a[1]; // e1e123 = e23
        res[7] = a[0]; // scalar e123
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
        res[2] = -a[5];
        res[3] = -a[4];
        res[4] = a[3];
        res[5] = a[2];
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
        res[0] = a[0]*b[0] // scalar/scalar
            +a[1]*b[1]+a[2]*b[2]+a[3]*b[3] // vector dot product
            -a[4]*b[4]-a[5]*b[5]-a[6]*b[6] // bivector dot produt
            -a[7]*b[7]; // trivector/trivector
        res[1] = (a[0]*b[1]+a[1]*b[0]) // scalar/e1
            +(a[4]*b[2]-a[2]*b[4]) // e12e2
            +(a[3]*b[5]-a[5]*b[3]) // e3e13
            +(-a[6]*b[7]-a[7]*b[6]); // e23e123
        res[2] = (a[0]*b[2]+a[2]*b[0]) // scalar/e2
            +(a[1]*b[4]-a[4]*b[1]) // e1e12
            +(a[6]*b[3]-a[3]*b[6]) // e23e3
            +(-a[5]*b[7]-a[7]*b[5]); // e31e123
        res[3] = (a[0]*b[3]+a[3]*b[0]) // scalar/e3
            +(a[5]*b[1]-a[1]*b[5]) // e31e1
            +(a[2]*b[6]-a[6]*b[2]) // e2e23
            +(-a[4]*b[7]-a[7]*b[4]); //e12e123
        res[4] = (a[0]*b[4]+a[4]*b[0]) // scalar/e12
            +(a[1]*b[2]-a[2]*b[1]) // e1e2
            +(a[5]*b[6]-a[6]*b[5]) // e31e23
            +(a[3]*b[7]+a[7]*b[3]); // e3e123
        res[5] = (a[0]*b[5]+a[5]*b[0]) // scalar/e31
            +(a[3]*b[1]-a[1]*b[3]) // e3e1
            +(a[6]*b[4]-a[4]*b[6]) // e23e12
            +(a[2]*b[7]+a[7]*b[2]); // e2e123
        res[6] = (a[0]*b[6]+a[6]*b[0]) // scalar/e23
            +(a[2]*b[3]-a[3]*b[2]) // e2e3
            +(a[4]*b[5]-a[5]*b[4]) // e12e31
            +(a[1]*b[7]+a[7]*b[1]); //e1e123
        res[7] = a[7] * b[0]
            + a[6] * b[1]
            + a[5] * b[2]
            + a[4] * b[3]
            + a[3] * b[4]
            + a[2] * b[5]
            + a[1] * b[6]
            + a[0] * b[7];
        res
    }
}

#[cfg(test)]
mod mvec_mvec_mul {
    // https://bivector.net/tools.html?p=3&q=0&r=0
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn mvec_mvec_mul() {
        let mvec1 = GaMultivector::new_mvec(6.0, 9.0, 7.0, 4.0, 7.0, 4.0, 8.0, 7.0);
        let mvec2 = GaMultivector::new_mvec(5.0, 8.0, 7.0, 3.0, 2.0, 8.0, 2.0, 1.0);
        let mvec_res = mvec1 * mvec2;
        // 94+126e1​−5e2​−65e3​+23e12​−131e13​+158e23​+236e123
        assert_relative_eq!(mvec_res[0], 94.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], 126.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], -5.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], -65.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 23.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 131.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 158.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 236.0, max_relative = 0.000001);
    }
    #[test]
    fn negetive_mvec_mvec_mul() {
        // let mvec1 = GaMultivector::new_mvec(-6.0, -8.0, -4.0, -1.0, -6.0, -4.0, -8.0, -5.0);
        let mvec1 = GaMultivector::new_mvec(-4.0, -1.0, -3.0, -2.0, -9.0, -6.0, -3.0, -10.0);
        let mvec2 = GaMultivector::new_mvec(-4.0, -2.0, -4.0, -9.0, -2.0, -1.0, -7.0, -1.0);
        let mvec_res = mvec1 * mvec2;
        // −7−83e1​+9e2​+35e3​+173e12​−9e13​+77e23​+169e123
        assert_relative_eq!(mvec_res[0], -7.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], -83.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], 9.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], 35.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 173.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 9.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 77.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 169.0, max_relative = 0.000001);
    }
}

// Wedge
// The outer product. (MEET)
impl BitXor for GaMultivector {
    type Output = GaMultivector;

    fn bitxor(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] * b[0]; // scalar/scalar
                              // +a[1]*b[1]+a[2]*b[2]+a[3]*b[3] // vector dot product
                              // -a[4]*b[4]-a[5]*b[5]-a[6]*b[6] // bivector dot produt
                              // -a[7]*b[7]; // trivector/trivector
        res[1] = a[0] * b[1] + a[1] * b[0]; // scalar/e1
                                            // +(a[4]*b[2]-a[2]*b[4]) // e12e2
                                            // +(a[3]*b[5]-a[5]*b[3]) // e3e13
                                            // +(-a[6]*b[7]-a[7]*b[6]); // e23e123
        res[2] = a[0] * b[2] + a[2] * b[0]; // scalar/e2
                                            // +(a[1]*b[4]-a[4]*b[1]) // e1e12
                                            // +(a[6]*b[3]-a[3]*b[6]) // e23e3
                                            // +(-a[5]*b[7]-a[7]*b[5]); // e31e123
        res[3] = a[0] * b[3] + a[3] * b[0]; // scalar/e3
                                            // +(a[5]*b[1]-a[1]*b[5]) // e31e1
                                            // +(a[2]*b[6]-a[6]*b[2]) // e2e23
                                            // +(-a[4]*b[7]-a[7]*b[4]); //e12e123
        res[4] = (a[0]*b[4]+a[4]*b[0]) // scalar/e12
            +(a[1]*b[2]-a[2]*b[1]); // e1e2
                                    // +(a[5]*b[6]-a[6]*b[5]) // e31e23
                                    // +(a[3]*b[7]+a[7]*b[3]); // e3e123
        res[5] = (a[0]*b[5]+a[5]*b[0]) // scalar/e31
            +(a[3]*b[1]-a[1]*b[3]); // e3e1
                                    // +(a[6]*b[4]-a[4]*b[6]) // e23e12
                                    // +(a[2]*b[7]+a[7]*b[2]); // e2e123
        res[6] = (a[0]*b[6]+a[6]*b[0]) // scalar/e23
            +(a[2]*b[3]-a[3]*b[2]); // e2e3
                                    // +(a[4]*b[5]-a[5]*b[4]) // e12e31
                                    // +(a[1]*b[7]+a[7]*b[1]); //e1e123
        res[7] = a[7] * b[0]
            + a[6] * b[1]
            + a[5] * b[2]
            + a[4] * b[3]
            + a[3] * b[4]
            + a[2] * b[5]
            + a[1] * b[6]
            + a[0] * b[7];
        res
    }
}

// scalar/rotor wedge
impl BitXor<GaMultivector> for f32 {
    type Output = GaMultivector;

    fn bitxor(self: f32, b: GaMultivector) -> GaMultivector {
        let mut a = GaMultivector::zero();
        a[0] = self;
        a ^ b
    }
}

// rotor/scalar multiplication
impl BitXor<f32> for GaMultivector {
    type Output = GaMultivector;

    fn bitxor(self: GaMultivector, scalar: f32) -> GaMultivector {
        let mut b = GaMultivector::zero();
        b[0] = scalar;
        self ^ b
    }
}

#[cfg(test)]
mod mvec_wedge {
    // https://bivector.net/tools.html?p=3&q=0&r=0
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn mvec_mvec_wedge() {
        // 6.0+ 9.0e1 +7.0e2+ 4.0e3+ 7.0e12 + 4.0e31 +8.0e23+ 7.0e123
        let mvec1 = GaMultivector::new_mvec(6.0, 9.0, 7.0, 4.0, 7.0, 4.0, 8.0, 7.0);
        // 5.0+ 8.0e1+ 7.0e2+ 3.0e3+ 2.0e12+ 8.0e31+ 2.0e23+ 1.0e123
        let mvec2 = GaMultivector::new_mvec(5.0, 8.0, 7.0, 3.0, 2.0, 8.0, 2.0, 1.0);
        let mvec_res = mvec1 ^ mvec2;
        // 30+93e1​+77e2​+38e3​+54e12​−73e13​+45e23​+236e123
        assert_relative_eq!(mvec_res[0], 30.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], 93.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], 77.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], 38.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 54.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 73.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 45.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 236.0, max_relative = 0.000001);
    }
    #[test]
    fn negetive_mvec_mvec_wedge() {
        // -4.0-1.0e1 -3.0e2 -2.0e3 -9.0e12 -6.0e31 -3.0e23 -10.0e123
        let mvec1 = GaMultivector::new_mvec(-4.0, -1.0, -3.0, -2.0, -9.0, -6.0, -3.0, -10.0);
        // -4.0 -2.0e1 -4.0e2 -9.0e3 -2.0e12 -1.0e31 -7.0e23 -1.0e123
        let mvec2 = GaMultivector::new_mvec(-4.0, -2.0, -4.0, -9.0, -2.0, -1.0, -7.0, -1.0);
        let mvec_res = mvec1 ^ mvec2;
        // 16+12e1​+28e2​+44e3​+42e12​−23e13​+59e23​+169e123
        assert_relative_eq!(mvec_res[0], 16.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], 12.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], 28.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], 44.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 42.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 23.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 59.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 169.0, max_relative = 0.000001);
    }
}

// Vee
// \[ A \vee B = ( A \mathrm{I} \wedge B \mathrm{I} ) \mathrm{I}^{-1} \]
// \[\mathrm{I}^{-1}=-\mathrm{I}\]
// The regressive product. (JOIN)
impl BitAnd for GaMultivector {
    type Output = GaMultivector;

    fn bitand(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        -1.0 * (self.Dual() ^ b.Dual()).Dual()
    }
}

// Dot
// The inner product.
// Grade lowering by the abs of the two grades subtracted
// \[A \cdot B = \left< AB \right>_{|a-b|}\]
impl BitOr for GaMultivector {
    type Output = GaMultivector;

    fn bitor(self: GaMultivector, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;

        res[0] = a[0]*b[0] // scalar/scalar
            +a[1]*b[1]+a[2]*b[2]+a[3]*b[3] // vector dot product
            -a[4]*b[4]-a[5]*b[5]-a[6]*b[6] // bivector dot produt
            -a[7]*b[7]; // trivector/trivector
        res[1] = (a[0]*b[1]+a[1]*b[0]) // scalar/e1
            +(a[4]*b[2]-a[2]*b[4]) // e12e2
            +(a[3]*b[5]-a[5]*b[3]) // e3e13
            +(-a[6]*b[7]-a[7]*b[6]); // e23e123
        res[2] = (a[0]*b[2]+a[2]*b[0]) // scalar/e2
            +(a[1]*b[4]-a[4]*b[1]) // e1e12
            +(a[6]*b[3]-a[3]*b[6]) // e23e3
            +(-a[5]*b[7]-a[7]*b[5]); // e31e123
        res[3] = (a[0]*b[3]+a[3]*b[0]) // scalar/e3
            +(a[5]*b[1]-a[1]*b[5]) // e31e1
            +(a[2]*b[6]-a[6]*b[2]) // e2e23
            +(-a[4]*b[7]-a[7]*b[4]); //e12e123
        res[4] = (a[0]*b[4]+a[4]*b[0]) // scalar/e12
            // +(a[1]*b[2]-a[2]*b[1]) // e1e2
            // +(a[5]*b[6]-a[6]*b[5]) // e31e23
            +(a[3]*b[7]+a[7]*b[3]); // e3e123
        res[5] = (a[0]*b[5]+a[5]*b[0]) // scalar/e31
            // +(a[3]*b[1]-a[1]*b[3]) // e3e1
            // +(a[6]*b[4]-a[4]*b[6]) // e23e12
            +(a[2]*b[7]+a[7]*b[2]); // e2e123
        res[6] = (a[0]*b[6]+a[6]*b[0]) // scalar/e23
            // +(a[2]*b[3]-a[3]*b[2]) // e2e3
            // +(a[4]*b[5]-a[5]*b[4]) // e12e31
            +(a[1]*b[7]+a[7]*b[1]); //e1e123
        res[7] = a[7] * b[0]
            // + a[6] * b[1]
            // + a[5] * b[2]
            // + a[4] * b[3]
            // + a[3] * b[4]
            // + a[2] * b[5]
            // + a[1] * b[6]
            + a[0] * b[7];
        res
    }
}

#[cfg(test)]
mod mvec_dot {
    // https://bivector.net/tools.html?p=3&q=0&r=0
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn mvec_mvec_dot() {
        // 6.0+ 9.0e1 +7.0e2+ 4.0e3+ 7.0e12 + 4.0e31 +8.0e23+ 7.0e123
        let mvec1 = GaMultivector::new_mvec(6.0, 9.0, 7.0, 4.0, 7.0, 4.0, 8.0, 7.0);
        // 5.0+ 8.0e1+ 7.0e2+ 3.0e3+ 2.0e12+ 8.0e31+ 2.0e23+ 1.0e123
        let mvec2 = GaMultivector::new_mvec(5.0, 8.0, 7.0, 3.0, 2.0, 8.0, 2.0, 1.0);
        let mvec_res = mvec1 | mvec2;
        // 94+126e1​−5e2​−65e3​+72e12​−124e13​+117e23​+41e123
        assert_relative_eq!(mvec_res[0], 94.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], 126.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], -5.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], -65.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 72.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 124.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 117.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 41.0, max_relative = 0.000001);
    }
    #[test]
    fn negetive_mvec_mvec_dot() {
        // let mvec1 = GaMultivector::new_mvec(-6.0, -8.0, -4.0, -1.0, -6.0, -4.0, -8.0, -5.0);
        let mvec1 = GaMultivector::new_mvec(-4.0, -1.0, -3.0, -2.0, -9.0, -6.0, -3.0, -10.0);
        let mvec2 = GaMultivector::new_mvec(-4.0, -2.0, -4.0, -9.0, -2.0, -1.0, -7.0, -1.0);
        let mvec_res = mvec1 | mvec2;
        // −7−83e1​+9e2​+35e3​+136e12​−71e13​+61e23​+44e123
        assert_relative_eq!(mvec_res[0], -7.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[1], -83.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[2], 9.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[3], 35.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[4], 136.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[5], 71.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[6], 61.0, max_relative = 0.000001);
        assert_relative_eq!(mvec_res[7], 44.0, max_relative = 0.000001);
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

// scalar/multivector subtraction
impl Sub<GaMultivector> for f32 {
    type Output = GaMultivector;

    fn sub(self: f32, b: GaMultivector) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a - b[0];
        res[1] = -b[1];
        res[2] = -b[2];
        res[3] = -b[3];
        res[4] = -b[4];
        res[5] = -b[5];
        res[6] = -b[6];
        res[7] = -b[7];
        res
    }
}

// multivector/scalar subtraction
impl Sub<f32> for GaMultivector {
    type Output = GaMultivector;

    fn sub(self: GaMultivector, b: f32) -> GaMultivector {
        let mut res = GaMultivector::zero();
        let a = self;
        res[0] = a[0] - b;
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
