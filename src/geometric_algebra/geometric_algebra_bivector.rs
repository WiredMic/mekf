#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
// #![feature(const_slice_len)]

use core::ops::{Add, BitAnd, BitOr, BitXor, Div, Index, IndexMut, Mul, Not, Sub};
use libm::{cosf, powf, sinf, sqrtf};

use super::geometric_algebra_multivector::basis_count;
use super::geometric_algebra_multivector::GaMultivector;
use super::geometric_algebra_vector::GaVector;

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct GaBivector {
    pub mvec: GaMultivector,
}

impl GaBivector {
    pub fn zero() -> Self {
        Self {
            mvec: GaMultivector::zero(),
        }
    }

    pub fn new(e1e2: f32, e1e3: f32, e2e3: f32) -> Self {
        let mut ret = Self::zero();
        ret.mvec[4] = e1e2;
        ret.mvec[5] = e1e3;
        ret.mvec[6] = e2e3;
        ret
    }

    pub fn new_from_vec(vector1: GaVector, vector2: GaVector) -> GaBivector {
        vector1 ^ vector2
    }
}

impl Index<usize> for GaBivector {
    type Output = f32;

    fn index<'a>(&'a self, index: usize) -> &'a Self::Output {
        &self.mvec[index]
    }
}

impl IndexMut<usize> for GaBivector {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Self::Output {
        &mut self.mvec[index]
    }
}

// Reverse
// Reverse the order of the basis blades.
impl GaBivector {
    pub fn Reverse(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Reverse(),
        }
    }
}

// Dual
// \[\vec{B}\mathrm{I}=-\vec{v}\]
// Poincare duality operator.
impl GaBivector {
    pub fn Dual(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Dual(),
        }
    }
}

impl Not for GaBivector {
    type Output = GaBivector;

    fn not(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector { mvec: !self.mvec }
    }
}

// Conjugate
// Clifford Conjugation
impl GaBivector {
    pub fn Conjugate(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Conjugate(),
        }
    }
}

// Involute
// Main involution
impl GaBivector {
    pub fn Involute(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Involute(),
        }
    }
}

// Mul
// \[ \vec{B}\vec{C} = \left<\vec{B}\vec{C}\right>_{2-2} \]
// The geometric product.
impl Mul for GaBivector {
    type Output = f32;

    fn mul(self: GaBivector, b: GaBivector) -> f32 {
        (self.mvec * b.mvec)[0]
    }
}

// scalar/bivector multiplication
impl Mul<GaBivector> for f32 {
    type Output = GaBivector;

    fn mul(self: f32, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self * b.mvec,
        }
    }
}

// bivector/scalar multiplication
impl Mul<f32> for GaBivector {
    type Output = GaBivector;

    fn mul(self: GaBivector, b: f32) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec * b,
        }
    }
}

// bivector/multivector multiplication
impl Mul<GaMultivector> for GaBivector {
    type Output = GaMultivector;

    fn mul(self: GaBivector, b: GaMultivector) -> GaMultivector {
        self.mvec * b
    }
}

// multivector/bivector multiplication
impl Mul<GaBivector> for GaMultivector {
    type Output = GaMultivector;

    fn mul(self: GaMultivector, b: GaBivector) -> GaMultivector {
        self * b.mvec
    }
}

// bivector/vector multiplication
impl Mul<GaVector> for GaBivector {
    type Output = GaMultivector;

    fn mul(self: GaBivector, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector multiplication
impl Mul<GaBivector> for GaVector {
    type Output = GaMultivector;

    fn mul(self: GaVector, b: GaBivector) -> GaMultivector {
        self * b.mvec
    }
}

#[cfg(test)]
mod bivector_mul {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn bivector_bivector_mul() {
        let bivector1 = GaBivector::new(3.0, 5.0, 4.0);
        let bivector2 = GaBivector::new(2.0, 1.0, 6.0);
        let scalar = bivector1 * bivector2;
        assert_relative_eq!(scalar, -35.0, max_relative = 0.000001);
    }

    // #[test]
    // fn vector_mvec_wedge() {}

    // #[test]
    // fn mvec_vector_wedge() {}
}

// Wedge
// \[ \vec{B} \wedge \vec{C} = \left<\vec{B}\vec{C}\right>_{2+2} \]
// The outer product. (MEET)
impl BitXor for GaBivector {
    type Output = f32;

    fn bitxor(self: GaBivector, _b: GaBivector) -> f32 {
        0.0
    }
    // type Output = GaBivector;

    // fn bitxor(self: GaBivector, b: GaBivector) -> GaBivector {
    //     super::geometric_algebra_bivector::GaBivector {
    //         mvec: self.mvec ^ b.mvec,
    //     }
    // }
}

// bivector/multivector wedge
// \[ \vec{B}\wedge M \]
impl BitXor<GaMultivector> for GaBivector {
    type Output = GaMultivector;

    fn bitxor(self: GaBivector, b: GaMultivector) -> GaMultivector {
        self.mvec ^ b
    }
}

// multivector/bivector wedge
// \[ M\wedge \vec{B} \]
impl BitXor<GaBivector> for GaMultivector {
    type Output = GaMultivector;

    fn bitxor(self: GaMultivector, b: GaBivector) -> GaMultivector {
        self ^ b.mvec
    }
}

// bivector/bivector wedge
// \[ \vec{B}\wedge \vec{v} \]
impl BitXor<GaVector> for GaBivector {
    type Output = GaMultivector;

    fn bitxor(self: GaBivector, b: GaVector) -> GaMultivector {
        self.mvec ^ b
    }
}

// vector/bivector wedge
// \[ \vec{v} \wedge \vec{B} \]
impl BitXor<GaBivector> for GaVector {
    type Output = GaMultivector;

    fn bitxor(self: GaVector, b: GaBivector) -> GaMultivector {
        self ^ b.mvec
    }
}

#[cfg(test)]
mod bivector_wedge {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn bivector_vector_wedge() {
        let bivector = GaBivector::new(3.0, 5.0, 4.0);
        let vector = GaVector::new(2.0, 1.0, 6.0);
        let trivector = bivector ^ vector;
        assert_relative_eq!(trivector[7], 21.0, max_relative = 0.000001);
    }

    // #[test]
    // fn vector_bivector_wedge() {}
}

// Vee
// The regressive product. (JOIN)
impl BitAnd for GaBivector {
    type Output = GaBivector;

    fn bitand(self: GaBivector, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec & b.mvec,
        }
    }
}

// Dot
// The inner product.
impl BitOr for GaBivector {
    type Output = f32;

    fn bitor(self: GaBivector, b: GaBivector) -> f32 {
        (self.mvec | b.mvec)[0]
    }
}

// bivector/multivector dot
// \[ \vec{B}\cdot M \]
impl BitOr<GaMultivector> for GaBivector {
    type Output = GaMultivector;

    fn bitor(self: GaBivector, b: GaMultivector) -> GaMultivector {
        self.mvec | b
    }
}

// multivector/bivector dot
// \[ M \cdot \vec{B} \]
impl BitOr<GaBivector> for GaMultivector {
    type Output = GaMultivector;

    fn bitor(self: GaMultivector, b: GaBivector) -> GaMultivector {
        self | b.mvec
    }
}

// bivector/vector dot
// \[ \vec{B}\cdot \vec{v} \]
impl BitOr<GaVector> for GaBivector {
    type Output = GaVector;

    fn bitor(self: GaBivector, b: GaVector) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec | b,
        }
    }
}

// vector/bivector dot
// \[ \vec{v} \wedge \vec{B} \]
impl BitOr<GaBivector> for GaVector {
    type Output = GaVector;

    fn bitor(self: GaVector, b: GaBivector) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self ^ b.mvec,
        }
    }
}

#[cfg(test)]
mod bivector_dot {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn bivector_bivector_dot() {
        let bivector1 = GaBivector::new(3.0, 5.0, 4.0);
        let bivector2 = GaBivector::new(2.0, 1.0, 6.0);
        let scalarMul = bivector1 * bivector2;
        let scalarDot = bivector1 | bivector2;
        assert_relative_eq!(scalarMul, scalarDot, max_relative = 0.000001);
    }

    #[test]
    fn bivector_vector_dot() {
        let bivector = GaBivector::new(3.0, 5.0, 4.0);
        let vector = GaVector::new(2.0, 1.0, 6.0);
        let vectorRes = bivector | vector;
        assert_relative_eq!(vectorRes[1], 33.0, max_relative = 0.000001);
        assert_relative_eq!(vectorRes[2], 18.0, max_relative = 0.000001);
        assert_relative_eq!(vectorRes[3], -14.0, max_relative = 0.000001);
    }

    // #[test]
    // fn vector_bivector_wedge() {}
}

// Add
// Vector addition
impl Add for GaBivector {
    type Output = GaBivector;

    fn add(self: GaBivector, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec + b.mvec,
        }
    }
}

// Sub
// Vector subtraction
impl Sub for GaBivector {
    type Output = GaBivector;

    fn sub(self: GaBivector, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec - b.mvec,
        }
    }
}

// scalar/bivector addition
impl Add<GaBivector> for f32 {
    type Output = GaBivector;

    fn add(self: f32, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self + b.mvec,
        }
    }
}

// bivector/scalar addition
impl Add<f32> for GaBivector {
    type Output = GaBivector;

    fn add(self: GaBivector, b: f32) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec + b,
        }
    }
}

// multivector/bivector addition
impl Add<GaBivector> for GaMultivector {
    type Output = GaMultivector;

    fn add(self: GaMultivector, b: GaBivector) -> GaMultivector {
        self + b.mvec
    }
}

// bivector/multivector addition
impl Add<GaMultivector> for GaBivector {
    type Output = GaMultivector;

    fn add(self: GaBivector, b: GaMultivector) -> GaMultivector {
        self.mvec + b
    }
}

// bivector/vector multiplication
impl Add<GaVector> for GaBivector {
    type Output = GaMultivector;

    fn add(self: GaBivector, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector addition
impl Add<GaBivector> for GaVector {
    type Output = GaMultivector;

    fn add(self: GaVector, b: GaBivector) -> GaMultivector {
        self * b.mvec
    }
}

// scalar/bivector addition
impl Sub<GaBivector> for f32 {
    type Output = GaBivector;

    fn sub(self: f32, b: GaBivector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self - b.mvec,
        }
    }
}

// bivector/scalar subtraction
impl Sub<f32> for GaBivector {
    type Output = GaBivector;

    fn sub(self: GaBivector, b: f32) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec - b,
        }
    }
}

// multivector/vector subtraction
impl Sub<GaBivector> for GaMultivector {
    type Output = GaMultivector;

    fn sub(self: GaMultivector, b: GaBivector) -> GaMultivector {
        self - b.mvec
    }
}

// vector/multivector subtraction
impl Sub<GaMultivector> for GaBivector {
    type Output = GaMultivector;

    fn sub(self: GaBivector, b: GaMultivector) -> GaMultivector {
        self.mvec - b
    }
}

// bivector/vector subtraction
impl Sub<GaVector> for GaBivector {
    type Output = GaMultivector;

    fn sub(self: GaBivector, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector subtration
impl Sub<GaBivector> for GaVector {
    type Output = GaMultivector;

    fn sub(self: GaVector, b: GaBivector) -> GaMultivector {
        self * b.mvec
    }
}

// the norm of a multivector |A|
// \[|A|^2=\left< A\^dag A \right>_0\]
impl GaBivector {
    pub fn Norm(self: Self) -> f32 {
        self.mvec.Norm()
    }
}

// the inverse of a multivector A^-1
// \[A^{-1}=\frac{A^\dag}{|A|^2}\]
impl GaBivector {
    pub fn Inverse(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Inverse(),
        }
    }
}
