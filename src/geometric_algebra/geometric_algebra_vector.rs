#![allow(unused_imports)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]
// #![feature(const_slice_len)]

use core::ops::{Add, BitAnd, BitOr, BitXor, Div, Index, IndexMut, Mul, Not, Sub};
use libm::{cosf, powf, sinf, sqrtf};

use super::geometric_algebra_bivector::GaBivector;
use super::geometric_algebra_multivector::basis_count;
use super::geometric_algebra_multivector::GaMultivector;

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct GaVector {
    pub mvec: GaMultivector,
}

impl GaVector {
    pub fn zero() -> Self {
        Self {
            mvec: GaMultivector::zero(),
        }
    }

    pub fn new(e1: f32, e2: f32, e3: f32) -> Self {
        let mut ret = Self::zero();
        ret.mvec[1] = e1;
        ret.mvec[2] = e2;
        ret.mvec[3] = e3;
        ret
    }
}

impl Index<usize> for GaVector {
    type Output = f32;

    fn index<'a>(&'a self, index: usize) -> &'a Self::Output {
        &self.mvec[index]
    }
}

impl IndexMut<usize> for GaVector {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Self::Output {
        &mut self.mvec[index]
    }
}

// Reverse
// \[ \vec{v}^\dag\]
// Reverse the order of the basis blades.
impl GaVector {
    pub fn Reverse(self: Self) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec.Reverse(),
        }
    }
}

// Dual
// \[ \vec{v}\mathrm{I} = \vec{B}\]
// Poincare duality operator.
impl GaVector {
    pub fn Dual(self: Self) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec.Dual(),
        }
    }
}

#[cfg(test)]
mod vector_dual {
    use super::*;
    #[test]
    fn vector_to_bivector() {
        let vector: GaVector = GaVector::new(1.0, 2.0, 3.0);
        let bivector: GaBivector = vector.Dual();
        assert_eq!(vector[1], bivector[6]);
        assert_eq!(vector[2], bivector[5]);
        assert_eq!(vector[3], bivector[4]);
    }
}

impl Not for GaVector {
    type Output = GaVector;

    fn not(self: Self) -> GaVector {
        super::geometric_algebra_vector::GaVector { mvec: !self.mvec }
    }
}

// Conjugate
// \[ \bar{\vec{v}}\]
// Clifford Conjugation
impl GaVector {
    pub fn Conjugate(self: Self) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec.Conjugate(),
        }
    }
}

// Involute
// Main involution
impl GaVector {
    pub fn Involute(self: Self) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec.Involute(),
        }
    }
}

// Mul
// \[ \vec{u}\vec{v}= \vec{u}\cdot\vec{v}+\vec{u}\wedge\vec{v}\]
// The geometric product.
impl Mul for GaVector {
    type Output = GaMultivector;

    fn mul(self: GaVector, b: GaVector) -> GaMultivector {
        self.mvec * b.mvec
    }
}

// scalar/vector multiplication
// \[ a\vec{v}\]
impl Mul<GaVector> for f32 {
    type Output = GaVector;

    fn mul(self: f32, b: GaVector) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self * b.mvec,
        }
    }
}

// vector/scalar multiplication
// \[ \vec{v}a\]
impl Mul<f32> for GaVector {
    type Output = GaVector;

    fn mul(self: GaVector, b: f32) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec * b,
        }
    }
}

// vector/multivector multiplication
// \[ \vec{v}M\]
impl Mul<GaMultivector> for GaVector {
    type Output = GaMultivector;

    fn mul(self: GaVector, b: GaMultivector) -> GaMultivector {
        self.mvec * b
    }
}

// multivector/vector multiplication
// \[ M\vec{v} \]
impl Mul<GaVector> for GaMultivector {
    type Output = GaMultivector;

    fn mul(self: GaMultivector, b: GaVector) -> GaMultivector {
        self * b.mvec
    }
}

// Wedge
// \[ \vec{u}\wedge\vec{v} = \vec{B} \]
// The outer product. (MEET)
impl BitXor for GaVector {
    type Output = GaBivector;

    fn bitxor(self: GaVector, b: GaVector) -> GaBivector {
        super::geometric_algebra_bivector::GaBivector {
            mvec: self.mvec ^ b.mvec,
        }
    }
}

// vector/multivector wedge
// \[ \vec{v}\wedge M \]
impl BitXor<GaMultivector> for GaVector {
    type Output = GaMultivector;

    fn bitxor(self: GaVector, b: GaMultivector) -> GaMultivector {
        self.mvec ^ b
    }
}

// multivector/vector wedge
// \[ M\wedge \vec{v} \]
impl BitXor<GaVector> for GaMultivector {
    type Output = GaMultivector;

    fn bitxor(self: GaMultivector, b: GaVector) -> GaMultivector {
        self ^ b.mvec
    }
}

#[cfg(test)]
mod vector_wedge {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn vector_vector_wedge() {
        // 3e1+5e2+4e3
        let vector1 = GaVector::new(3.0, 5.0, 4.0);
        // 2e1+1e2+6e3
        let vector2 = GaVector::new(2.0, 1.0, 6.0);
        let bivector = vector1 ^ vector2;
        // −7e12​-10e31​+26e23
        assert_relative_eq!(bivector[4], -7.0, max_relative = 0.000001);
        assert_relative_eq!(bivector[5], -10.0, max_relative = 0.000001);
        assert_relative_eq!(bivector[6], 26.0, max_relative = 0.000001);
    }
}

// Vee
// The regressive product. (JOIN)
impl BitAnd for GaVector {
    type Output = GaMultivector;

    fn bitand(self: GaVector, b: GaVector) -> GaMultivector {
        self.mvec & b.mvec
    }
}

// Dot
// \[ \vec{u} \cdot \vec{v} = a \]
// The inner product.
impl BitOr for GaVector {
    type Output = f32;

    fn bitor(self: GaVector, b: GaVector) -> f32 {
        (self.mvec | b.mvec)[0]
    }
}

// vector/multivector dot
// \[ \vec{v} \cdot M \]
impl BitOr<GaMultivector> for GaVector {
    type Output = GaMultivector;

    fn bitor(self: GaVector, b: GaMultivector) -> GaMultivector {
        self.mvec | b
    }
}

// multivector/vector dot
// \[ M \cdot \vec{v}\]
impl BitOr<GaVector> for GaMultivector {
    type Output = GaMultivector;

    fn bitor(self: GaMultivector, b: GaVector) -> GaMultivector {
        self | b.mvec
    }
}

#[cfg(test)]
mod vector_dot {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn vector_vector_dot() {
        // 3e1+5e2+4e3
        let vector1: GaVector = GaVector::new(3.0, 5.0, 4.0);
        // 2e1+1e2+6e3
        let vector2: GaVector = GaVector::new(2.0, 1.0, 6.0);
        let scalar = vector1 | vector2;
        assert_relative_eq!(scalar, 35.0, max_relative = 0.000001);
    }
}

// Add
// Vector addition
impl Add for GaVector {
    type Output = GaVector;

    fn add(self: GaVector, b: GaVector) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec + b.mvec,
        }
    }
}

// scalar/vector addition
impl Add<GaVector> for f32 {
    type Output = GaMultivector;

    fn add(self: f32, b: GaVector) -> GaMultivector {
        self + b.mvec
    }
}

// vector/scalar addition
impl Add<f32> for GaVector {
    type Output = GaMultivector;

    fn add(self: GaVector, b: f32) -> GaMultivector {
        self.mvec + b
    }
}

// multivector/vector addition
impl Add<GaVector> for GaMultivector {
    type Output = GaMultivector;

    fn add(self: GaMultivector, b: GaVector) -> GaMultivector {
        self + b.mvec
    }
}

// vector/multivector addition
impl Add<GaMultivector> for GaVector {
    type Output = GaMultivector;

    fn add(self: GaVector, b: GaMultivector) -> GaMultivector {
        self.mvec + b
    }
}

// Sub
// Vector subtraction
impl Sub for GaVector {
    type Output = GaVector;

    fn sub(self: GaVector, b: GaVector) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec - b.mvec,
        }
    }
}

// scalar/vector subtraction
impl Sub<GaVector> for f32 {
    type Output = GaMultivector;

    fn sub(self: f32, b: GaVector) -> GaMultivector {
        self - b.mvec
    }
}

// vector/scalar subtraction
impl Sub<f32> for GaVector {
    type Output = GaMultivector;

    fn sub(self: GaVector, b: f32) -> GaMultivector {
        self.mvec - b
    }
}

// multivector/vector subtraction
impl Sub<GaVector> for GaMultivector {
    type Output = GaMultivector;

    fn sub(self: GaMultivector, b: GaVector) -> GaMultivector {
        self - b.mvec
    }
}

// vector/multivector subtraction
impl Sub<GaMultivector> for GaVector {
    type Output = GaMultivector;

    fn sub(self: GaVector, b: GaMultivector) -> GaMultivector {
        self.mvec - b
    }
}

// the norm of a multivector |A|
// \[ |A|^2=\left< A^\dag A \right>_0\]
impl GaVector {
    pub fn Norm(self: Self) -> f32 {
        self.mvec.Norm()
    }
}

// the inverse of a multivector A^-1
// \[A^{-1}=\frac{A^\dag}{|A|^2}\]
impl GaVector {
    pub fn Inverse(self: Self) -> GaVector {
        super::geometric_algebra_vector::GaVector {
            mvec: self.mvec.Inverse(),
        }
    }
}
