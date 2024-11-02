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
use super::geometric_algebra_vector::GaVector;

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub struct GaRotor {
    mvec: GaMultivector,
}

impl GaRotor {
    pub fn zero() -> Self {
        Self {
            mvec: GaMultivector::zero(),
        }
    }

    pub fn new(angle_radians: f32, bivector: GaBivector) -> Self {
        let mut res = Self::zero();

        let bivector_norm = bivector.Norm();

        res[0] = cosf(angle_radians / 2.0);
        res[4] = sinf(angle_radians / 2.0) * bivector[4] / bivector_norm;
        res[5] = sinf(angle_radians / 2.0) * bivector[5] / bivector_norm;
        res[6] = sinf(angle_radians / 2.0) * bivector[6] / bivector_norm;
        res
    }

    pub fn new_vectors(angle_radians: f32, vector1: GaVector, vector2: GaVector) -> Self {
        let mut res = Self::zero();
        let bivector = vector1 ^ vector2;
        let bivector_norm = bivector.Norm();

        res[0] = cosf(angle_radians / 2.0);
        res[4] = sinf(angle_radians / 2.0) * bivector[4] / bivector_norm;
        res[5] = sinf(angle_radians / 2.0) * bivector[5] / bivector_norm;
        res[6] = sinf(angle_radians / 2.0) * bivector[6] / bivector_norm;
        res
    }
}

impl Index<usize> for GaRotor {
    type Output = f32;

    fn index<'a>(&'a self, index: usize) -> &'a Self::Output {
        &self.mvec[index]
    }
}

impl IndexMut<usize> for GaRotor {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Self::Output {
        &mut self.mvec[index]
    }
}

// Reverse
// \[ R^\dag \]
// Reverse the order of the basis blades.
impl GaRotor {
    pub fn Reverse(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec.Reverse(),
        }
    }
}

#[cfg(test)]
mod rotor_reverse {
    use super::*;
    use approx::assert_relative_eq;
    // The reverse of the geometric product of to rotors is the geometric product of the reverse rotors flipped
    // \[ (R_1R_2)^\dag = R_2^\dag R_1^\dag\]
    #[test]
    fn rotor_rotor_reverse() {
        let angle1 = PI / 4.0;
        let rotation_plane = GaBivector::new(3.0, 2.0, 10.0);
        let rotor1 = GaRotor::new(angle1, rotation_plane);
        let angle2 = PI / 2.0;
        let vector1 = GaVector::new(2.0, -3.0, -1.0);
        let vector2 = GaVector::new(-5.0, 3.0, 4.0);
        let rotor2 = GaRotor::new_vectors(angle2, vector1, vector2);

        assert_relative_eq!(
            (rotor1 * rotor2).Reverse()[0],
            (rotor2.Reverse() * rotor1.Reverse())[0],
            max_relative = 0.000001
        );
        assert_relative_eq!(
            (rotor1 * rotor2).Reverse()[4],
            (rotor2.Reverse() * rotor1.Reverse())[4],
            max_relative = 0.000001
        );
        assert_relative_eq!(
            (rotor1 * rotor2).Reverse()[5],
            (rotor2.Reverse() * rotor1.Reverse())[5],
            max_relative = 0.000001
        );
        assert_relative_eq!(
            (rotor1 * rotor2).Reverse()[6],
            (rotor2.Reverse() * rotor1.Reverse())[6],
            max_relative = 0.000001
        );
    }
}

// Dual
// Poincare duality operator.
impl GaRotor {
    pub fn Dual(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec.Dual(),
        }
    }
}

impl Not for GaRotor {
    type Output = GaRotor;

    fn not(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor { mvec: !self.mvec }
    }
}

// Conjugate
// Clifford Conjugation
impl GaRotor {
    pub fn Conjugate(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec.Conjugate(),
        }
    }
}

// Involute
// Main involution
impl GaRotor {
    pub fn Involute(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec.Involute(),
        }
    }
}

// Mul
// The geometric product.
impl Mul for GaRotor {
    type Output = GaRotor;

    fn mul(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec * b.mvec,
        }
    }
}

// Wedge
// The outer product. (MEET)
impl BitXor for GaRotor {
    type Output = GaRotor;

    fn bitxor(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec ^ b.mvec,
        }
    }
}

// Vee
// The regressive product. (JOIN)
impl BitAnd for GaRotor {
    type Output = GaRotor;

    fn bitand(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec & b.mvec,
        }
    }
}

// Dot
// The inner product.
impl BitOr for GaRotor {
    type Output = GaRotor;

    fn bitor(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec | b.mvec,
        }
    }
}

// Add
// Vector addition
impl Add for GaRotor {
    type Output = GaRotor;

    fn add(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec + b.mvec,
        }
    }
}

// Sub
// Vector subtraction
impl Sub for GaRotor {
    type Output = GaRotor;

    fn sub(self: GaRotor, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec - b.mvec,
        }
    }
}

// scalar/bivector multiplication
impl Mul<GaRotor> for f32 {
    type Output = GaRotor;

    fn mul(self: f32, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self * b.mvec,
        }
    }
}

// bivector/scalar multiplication
impl Mul<f32> for GaRotor {
    type Output = GaRotor;

    fn mul(self: GaRotor, b: f32) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec * b,
        }
    }
}

// bivector/multivector multiplication
impl Mul<GaMultivector> for GaRotor {
    type Output = GaMultivector;

    fn mul(self: GaRotor, b: GaMultivector) -> GaMultivector {
        self.mvec * b
    }
}

// multivector/bivector multiplication
impl Mul<GaRotor> for GaMultivector {
    type Output = GaMultivector;

    fn mul(self: GaMultivector, b: GaRotor) -> GaMultivector {
        self * b.mvec
    }
}

// bivector/vector multiplication
impl Mul<GaVector> for GaRotor {
    type Output = GaMultivector;

    fn mul(self: GaRotor, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector multiplication
impl Mul<GaRotor> for GaVector {
    type Output = GaMultivector;

    fn mul(self: GaVector, b: GaRotor) -> GaMultivector {
        self * b.mvec
    }
}

// scalar/bivector addition
impl Add<GaRotor> for f32 {
    type Output = GaRotor;

    fn add(self: f32, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self + b.mvec,
        }
    }
}

// bivector/scalar addition
impl Add<f32> for GaRotor {
    type Output = GaRotor;

    fn add(self: GaRotor, b: f32) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec + b,
        }
    }
}

// multivector/bivector addition
impl Add<GaRotor> for GaMultivector {
    type Output = GaMultivector;

    fn add(self: GaMultivector, b: GaRotor) -> GaMultivector {
        self + b.mvec
    }
}

// bivector/multivector addition
impl Add<GaMultivector> for GaRotor {
    type Output = GaMultivector;

    fn add(self: GaRotor, b: GaMultivector) -> GaMultivector {
        self.mvec + b
    }
}

// bivector/vector multiplication
impl Add<GaVector> for GaRotor {
    type Output = GaMultivector;

    fn add(self: GaRotor, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector addition
impl Add<GaRotor> for GaVector {
    type Output = GaMultivector;

    fn add(self: GaVector, b: GaRotor) -> GaMultivector {
        self * b.mvec
    }
}

// scalar/bivector addition
impl Sub<GaRotor> for f32 {
    type Output = GaRotor;

    fn sub(self: f32, b: GaRotor) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self - b.mvec,
        }
    }
}

// bivector/scalar subtraction
impl Sub<f32> for GaRotor {
    type Output = GaRotor;

    fn sub(self: GaRotor, b: f32) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec - b,
        }
    }
}

// multivector/vector subtraction
impl Sub<GaRotor> for GaMultivector {
    type Output = GaMultivector;

    fn sub(self: GaMultivector, b: GaRotor) -> GaMultivector {
        self - b.mvec
    }
}

// vector/multivector subtraction
impl Sub<GaMultivector> for GaRotor {
    type Output = GaMultivector;

    fn sub(self: GaRotor, b: GaMultivector) -> GaMultivector {
        self.mvec - b
    }
}

// bivector/vector subtraction
impl Sub<GaVector> for GaRotor {
    type Output = GaMultivector;

    fn sub(self: GaRotor, b: GaVector) -> GaMultivector {
        self.mvec * b
    }
}

// vector/bivector subtration
impl Sub<GaRotor> for GaVector {
    type Output = GaMultivector;

    fn sub(self: GaVector, b: GaRotor) -> GaMultivector {
        self * b.mvec
    }
}

// the norm of a multivector |A|
// \[|A|^2=\left< A\^dag A \right>_0\]
impl GaRotor {
    pub fn Norm(self: Self) -> f32 {
        self.mvec.Norm()
    }
}

// the inverse of a multivector A^-1
// \[A^{-1}=\frac{A^\dag}{|A|^2}\]
impl GaRotor {
    pub fn Inverse(self: Self) -> GaRotor {
        super::geometric_algebra_rotor::GaRotor {
            mvec: self.mvec.Inverse(),
        }
    }
}
