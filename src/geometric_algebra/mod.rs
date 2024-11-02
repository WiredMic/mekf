#![allow(dead_code)]
#![allow(unused_imports)]
// This is an implemention of 3D vectorspace geometric algebra
pub mod geometric_algebra_bivector;
pub mod geometric_algebra_multivector;
pub mod geometric_algebra_rotor;
pub mod geometric_algebra_vector;

pub mod rotations {
    use super::geometric_algebra_multivector::GaMultivector;
    use super::geometric_algebra_rotor::GaRotor;
    use super::geometric_algebra_vector::GaVector;

    // \[ R^\dag \vec{v} R \]
    // A rotation in g3 is a sandwitch product of a rotor ( R ) and a vector ( v )
    pub fn ga_rotation(rotor: GaRotor, vector: GaVector) -> GaVector {
        // normilise roter
        let rotor = rotor * (1.0 / rotor.Norm());
        let norm = vector.Norm();

        super::geometric_algebra_vector::GaVector {
            mvec: (rotor.Reverse() * vector * rotor) * (1.0 / norm) * norm,
        }
    }

    #[cfg(test)]
    mod ga_rotation_tests {

        use crate::geometric_algebra::geometric_algebra_bivector::GaBivector;

        use super::*;
        use approx::assert_relative_eq;
        use core::f32::consts::PI;

        #[test]
        fn rotor() {
            let angle: f32 = PI / 2.0;
            let bivector = GaBivector::new(1.0, 0.0, 0.0);
            let rotor = GaRotor::new(angle, bivector);
            assert_relative_eq!(rotor[0], 0.7071067, max_relative = 0.000001);
            assert_relative_eq!(rotor[4], 0.7071067, max_relative = 0.000001);
        }

        #[test]
        fn vec_rot_quarter_e1e2() {
            let vector = GaVector::new(3.0, 0.0, 0.0);
            let angle: f32 = PI / 2.0;
            let bivector = GaBivector::new(1.0, 0.0, 0.0);
            let rotor = GaRotor::new(angle, bivector);
            let vector_rot = ga_rotation(rotor, vector);

            assert_relative_eq!(vector_rot[1], 0.0, max_relative = 0.000001);
            assert_relative_eq!(vector_rot[2], 3.0, max_relative = 0.000001);
            assert_relative_eq!(vector_rot[3], 0.0, max_relative = 0.000001);
        }
    }
}
