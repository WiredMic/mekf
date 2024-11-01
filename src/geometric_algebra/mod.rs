#![allow(dead_code)]

// This is an implemention of 3D vectorspace geometric algebra
mod geometric_algrbra_multivector;

pub mod rotations {
    use super::geometric_algrbra_multivector::GaMultivector;
    // \[ R^\dag \vec{v} R \]
    // A rotation in g3 is a sandwitch product of a rotor ( R ) and a vector ( v )
    pub fn ga_rotation(rotor: GaMultivector, vector: GaMultivector) -> GaMultivector {
        // normilise roter
        let rotor = rotor * (1.0 / rotor.Norm());

        rotor.Reverse() * vector * rotor
    }

    #[cfg(test)]
    mod ga_rotation_tests {

        use super::GaMultivector;
        use crate::geometric_algebra::rotations::ga_rotation;
        use approx::assert_relative_eq;
        use core::f32::consts::PI;

        #[test]
        fn rotor() {
            let angle: f32 = PI / 2.0;
            let rotor = GaMultivector::new_rotor(angle, 1.0, 0.0, 0.0);

            assert_relative_eq!(rotor[0], 0.7071067, max_relative = 0.000001);
            assert_relative_eq!(rotor[4], 0.7071067, max_relative = 0.000001);
        }

        #[test]
        fn vec_rot_quarter_e1e2() {
            let vector = GaMultivector::new_vector(3.0, 0.0, 0.0);
            let angle: f32 = PI / 2.0;
            let rotor = GaMultivector::new_rotor(angle, 1.0, 0.0, 0.0);
            let vector_rot = ga_rotation(rotor, vector);

            assert_relative_eq!(vector_rot[1], 0.0, max_relative = 0.000001);
            assert_relative_eq!(vector_rot[2], 3.0, max_relative = 0.000001);
            assert_relative_eq!(vector_rot[3], 0.0, max_relative = 0.000001);
        }
    }
}
