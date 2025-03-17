extern crate plonky2;
extern crate anyhow;

use anyhow::Result;
use plonky2::field::types::{Field, Sample};
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};

/// An example of using Plonky2 to prove a statement of the form
/// "I know x such that y = isZero(x), where y = 1 if x = 0, else y = 0".
pub fn test_is_zero() -> Result<()> {
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    // The secret value.
    let x = builder.add_virtual_target();

    // The public values.
    let zero = builder.constant(F::ZERO);
    let y = builder.is_equal(x, zero);

    // Registered as a public input.
    builder.register_public_input(x);
    builder.register_public_input(y.target);

    let mut pw = PartialWitness::new();
    // pw.set_target(x, F::from_canonical_usize(42))?;
    pw.set_target(x, F::from_canonical_usize(0))?;

    let data = builder.build::<C>();
    let proof = data.prove(pw)?;

    println!(
        "x =  {}, y = {}",
        proof.public_inputs[0], &proof.public_inputs[1],
    );

    data.verify(proof)
}