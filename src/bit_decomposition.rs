extern crate plonky2;
extern crate anyhow;

use anyhow::Result;
use plonky2::field::types::{Field, Sample};
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};

// /// A generator used by the prover to calculate the `bit_decomposition` of a 
// /// given `value`, outside of the circuit, in order to supply it as an additional public input.
// #[derive(Debug, Default)]
// struct BitDecompositionGenerator<F: RichField + Extendable<D>, const D: usize> {
//     value: Target,
//     bit_decomposition: Vec<Target>,
//     _phantom: PhantomData<F>,
// }

// impl<F: RichField + Extendable<D>, const D: usize> SimpleGenerator<F, D>
//     for BitDecompositionGenerator<F, D>
// {
//     fn id(&self) -> String {
//         "BitDecompositionGenerator".to_string()
//     }

//     fn dependencies(&self) -> Vec<Target> {
//         vec![self.value]
//     }

//     fn run_once(
//         &self,
//         witness: &PartitionWitness<F>,
//         out_buffer: &mut GeneratedValues<F>,
//     ) -> Result<()> {
//         let value = witness.get_target(self.value);
//         let mut bit_decomposition = Vec::new();
//         let mut tmp_value = value.0;
//         for _ in 0..F::BITS {
//             let bit = ;
//             bit_decomposition.push(value & F::ONE == F::ONE);
//             tmp_value >>= 1;
//         }
//         println!("Bit Decomposition: {bit_decomposition}");

//         out_buffer.set_target(self.bit_decomposition, bit_decomposition)
//     }

//     fn serialize(&self, dst: &mut Vec<u8>, _common_data: &CommonCircuitData<F, D>) -> IoResult<()> {
//         dst.write_target(self.value)?;
//         dst.write_target(self.bit_decomposition)
//     }
//     fn deserialize(src: &mut Buffer, _common_data: &CommonCircuitData<F, D>) -> IoResult<Self> {
//         let value = src.read_target()?;
//         let bit_decomposition = src.read_target()?;
//         Ok(Self {
//             value,
//             bit_decomposition,
//             _phantom: PhantomData,
//         })
//     }
// }


/// An example of using Plonky2 to prove a statement of the form
/// "I know bit decomposition of value".
pub fn test_bit_decomposition() -> Result<()> {
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    // The secret value.
    let value = builder.add_virtual_target();
    let bit_decomposition = builder.split_le(value, F::BITS);

    // Registered as a public input.
    builder.register_public_input(value);
    for bit in bit_decomposition.iter() {
        builder.register_public_input(bit.target);
    }

    let mut pw = PartialWitness::new();
    pw.set_target(value, F::from_canonical_usize(42))?;
    // pw.set_target(value, F::rand())?;


    let data = builder.build::<C>();
    let proof = data.prove(pw)?;

    println!(
        "Bit Decomposition of Value {} is {:?}",
        proof.public_inputs[0], &proof.public_inputs[1..],
    );

    data.verify(proof)
}
