extern crate plonky2;
extern crate anyhow;
extern crate plonky2_programs;

use anyhow::Result;
use plonky2_programs::test_merkle_inclusion; // Import from your library
use plonky2::hash::poseidon::Poseidon;
use plonky2::field::goldilocks_field::GoldilocksField;

fn main() -> Result<()> {
    test_merkle_inclusion::<plonky2::hash::poseidon::PoseidonHash>(5)
    // test_merkle_inclusion::<GoldilocksField>(2)
}