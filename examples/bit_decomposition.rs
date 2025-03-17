extern crate anyhow;
extern crate plonky2_programs;

use anyhow::Result;
use plonky2_programs::test_bit_decomposition; // Import from your library

fn main() -> Result<()> {
    test_bit_decomposition()
}