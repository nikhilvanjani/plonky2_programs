extern crate anyhow;
extern crate plonky2_programs;

use anyhow::Result;
use plonky2_programs::test_is_zero; // Import from your library

fn main() -> Result<()> {
    test_is_zero()
}