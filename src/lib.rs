extern crate anyhow;
extern crate plonky2;

use anyhow::Result;
pub mod bit_decomposition;
pub mod is_zero;
pub mod merkle_inclusion;

pub use bit_decomposition::test_bit_decomposition;
pub use is_zero::test_is_zero;
pub use merkle_inclusion::test_merkle_inclusion;

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use bit_dec::test_bit_decomposition;
//     use anyhow::Result;

//     #[test]
//     fn test_bit_decomposition_works() -> Result<()> {
//         test_bit_decomposition()
//     }
// }