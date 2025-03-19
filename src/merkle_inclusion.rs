extern crate plonky2;
extern crate anyhow;

use anyhow::Result;
use plonky2::field::types::{Field, Sample};
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
use plonky2::hash::poseidon::Poseidon;
use plonky2::iop::generator::SimpleGenerator;
use plonky2::hash::hash_types::RichField;
use plonky2::field::extension::Extendable;
use plonky2::iop::witness::PartitionWitness;
use plonky2::iop::generator::GeneratedValues;
use plonky2::plonk::circuit_data::CommonCircuitData;
use plonky2::util::serialization::IoResult;
use plonky2::util::serialization::Buffer;
use plonky2::iop::target::Target;
use plonky2::iop::target::BoolTarget;
use std::convert::TryInto;
use plonky2::plonk::config::Hasher;
use plonky2::field::goldilocks_field::GoldilocksField;
use plonky2::iop::witness::Witness;
use plonky2::util::serialization::Write;
use plonky2::util::serialization::Read;
use plonky2::field::types::PrimeField64;
use plonky2::hash::hash_types::HashOut;
use plonky2::hash::hash_types::HashOutTarget;

/// F is the native field we use everywhere.  Currently it's Goldilocks from plonky2
pub type F = GoldilocksField;
/// C is the Plonky2 config used to work with Plonky2 recursion.
pub type C = PoseidonGoldilocksConfig;
/// D defines the extension degree of the field used in the Plonky2 proofs (quadratic extension).
pub const D: usize = 2;

/// An example of using Plonky2 to prove a statement of the form
/// "I know x such that y = MerkleHash(..., x, ...)".
/// This is a simple example of a Merkle inclusion proof.
pub fn test_merkle_inclusion<H: Hasher<F, Hash = HashOut<GoldilocksField>> + std::marker::Sync + std::marker::Send + 'static>(num_leaf: usize) -> Result<()> {

    // Compute the depth d of merkle tree from the number of leaves num_leaf.
    let d = 1 + (num_leaf as f64).log2().ceil() as usize;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    // The secret value.
    let depth = builder.add_virtual_target();
    println!("Depth: {:?}", depth);
    let x = builder.add_virtual_target();
    println!("x: {:?}", x);
    let x_index = builder.add_virtual_target();
    println!("x_index: {:?}", x_index);
    let root_hash = builder.add_virtual_hash();
    println!("root_hash: {:?}", root_hash);
    let sibling_hashes = (0..d-1).map(|_| builder.add_virtual_hash()).collect::<Vec<_>>();
    println!("sibling_hashes: {:?}", sibling_hashes);
    let local_hashes = (0..d).map(|_| builder.add_virtual_hash()).collect::<Vec<_>>();
    println!("local_hashes: {:?}", local_hashes);

    // Constraints.
    let x_index_bits = builder.split_le(x_index, (d-1).try_into().unwrap());
    println!("x_index_bits: {:?}", x_index_bits);
    builder.connect_array(root_hash.elements, local_hashes[0].elements);
    // The arithmetic circuit.
    builder.add_simple_generator(MerkleInclusionGenerator::<H> {
        depth, 
        x, 
        x_index_bits, 
        sibling_hashes: sibling_hashes.clone(), 
        local_hashes: local_hashes.clone(),
        _phantom: std::marker::PhantomData,
    });

    // Registered as a public input.
    builder.register_public_input(x);
    builder.register_public_inputs(&root_hash.elements);
    
    // Generate the inputs.
    let (x_val, x_index_val, sibling_hashes_val, root_hash_val) = generate_merkle_inclusion_inputs::<F, H>(num_leaf);

    let mut pw = PartialWitness::new();
    pw.set_target(depth, F::from_canonical_usize(d.try_into().unwrap()));
    pw.set_target(x, F::from_canonical_usize(x_val.try_into().unwrap()));
    pw.set_target(x_index, F::from_canonical_usize(x_index_val.try_into().unwrap()));
    pw.set_hash_target(root_hash, root_hash_val);
    for i in 0..(d-1).try_into().unwrap() {
        pw.set_hash_target(sibling_hashes[i], sibling_hashes_val[i])?;
    }

    let data = builder.build::<C>();
    println!("data SET");
    println!("pw: {:?}", pw);
    let proof = data.prove(pw)?;
    println!("proof SET");

    println!(
        "Merkle inclusion proof for x = {}, root = {}",
        proof.public_inputs[0], proof.public_inputs[1],
    );

    data.verify(proof)
}

fn generate_merkle_inclusion_inputs<F: RichField, H: Hasher<F>>(num_leaf: usize) -> (usize, usize, Vec<H::Hash>, H::Hash) {
    // Depth is 1 more than usual since leaf nodes are hashed first at depth 0 individually.
    // Eg: If num_leaf = 8, then depth = 4.
    let depth = 1 + (num_leaf as f64).log2().ceil() as i64;

    // leaf[i] = i for all i in [0, num_leaf)
    let leaf = (0..num_leaf).collect::<Vec<_>>();
    let index = 2;
    let index_bits = (0..depth-1).rev().map(|i| (index >> i) & 1).collect::<Vec<_>>();
    println!("Index bits: {:?}", index_bits);

    let mut own_indices = vec![0; (depth - 1) as usize];
    for j in 0..(depth - 1) as usize {
        if j == 0 {
            own_indices[0] = 1 + index_bits[0];
        } else {
            own_indices[j] = 2 * own_indices[j - 1] + 1 + index_bits[j];
        }
    }

    let mut sibling_indices = vec![0; (depth - 1) as usize];
    for j in 0..(depth - 1) as usize {
        sibling_indices[j] = if index_bits[j] == 0 { own_indices[j] + 1 } else { own_indices[j] - 1 };
    }
    println!("Own indices: {:?}", own_indices);
    println!("Sibling indices: {:?}", sibling_indices);

    // let mut all_hashes = vec![H::Hash; 2 * num_leaf - 1];
    // let mut sibling_hashes = vec![H::Hash; (depth - 1) as usize];
    // let mut all_hashes = vec![H::hash_no_pad(&[F::ZERO]); 2 * num_leaf - 1];
    let mut all_hashes = vec![H::hash_no_pad(&[F::ZERO]); 2_usize.pow(depth.try_into().unwrap()) - 1];
    let mut sibling_hashes = vec![H::hash_no_pad(&[F::ZERO]); (depth - 1) as usize];
    // let mut all_hashes = Vec::with_capacity(2 * num_leaf - 1);
    // let mut sibling_hashes = Vec::with_capacity((depth - 1) as usize);
    println!("depth: {:?}", depth);
    println!("leaf: {:?}", leaf);
    for j in (0..depth).rev() {
        for i in (2_usize.pow(j as u32) - 1)..(2_usize.pow(j as u32 + 1) - 1) {
            if j == depth - 1 {
                println!("i: {}", i);
                if i - (2_usize.pow(j as u32) - 1) < num_leaf {
                    all_hashes[i] = H::hash_no_pad(&[F::from_canonical_usize(leaf[i - (2_usize.pow(j as u32) - 1)])]);
                } else {
                    all_hashes[i] = H::hash_no_pad(&[F::ZERO]);
                }
            } else {
                all_hashes[i] = H::two_to_one(all_hashes[2 * i + 1], all_hashes[2 * i + 2]);
            }
        }
        if j != depth - 1 {
            sibling_hashes[j as usize] = all_hashes[sibling_indices[j as usize]];
        }
    }
    for i in 0..(2_usize.pow(depth as u32) - 1) {
        println!("all_hashes[{}]: {:?}", i, all_hashes[i]);
    }
    for i in 0..(depth - 1) as usize {
        println!("sibling_hashes[{}]: {:?}", i, sibling_hashes[i]);
    }

    (leaf[index], index,  sibling_hashes, all_hashes[0])
}

/// A generator for Merkle inclusion proofs.
// #[derive(Default)]
#[derive(Debug, Default)]
struct MerkleInclusionGenerator<H: Hasher<F, Hash = HashOut<GoldilocksField>>> {
    depth: Target,
    x: Target,
    x_index_bits: Vec<BoolTarget>,
    sibling_hashes: Vec<HashOutTarget>,
    local_hashes: Vec<HashOutTarget>,
    _phantom: std::marker::PhantomData<H>,
}

/// Implementation of `SimpleGenerator` for `MerkleInclusionGenerator`.
/// This generator is responsible for generating the witness values for the Merkle inclusion proof.
// impl<H: Hasher<GoldilocksField, Hash = GoldilocksField> + std::marker::Sync + std::marker::Send, const D: usize> SimpleGenerator<GoldilocksField, D> for MerkleInclusionGenerator<H> {
impl<H: Hasher<GoldilocksField, Hash = HashOut<GoldilocksField>> + std::marker::Sync + std::marker::Send + 'static> SimpleGenerator<GoldilocksField, 2> for MerkleInclusionGenerator<H> {
    fn id(&self) -> String {
        "MerkleInclusionGenerator".to_string()
    }

    fn dependencies(&self) -> Vec<Target> {
        let mut deps = Vec::new();
        deps.push(self.depth);
        deps.push(self.x);
        for i in 0..self.x_index_bits.len() {
            deps.push(self.x_index_bits[i].target);
        }
        for i in 0..self.sibling_hashes.len() {
            for j in 0..self.sibling_hashes[i].elements.len() {
                deps.push(self.sibling_hashes[i].elements[j]);
            }
            // deps.push(self.sibling_hashes[i]);
        }
        // for i in 0..self.local_hashes.len() {
        //     for j in 0..self.local_hashes[i].elements.len() {
        //         deps.push(self.local_hashes[i].elements[j]);
        //     }
        //     // deps.push(self.local_hashes[i]);
        // }
        println!("MerkleInclusionGenerator dependencies: {:?}", deps);
        deps
    }

    fn run_once(
        &self,
        witness: &PartitionWitness<F>,
        out_buffer: &mut GeneratedValues<F>,
    ) -> Result<()> {
        println!("Running MerkleInclusionGenerator::run_once");
        let depth = witness.get_target(self.depth);
        let d : usize = depth.to_canonical_u64() as usize;

        let x = witness.get_target(self.x);
        // self.x_index_bits is in little endian notation, eg, 3 is 110. But for computation here, we want it in big endian notation, so we reverse it using rev().
        let x_index_bits = self.x_index_bits.iter().rev().map(|&bit| witness.get_target(bit.target).to_canonical_u64()).collect::<Vec<_>>();
        println!("x_index_bits: {:?}", x_index_bits);
        let sibling_hashes = self.sibling_hashes.iter().map(|&hash| witness.get_hash_target(hash)).collect::<Vec<_>>();
        let mut local_hashes = vec![H::Hash::default(); d];

        for j in (0..d).rev() {
            if (j == d - 1) {
                local_hashes[j] = H::hash_no_pad(&[x]);
            }
            else {
                let left = if x_index_bits[j] == 0 { local_hashes[j+1] } else { sibling_hashes[j] };
                let right = if x_index_bits[j] == 0 { sibling_hashes[j] } else { local_hashes[j+1] };
                println!("j: {} =======", j);
                println!("left: {:?}", left);
                println!("right: {:?}", right);
                local_hashes[j] = H::two_to_one(left, right);
            }
        }
        for i in 0..d {
            println!("local_hashes[{}]: {:?}", i, local_hashes[i]);
            out_buffer.set_hash_target(self.local_hashes[i], local_hashes[i]);
        }
        println!("out_buffer: {:?}", out_buffer);
        Ok(())
    }

    fn serialize(&self, dst: &mut Vec<u8>, _common_data: &CommonCircuitData<F, 2>) -> IoResult<()> {
        dst.write_target(self.depth)?;
        dst.write_target(self.x)?;
        dst.write_target_vec(&self.x_index_bits.iter().map(|b| b.target).collect::<Vec<Target>>())?;
        for i in 0..self.sibling_hashes.len() {
            dst.write_target_hash(&self.sibling_hashes[i])?;
            // for j in 0..self.sibling_hashes[i].elements.len() {
            //     dst.write_target(self.sibling_hashes[i].elements[j])?;
            // }
        }
        // dst.write_target_vec(&self.sibling_hashes)?;
        for i in 0..self.local_hashes.len() {
            dst.write_target_hash(&self.local_hashes[i])?;
            // for j in 0..self.local_hashes[i].elements.len() {
            //     dst.write_target(self.local_hashes[i].elements[j])?;
            // }
        }
        // dst.write_target_vec(&self.local_hashes)
        Ok(())
    }

    fn deserialize(src: &mut Buffer, _common_data: &CommonCircuitData<F, 2>) -> IoResult<Self> {
        let depth = src.read_target()?;
        let x = src.read_target()?;
        let x_index_bits_tmp = src.read_target_vec()?;
        let x_index_bits = x_index_bits_tmp.iter().map(|b| BoolTarget::new_unsafe(*b)).collect::<Vec<BoolTarget>>();
        let d = x_index_bits.len() + 1; 
        // let d = depth.to_canonical_u64() as usize;
        let mut sibling_hashes = Vec::with_capacity(d - 1);
        let mut local_hashes = Vec::with_capacity(d);
        for i in 0..d-1 {
            sibling_hashes[i] = src.read_target_hash()?;
        }
        for i in 0..d {
            local_hashes[i] = src.read_target_hash()?;
        }
        Ok(Self {
            depth,
            x,
            x_index_bits,
            sibling_hashes,
            local_hashes,
            _phantom: std::marker::PhantomData,
        })
    }
}