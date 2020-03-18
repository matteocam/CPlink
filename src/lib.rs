#[macro_use]

extern crate cute;

mod matrix;
mod snark;

extern crate bn;
extern crate rand;

extern crate rug;
extern crate algebra;

extern crate rand_xorshift;

use bn::*;

use crate::matrix::*;
use crate::snark::*;

use cpsnarks_set::commitments::{
    integer::IntegerCommitment, pedersen::PedersenCommitment, Commitment,
};

pub fn test() {
    
    use rug::Integer;
    use algebra::jubjub::JubJubProjective;
    use rand_xorshift::XorShiftRng;
    use rand::SeedableRng;

    let rng = rand::thread_rng();

    let mut pp = PP {l:1, t: 2, rng:rng};

    let m = Matrix::new(pp.l, pp.t, &vec![G1::one(), G1::one()]);

    let x:Vec<Fr> = vec![Fr::one(), Fr::zero()];

    let y:VecG = vec![G1::one()];

    let (ek, vk) = keygen(&mut pp, m);

    let pi = prove(&mut pp, &ek, &x);

    let b = verify(&pp, &vk, &y, &pi);

    println!("Result is {}.", b);

}
