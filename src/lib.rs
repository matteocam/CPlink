#![feature(associated_type_defaults)]
#![allow(non_snake_case)]
#[macro_use]

extern crate cute;

mod matrix;
mod subspace_snark;
mod linking_snark;

extern crate rand;

extern crate rug;
extern crate algebra;

extern crate rand_xorshift;

//use bn::*;

use crate::matrix::*;
use crate::subspace_snark::*;
use crate::linking_snark::*;

use rand::Rng;

use cpsnarks_set::commitments::{
    integer::IntegerCommitment, pedersen::PedersenCommitment, Commitment,
};



use algebra_core::test_rng;
use std::ops::Add;
use algebra::{
    Zero,
    One,
    ProjectiveCurve,
    UniformRand,
    bls12_381::{
        g1, g2, Bls12_381, Fq, Fq12, Fq2, Fr, G1Affine, G1Projective, G2Affine, G2Projective,
    }
};

impl SparseLinAlgebra<Bls12_381> for Bls12_381 { } 


pub fn test() {
    
    use rug::Integer;
    use algebra::jubjub::JubJubProjective;
    //use rand_xorshift::XorShiftRng;
    //use rand::SeedableRng;

    let mut rng = test_rng();
    let g1 = G1Projective::rand(&mut rng);
    let g2 = G2Projective::rand(&mut rng);
    
    let mut pp =  PP::<G1Projective, G2Projective> {l:1, t: 2, g1, g2};

    let mut m = SparseMatrix::new(1, 2);
    m.insert_row_slice(0, 0, &vec![g1, g1]);

    let x:Vec<Fr> = vec![Fr::one(), Fr::zero()];

    let x_bad:Vec<Fr> = vec![Fr::one(), Fr::one()];


    let y:Vec<G1Projective> = vec![g1];

    let (ek, vk) = Bls12_381::keygen(&mut rng, &pp, m);

    let pi = Bls12_381::prove(&mut pp, &ek, &x);
    let pi_bad = Bls12_381::prove(&mut pp, &ek, &x_bad);

    let b = Bls12_381::verify(&pp, &vk, &y, &pi);
    let b_bad = Bls12_381::verify(&pp, &vk, &y, &pi_bad);


    println!("Result is {}.", b);
    println!("Result is {}.", b_bad);


    //let val:G1Projective = Bls12_381::inner_product(&vec![s, s], &vec![a, a]);



    /*

    let rng = rand::thread_rng();

    let mut pp = PP {l:1, t: 2, rng:rng};

    let m = Matrix::new(pp.l, pp.t, &vec![G1::one(), G1::one()]);

    let x:Vec<Fr> = vec![Fr::one(), Fr::zero()];

    let y:Vec<G1> = vec![G1::one()];

    let (ek, vk) = keygen(&mut pp, m);

    let pi = prove(&mut pp, &ek, &x);

    let b = verify(&pp, &vk, &y, &pi);

    println!("Result is {}.", b);

    */

}
