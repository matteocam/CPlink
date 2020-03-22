//extern crate bn;
extern crate rand;

//use bn::*;
use rand::Rng;
use algebra_core::{PairingEngine,test_rng};
use algebra::{
    One,
    Zero,
    ProjectiveCurve,
    UniformRand,
    bls12_381::{
        g1, g2, Bls12_381, Fq, Fq12, Fq2, Fr, G1Affine, G1Projective, G2Affine, G2Projective,
    }
};

use crate::matrix::*;


//pub type VecG = Vec<G1>;


pub struct PP<G1: Clone,G2: Clone>
{
    pub     l: usize, // # of rows
    pub     t: usize, // # of cols
    pub     g1: G1,
    pub     g2: G2,
}

impl<G1: Clone,G2: Clone> PP<G1,G2> {
    pub fn new(l: usize, t: usize, g1: &G1, g2: &G2) -> PP<G1,G2> {
        PP {l:l, t:t, g1: g1.clone(), g2: g2.clone() }
    }
}

pub struct EK<G1> {
    p: Vec<G1>,
}
pub struct VK<G2> {
    c: Vec<G2>,
    a: G2,
}

pub trait SubspaceSnark
{

    type KMtx;
    type InVec;
    type OutVec;

    type PP;

    type EK;
    type VK;
    type Crs = (Self::EK, Self::VK);

    type Proof;

    fn keygen<R:Rng>(rng: &mut R, pp: &Self::PP, m: Self::KMtx) -> Self::Crs;
    fn prove(pp : &Self::PP, ek: &Self::EK, x: &Self::InVec) -> Self::Proof;
    fn verify(pp : &Self::PP, vk: &Self::VK, y: &Self::OutVec, pi: &Self::Proof) -> bool;

}


fn vec_to_G2<PE: PairingEngine>(pp: &PP<PE::G1Projective,PE::G2Projective>, v: & Vec<PE::Fr>) ->  Vec<PE::G2Projective>
{
    c![ pp.g2.mul(*x), for x in v]
}

// NB: Now the system is for y = Mx
impl<PE> SubspaceSnark for PE
    where PE: PairingEngine + SparseLinAlgebra<PE>
{
    
    type KMtx = SparseMatrix<PE::G1Projective>;
    type InVec = Vec<PE::Fr>;
    type OutVec = Vec<PE::G1Projective>;

    type PP = PP<PE::G1Projective, PE::G2Projective>;

    type EK = EK<PE::G1Projective>;
    type VK = VK<PE::G2Projective>;

    type Proof = PE::G1Projective;


    fn keygen<R:Rng>(rng: &mut R, pp: &Self::PP, m: Self::KMtx) -> Self::Crs
    {
        let mut k: Vec<PE::Fr>  = Vec::with_capacity(pp.l);
        for _ in 0..pp.l {
            k.push(PE::Fr::rand(rng));
        }

        let a = PE::Fr::rand(rng);
        let mut p: Vec<PE::G1Projective> = Vec::with_capacity(pp.t);

        PE::sparse_vector_matrix_mult(&k, &m, &mut p);
        let mut c: Vec<PE::Fr> = Vec::with_capacity(pp.l);

        scalar_vector_mult::<PE>(&a, &k, &mut c);
        let ek = EK::<PE::G1Projective> {p:p};
        let vk = VK::<PE::G2Projective> {c: vec_to_G2::<PE>(pp, &c), a: pp.g2.mul(a)};
        (ek, vk)
    }
    
    fn prove(pp : &Self::PP, ek: &Self::EK, x: &Self::InVec) -> Self::Proof {
        assert_eq!(pp.t, x.len());
        inner_product::<PE>(x, &ek.p) as Self::Proof
    }

    fn verify(pp : &Self::PP, vk: &Self::VK, y: &Self::OutVec, pi: &Self::Proof) -> bool {
        assert_eq!(pp.l, y.len());

        let mut res = PE::Fqk::one();
        for i in 0..y.len() {
            res *= &PE::pairing(y[i], vk.c[i]);
        }
        res == PE::pairing(*pi, vk.a)
    }
    
}



/*

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_snark() {
        /*
        let rng = rand::thread_rng();

        let mut pp = PP {l:2, t: 2, rng:rng};

        /*
            m =  1  2
                 1  1
        */
        let m = Matrix::new(pp.l, pp.t, &vec![G1::one(), G1::one(), G1::one()+G1::one(), G1::one()]);

        /*
            x =  [0, 1]
        */
        let x:Vec<Fr> = vec![Fr::zero(), Fr::one()];

        /*
            y =  [1, 2]
        */
        let y:VecG = vec![G1::one()+G1::one(), G1::one(),];

        let (ek, vk) = keygen(&mut pp, m);

        let pi = prove(&mut pp, &ek, &x);

        let b = verify(&pp, &vk, &y, &pi);

        assert_eq!(b, true);
        */

    }
}
*/

