extern crate rand;

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
use crate::subspace_snark::*;

extern crate rug;

use rug::Integer;
use crate::utils::integer_to_bigint;



pub struct PedCommOut<G,F>
{
    c: G,  // commitment
    vs: Vec<F>, // committed values
    o: F, // opening
}

pub struct LinkingPP<SubspacePP>
{
    n: usize,
    ss_pp: SubspacePP
}

/* Links _two_ commitments; could be generalized to more. */
pub trait LinkingSnark
{
    type G;
    type F;

    type LinkingPP;

    type LinkingEK;
    type LinkingVK;
    type LinkingCrs = (LinkingEK, LinkingVK);

    type LinkingProof;

    fn lnk_keygen<R:Rng>(rng: &mut R, pp: &Self::LinkingPP, h0: &G, h: &Vec<G>, f0 :&G, f: &Vec<G>) -> Self::LinkingCrs;
    fn lnk_prove(pp : &Self::LinkingPP, ek: &Self::LinkingEK, comm_out: &PedCommOut<G,F>, comm_out_prime: &PedCommOut<G,F>) -> Self::LinkingProof;
    fn lnk_verify(pp : &Self::LinkingPP, vk: &Self::LinkingVK, c: &G, c_prime: &G, pf: &Self::LinkingProof) -> bool;

}

// XXX: *** Change type parametrization. Have it on parameters or something similar ***

impl<PE> LinkingSnark for PE
    where PE: PairingEngine  + SubspaceSnark  + SparseLinAlgebra<PE>
{
    type G = PE::G1Projective;
    type F = PE::Fr;
    type LinkingPP = LinkingPP<PE::PP>;

    type LinkingEK = PE::EK;
    type LinkingVK = PE::VK;

    type LinkingProof = PE::Proof;

    fn lnk_keygen<R:Rng>(rng: &mut R, pp: &Self::LinkingPP, n:usize, h0: &G, h: &Vec<G>, f0:&G, f: &Vec<G>) -> Self::LinkingCrs
    {
        let n = pp.n;
        
        assert_eq!(h.len(), n);
        assert_eq!(f.len(), n);

        let mut m = SparseMatrix::new(2, 2+2*n);
        m.insert_val(0, 0, h0);
        m.insert_row_slice(0, 2, h);

        m.insert_val(1, 1, f0);
        m.insert_row_slice(0, 2, f);
    
        let (ek, vk) = PE::keygen(&mut rng, &pp.ss_pp, &m);
    }

    fn lnk_prove(pp : &Self::LinkingPP, ek: &Self::LinkingEK, comm_out: &PedCommOut<G, F>, comm_out_prime: &PedCommOut<G,F>) -> Self::LinkingProof
    {
        assert_eq!(comm_out.vs, comm_out_prime.vs);

        let v_in = vec![comm_out.o, comm_out_prime.o].concat(comm_out.vs);
        
        let pi = PE::prove(&pp.ss_pp, ek, &v_in);
        pi
    }

    fn lnk_verify(pp : &Self::LinkingPP, vk: &Self::LinkingVK, c: &G, c_prime: &G, pf: &Self::LinkingProof) -> bool
    {
        let v_out:Vec<PE::G1Projective> = vec![c, c_prime];
        PE::verify(&pp.ss_pp, vk, &v_out, pf)
    }

}
