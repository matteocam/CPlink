use std::ops::{Mul,AddAssign};
use algebra_core::PairingEngine;
use algebra::ProjectiveCurve;
use algebra::Zero;


/*
    CoeffPos: A struct to help build sparse matrices.
*/
#[derive(Clone, Debug)]
struct CoeffPos<T>
{
    val:T,
    pos: usize
}

// a column is a vector of CoeffPos-s
type Col<T> = Vec<CoeffPos<T>>;


/* TODO: One could consider a cache-friendlier implementation for the 2-row case*/

/* Column-Major Sparse Matrix */
#[derive(Clone, Debug)]
pub struct SparseMatrix<T>
{
    cols: Vec<Col<T>>, // a vector of columns
    pub    nr: usize,
    pub    nc: usize,
    //pub    N: usize
}

impl<T: Copy> SparseMatrix<T> {

    // NB: Given column by column
    pub fn new(nr:usize, nc:usize) -> SparseMatrix<T>
    {
        SparseMatrix {
            cols: vec![vec![]; nc],
            nr: nr,
            nc: nc
        }
    }

    pub fn insert_val(&mut self, r: usize, c: usize, v: &T)
    {
        let coeff_pos = CoeffPos {pos: r, val: *v};
        self.cols[c].push(coeff_pos);
    } 

    // insert a continguous sequence of values at row r starting from c_offset
    pub fn insert_row_slice(&mut self, r: usize, c_offset: usize, vs: &Vec<T>)
    {
        // NB: could be improved in efficienc by first extending the vector
        for (i, x) in vs.iter().enumerate() {
            self.insert_val(r, c_offset+i, x);
        }
    }

    pub fn get_col(&self, c: usize) -> &Col<T>
    {
        &self.cols[c]
    }

}

pub trait SparseLinAlgebra<PE: PairingEngine>
{
    // this is basically a multi-exp
    fn sparse_inner_product(v: &Vec<PE::Fr>, w: &Col<PE::G1Projective>) -> PE::G1Projective
    {
        assert_eq!(v.len(), w.len());
        let mut res:PE::G1Projective = PE::G1Projective::zero();
        for coeffpos in w {
            let g = coeffpos.val;
            let i = coeffpos.pos;
            // XXX: Should this be optimized for special cases 
            //         (e.g. 0 or 1) or is this already in .mul?
            let tmp = g.mul(v[i]);
            res.add_assign(&tmp);
        }
        res
    }

    fn sparse_vector_matrix_mult(v: &Vec<PE::Fr>, m:&SparseMatrix<PE::G1Projective>, res: &mut Vec<PE::G1Projective>) {
        // the result should contain every column of m multiplied by v
        for c in 0..m.nc {
            res.push(
                Self::sparse_inner_product(
                    &v, &m.get_col(c)));
        }
    
    }


}

pub fn inner_product<PE: PairingEngine>(v: &Vec<PE::Fr>, w: &Vec<PE::G1Projective>) -> PE::G1Projective
{
    assert_eq!(v.len(), w.len());
    let mut res:PE::G1Projective = PE::G1Projective::zero();
    for i in 0..v.len() {
        let tmp = w[i].mul(v[i]);
        res.add_assign(&tmp);
    }
    res
}

pub fn scalar_vector_mult<PE: PairingEngine>(a: &PE::Fr, v: &Vec<PE::Fr>, res: &mut Vec<PE::Fr>)
{
    for i in 0..v.len() {
        let x:PE::Fr = a.mul(&v[i]);
        res.push(x);
    }
}

/*




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_col() {
        let m = Matrix::new(3, 2, &vec![(0,0), (1,0), (2,0), (0,1),  (1,1), (2,1) ]);
        //println!("{:?}", m.get_col(0));

        assert_eq!(m.get_col(0), vec![(0,0), (1,0), (2,0)]);
    }

    #[test]
    fn test_inner_product_nat() {
        let u = vec![0,1,2];
        let v = vec![1, 2, 3];
        let y = inner_product(&u,&v,0);
        assert_eq!(y, 8);
    }

    fn testEq(a:&G1, b:&G1) -> bool {
        (*a-*b).is_zero()
    }

    #[test]
    fn test_inner_product_G1() {
        let u = vec![G1::one(), G1::zero()];
        let v:Vec<Fr> = vec![Fr::one()+Fr::one(), Fr::zero()];
        let y = inner_product(&v,&u, G1::zero());

        assert!(y == G1::one()+G1::one());
    }

    #[test]
    fn test_scalar_vector_mult() {
        let a = 10;
        let u:Vec<i32> = (1..4).collect();
        let mut y = Vec::with_capacity(u.len());
        scalar_vector_mult(&a, &u, &mut y);
        assert_eq!(y, [10,20,30]);
    }

}

*/