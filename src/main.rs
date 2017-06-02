extern crate amcl;

use self::amcl::big::BIG;
use self::amcl::rom::{
    CURVE_GX,
    CURVE_GY,
    CURVE_ORDER,
    CURVE_PXA,
    CURVE_PYA,
    CURVE_PXB,
    CURVE_PYB
};
use self::amcl::ecp::ECP;
use self::amcl::ecp2::ECP2;
use self::amcl::fp2::FP2;
use self::amcl::pair::{ate, g1mul, g2mul, gtpow, fexp};
use self::amcl::rand::RAND;

extern crate rand;

use self::rand::os::{OsRng};
use self::rand::Rng;

fn get_group_order() -> BIG {
    BIG::new_ints(&CURVE_ORDER)
}

fn get_rng() -> RAND {
    // Seed generation
    let mut seed = vec![0; 96];
    let mut os_rng = OsRng::new().unwrap();
    os_rng.fill_bytes(&mut seed.as_mut_slice());
    let mut rng = RAND::new();
    rng.clean();
    rng.seed(96, &seed);
    rng
}

fn get_random_element_from_gr() -> BIG {
    BIG::randomnum(&get_group_order(), &mut get_rng())
}

fn get_generator_of_g1() -> ECP {
    let point_x = BIG::new_ints(&CURVE_GX);
    let point_y = BIG::new_ints(&CURVE_GY);
    ECP::new_bigs(&point_x, &point_y)
}

fn get_generator_of_g2() -> ECP2 {
    let point_xa = BIG::new_ints(&CURVE_PXA);
    let point_xb = BIG::new_ints(&CURVE_PXB);
    let point_ya = BIG::new_ints(&CURVE_PYA);
    let point_yb = BIG::new_ints(&CURVE_PYB);
    let point_x = FP2::new_bigs(&point_xa, &point_xb);
    let point_y = FP2::new_bigs(&point_ya, &point_yb);
    ECP2::new_fp2s(&point_x, &point_y)
}

fn get_point_from_g1() -> ECP {
    let mut gen_g1 = get_generator_of_g1();
    let mut a = get_random_element_from_gr();
    g1mul(&mut gen_g1, &mut a)
}

fn get_point_from_g2() -> ECP2 {
    let mut gen_g2 = get_generator_of_g2();
    let mut a = get_random_element_from_gr();
    g2mul(&mut gen_g2, &mut a)
}

fn main() {
    /// Order of the groups
    let group_order = get_group_order();
    println!("{0: <10}: {1:?}", "order", group_order);

    /// random element in 0, ..., group_order-1
    let mut a = get_random_element_from_gr();
    println!("{0: <10}: {1:?}", "a", a);

    let mut b = get_random_element_from_gr();
    println!("{0: <10}: {1:?}", "b", b);

    /// Generator of Group G1
    let mut gen_g1 = get_generator_of_g1();
    println!("{0: <10}: {1:?}", "gen_g1", gen_g1);

    /// Generator of Group G2
    let mut gen_g2 = get_generator_of_g2();
    println!("{0: <10}: {1:?}", "gen_g2", gen_g2);

    /// Calculating left part
    let mut p = g1mul(&mut gen_g1, &mut a);
    let mut q = g2mul(&mut gen_g2, &mut b);
    let mut left = ate(&mut q, &mut p);
    left = fexp(&left);

    /// Calculating right part
    let mut right = ate(&mut gen_g2, &mut gen_g1);
    right = fexp(&right);
    let mut ab = BIG::modmul(&mut a, &mut b, &group_order);
    println!("{0: <10}: {1:?}", "a * b", ab);
    right = gtpow(&mut right, &mut ab);

    println!("{0: <10}: {1:?}", "left", left);
    println!("{0: <10}: {1:?}", "right", right);

    /// check bilinearity
    assert_eq!(left, right);


    /// check infinity point
    let mut p = ECP::new();
    p.inf();
    p.add(&mut gen_g1);
    assert_eq!(p, gen_g1);

    /// check inverse for fp12(result of pairing)
    let mut p1 = get_point_from_g1();
    let mut q1 = get_point_from_g2();
    let mut pair1 = ate(&mut q1, &mut p1);
    pair1 = fexp(&pair1);
    println!("{0: <10}: {1:?}", "pair1", pair1);

    let mut p2 = get_point_from_g1();
    let mut q2 = get_point_from_g2();
    let mut pair2 = ate(&mut q2, &mut p2);
    pair2 = fexp(&pair2);
    println!("{0: <10}: {1:?}", "pair2", pair2);

    let mut result = pair1;
    result.mul(&mut pair2);
    println!("{0: <10}: {1:?}", "result", result);

    pair1.inverse();
    result.mul(&mut pair1);
    println!("{0: <10}: {1:?}", "result", result);
    
    pair2.reduce(); result.reduce();
    assert_eq!(pair2, result);
}