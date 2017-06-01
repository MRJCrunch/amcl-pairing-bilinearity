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
use self::amcl::pair::{ate, g1mul, g2mul, gtpow};
use self::amcl::rand::RAND;

extern crate rand;

use self::rand::os::{OsRng};
use self::rand::Rng;

fn main() {
    // Order of the groups
    let group_order = BIG::new_ints(&CURVE_ORDER);
    println!("{0: <10}: {1:?}", "order", group_order);

    // Seed generation
    let mut seed = vec![0; 96];
    let mut os_rng = OsRng::new().unwrap();
    os_rng.fill_bytes(&mut seed.as_mut_slice());
    let mut rng = RAND::new();
    rng.clean();
    rng.seed(96, &seed);

    // random element in 0, ..., group_order-1
    let mut a = BIG::randomnum(&group_order, &mut rng);
    println!("{0: <10}: {1:?}", "a", a);

    let mut b = BIG::randomnum(&group_order, &mut rng);
    println!("{0: <10}: {1:?}", "b", b);

    // Generator of Group G1
    let point_x = BIG::new_ints(&CURVE_GX);
    let point_y = BIG::new_ints(&CURVE_GY);
    let mut gen_g1 = ECP::new_bigs(&point_x, &point_y);
    println!("{0: <10}: {1:?}", "gen_g1", gen_g1);

    // Generator of Group G2
    let point_xa = BIG::new_ints(&CURVE_PXA);
    let point_xb = BIG::new_ints(&CURVE_PXB);
    let point_ya = BIG::new_ints(&CURVE_PYA);
    let point_yb = BIG::new_ints(&CURVE_PYB);
    let point_x = FP2::new_bigs(&point_xa, &point_xb);
    let point_y = FP2::new_bigs(&point_ya, &point_yb);
    let mut gen_g2 = ECP2::new_fp2s(&point_x, &point_y);
    println!("{0: <10}: {1:?}", "gen_g2", gen_g2);

    // Calculating left part
    let mut p = g1mul(&mut gen_g1, &mut a);
    let mut q = g2mul(&mut gen_g2, &mut b);
    let left = ate(&mut q, &mut p);

    // Calculating right part
    let mut right = ate(&mut gen_g2, &mut gen_g1);
    let mut ab = BIG::modmul(&mut a, &mut b, &group_order);
    println!("{0: <10}: {1:?}", "a * b", ab);
    right = gtpow(&mut right, &mut ab);

    println!("{0: <10}: {1:?}", "left", left);
    println!("{0: <10}: {1:?}", "right", right);

    assert_eq!(left, right);
}