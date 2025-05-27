use poly_interp::*;
#[test]
fn evaluate() {
    //4x^3+3x^2+2x+1
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(polynomial.evaluate(5.0), PointXY::new(5.0, 586.0));
}
#[test]
fn derivative() {
    //4x^3+3x^2+2x+1
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    //12x^2+6x+2
    assert_eq!(
        polynomial.derivative(),
        Polynomial::new(vec![2.0, 6.0, 12.0])
    );
}
#[test]
fn polynomial_add() {
    let polynomial_a = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    let polynomial_b = Polynomial::new(vec![5.0, 6.0, 7.0, 8.0]);
    assert_eq!(
        polynomial_a + polynomial_b,
        Polynomial::new(vec![6.0, 8.0, 10.0, 12.0])
    );
}
#[test]
fn polynomial_mul_f64() {
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(polynomial * 2.0, Polynomial::new(vec![2.0, 4.0, 6.0, 8.0]));
}
#[test]
fn polynomial_div_f64() {
    let polynomial = Polynomial::new(vec![2.0, 4.0, 6.0, 8.0]);
    assert_eq!(polynomial / 2.0, Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]));
}
