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
