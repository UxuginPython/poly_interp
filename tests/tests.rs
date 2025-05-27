use poly_interp::*;
#[test]
fn evaluate() {
    //4x^3+3x^2+2x+1
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(polynomial.evaluate(5.0), PointXY::new(5.0, 586.0));
}
#[test]
fn from_zeros() {
    let polynomial = Polynomial::from_zeros(vec![1.0, 2.0, 3.0]);
    assert_eq!(polynomial, Polynomial::new(vec![-6.0, 11.0, -6.0, 1.0]));
    assert_eq!(polynomial.evaluate(1.0), PointXY::new(1.0, 0.0));
    assert_eq!(polynomial.evaluate(2.0), PointXY::new(2.0, 0.0));
    assert_eq!(polynomial.evaluate(3.0), PointXY::new(3.0, 0.0));
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
fn integral() {
    //4x^3+3x^2+2x+1
    let polynomial = Polynomial::new(vec![4.0, 3.0, 3.0, 1.0]);
    //0.25x^4+x^3+1.5x^2+4x+5
    assert_eq!(
        polynomial.integral(5.0),
        Polynomial::new(vec![5.0, 4.0, 1.5, 1.0, 0.25])
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
#[test]
fn polynomial_mul() {
    let polynomial_a = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    let polynomial_b = Polynomial::new(vec![5.0, 6.0, 7.0, 8.0]);
    assert_eq!(
        polynomial_a * polynomial_b,
        Polynomial::new(vec![5.0, 16.0, 34.0, 60.0, 61.0, 52.0, 32.0])
    );
}
#[test]
fn interpolate() {
    let polynomial = Polynomial::interpolate(vec![
        PointXY::new(2.0, 3.0),
        PointXY::new(5.0, 7.0),
        PointXY::new(11.0, 13.0),
    ]);
    assert_eq!(polynomial.evaluate(2.0), PointXY::new(2.0, 3.0));
    assert_eq!(polynomial.evaluate(5.0), PointXY::new(5.0, 7.0));
    assert_eq!(polynomial.evaluate(11.0), PointXY::new(11.0, 13.0));
}
#[test]
fn newtons_method_() {
    fn square(x: f64) -> f64 {
        x.powi(2)
    }
    //Derivative of square
    fn double(x: f64) -> f64 {
        2.0 * x
    }
    assert_eq!(
        newtons_method(square, double, 9.0, 1.0, 1000),
        PointXY::new(3.0, 9.0)
    );
}
