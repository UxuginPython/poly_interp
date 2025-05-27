use approx::assert_abs_diff_eq;
use poly_interp::*;
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
#[test]
fn from_zeros() {
    let polynomial = Polynomial::from_zeros(vec![1.0, 2.0, 3.0]);
    assert_eq!(polynomial, Polynomial::new(vec![-6.0, 11.0, -6.0, 1.0]));
    assert_eq!(polynomial.evaluate(1.0), PointXY::new(1.0, 0.0));
    assert_eq!(polynomial.evaluate(2.0), PointXY::new(2.0, 0.0));
    assert_eq!(polynomial.evaluate(3.0), PointXY::new(3.0, 0.0));
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
fn polynomial_newtons_method() {
    //x^2
    let polynomial = Polynomial::new(vec![0.0, 0.0, 1.0]);
    assert_eq!(
        polynomial.newtons_method(9.0, 1.0, 1000),
        PointXY::new(3.0, 9.0)
    );
}
#[test]
fn polynomial_from_f64() {
    let polynomial = Polynomial::from(5.0);
    assert_eq!(polynomial.evaluate(0.0), PointXY::new(0.0, 5.0));
    assert_eq!(polynomial.evaluate(-1.0), PointXY::new(-1.0, 5.0));
    assert_eq!(polynomial.evaluate(1.0), PointXY::new(1.0, 5.0));
    assert_eq!(polynomial.evaluate(1000.0), PointXY::new(1000.0, 5.0));
}
#[test]
fn polynomial_neg() {
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(-polynomial, Polynomial::new(vec![-1.0, -2.0, -3.0, -4.0]));
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
fn polynomial_sub() {
    let polynomial_a = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    let polynomial_b = Polynomial::new(vec![8.0, 7.0, 6.0, 5.0]);
    assert_eq!(
        polynomial_a - polynomial_b,
        Polynomial::new(vec![-7.0, -5.0, -3.0, -1.0])
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
fn xyt_curve() {
    let xyt_curve = XYTCurve::new(vec![
        PointXYT::new(2.0, 3.0, 5.0),
        PointXYT::new(7.0, 11.0, 13.0),
        PointXYT::new(17.0, 19.0, 23.0),
    ]);
    let point = xyt_curve.evaluate(5.0);
    assert_abs_diff_eq!(point.x, 2.0);
    assert_abs_diff_eq!(point.y, 3.0);
    assert_abs_diff_eq!(point.t, 5.0);
    let point = xyt_curve.evaluate(13.0);
    assert_abs_diff_eq!(point.x, 7.0);
    assert_abs_diff_eq!(point.y, 11.0, epsilon = 0.00000000000001);
    assert_abs_diff_eq!(point.t, 13.0);
    let point = xyt_curve.evaluate(23.0);
    assert_abs_diff_eq!(point.x, 17.0);
    assert_abs_diff_eq!(point.y, 19.0);
    assert_abs_diff_eq!(point.t, 23.0);
}
#[test]
fn xyt_curve_derivative() {
    //This should make x=t and y=t^2
    let xyt_curve = XYTCurve::new(vec![
        PointXYT::new(0.0, 0.0, 0.0),
        PointXYT::new(1.0, 1.0, 1.0),
        PointXYT::new(2.0, 4.0, 2.0),
    ]);
    //This should be x=1 and y=2t
    let derivative = xyt_curve.derivative();
    assert_eq!(
        derivative,
        XYTCurve::new(vec![
            PointXYT::new(1.0, 0.0, 0.0),
            PointXYT::new(1.0, 2.0, 1.0),
            PointXYT::new(1.0, 4.0, 2.0),
        ])
    );
}
