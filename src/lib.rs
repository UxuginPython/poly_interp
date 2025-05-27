use core::cmp::max;
use core::ops::*;
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointXY {
    pub x: f64,
    pub y: f64,
}
impl PointXY {
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}
pub fn newtons_method<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(
    function: F,
    derivative: D,
    y: f64,
    x_guess: f64,
    iterations: u16,
) -> PointXY {
    let mut x = x_guess;
    for _ in 0..iterations {
        let function_evaluation_minus_y = function(x) - y;
        if function_evaluation_minus_y == 0.0 {
            break;
        }
        x -= function_evaluation_minus_y / derivative(x);
    }
    PointXY::new(x, function(x))
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Polynomial {
    coefficients: Vec<f64>,
}
impl Polynomial {
    #[inline]
    pub fn new(coefficients: Vec<f64>) -> Self {
        let mut coefficients = coefficients;
        while coefficients.last() == Some(&0.0) {
            coefficients.pop();
        }
        Self { coefficients }
    }
    pub fn from_zeros(zeros: Vec<f64>) -> Self {
        let mut new_self = Polynomial::new(vec![1.0]);
        for x in zeros {
            new_self = new_self * Polynomial::new(vec![-x, 1.0]);
        }
        new_self
    }
    pub fn interpolate(points: Vec<PointXY>) -> Self {
        //In the Newton polynomial form, coefficients are almost always seen with their
        //corresponding binomial factors, e.g., c3(x-x0)(x-x1)(x-x2). This generates these binomial
        //factors without needing to calculate the coefficients themselves yet.
        let mut zeroing_polynomials = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let mut new_zeros = Vec::with_capacity(i);
            for j in 0..i {
                new_zeros.push(points[j].x);
            }
            zeroing_polynomials.push(Self::from_zeros(new_zeros));
        }
        //Now we do actually calculate the coefficients. Here's an example showing the algebra used
        //to derive the math for this section:
        //y3 = c0 + c1(x3-x0) + c2(x3-x0)(x3-x1) + c3(x3-x0)(x3-x1)(x3-x2)
        //y3 - c0 - c1(x3-x0) - c2(x3-x0)(x3-x1) = c3(x3-x0)(x3-x1)(x3-x2)
        //[y3 - c0 - c1(x3-x0) - c2(x3-x0)(x3-x1)] / [(x3-x0)(x3-x1)(x3-x2)] = c3
        //Remember that we generated the binomial factors (for example (x3-x0)(x3-x1)(x3-x2))
        //before, so we new just fetch them from that Vec and they don't look very interesting
        //here.
        let mut coefficients = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let mut numerator = points[i].y;
            for j in 0..i {
                numerator -= zeroing_polynomials[j].evaluate(points[i].x).y * coefficients[j];
            }
            coefficients.push(numerator / zeroing_polynomials[i].evaluate(points[i].x).y);
        }
        //This is the final summation into, to continue our example,
        //y = c0 + c1(x-x0) + c2(x-x0)(x-x1) + c3(x-x0)(x-x1)(x-x2)
        let mut final_polynomial = Self::default();
        for (i, zeroing_polynomial) in zeroing_polynomials.into_iter().enumerate() {
            final_polynomial = final_polynomial + zeroing_polynomial * coefficients[i];
        }
        final_polynomial
    }
    pub fn evaluate(&self, x: f64) -> PointXY {
        let mut y = 0.0;
        let mut x_power = 1.0;
        for coefficient in &self.coefficients {
            y += coefficient * x_power;
            x_power *= x;
        }
        PointXY::new(x, y)
    }
    pub fn derivative(&self) -> Polynomial {
        let mut derivative_coefficients = Vec::with_capacity(self.coefficients.len() - 1);
        for (i, coefficient) in self.coefficients.iter().enumerate().skip(1) {
            derivative_coefficients.push(i as f64 * coefficient);
        }
        Self::new(derivative_coefficients)
    }
    pub fn integral(&self, c: f64) -> Polynomial {
        let mut integral_coefficients = Vec::with_capacity(self.coefficients.len() + 1);
        integral_coefficients.push(c);
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            integral_coefficients.push(coefficient / (i + 1) as f64)
        }
        Polynomial::new(integral_coefficients)
    }
    pub fn newtons_method(&self, y: f64, x_guess: f64, iterations: u16) -> PointXY {
        let derivative = self.derivative();
        newtons_method(
            |x| self.evaluate(x).y,
            |x| derivative.evaluate(x).y,
            y,
            x_guess,
            iterations,
        )
    }
}
impl From<f64> for Polynomial {
    fn from(was: f64) -> Self {
        Self::new(vec![was])
    }
}
impl Neg for Polynomial {
    type Output = Self;
    fn neg(self) -> Self {
        let mut coefficients = self.coefficients;
        for coefficient in &mut coefficients {
            *coefficient *= -1.0;
        }
        Self::new(coefficients)
    }
}
impl Add for Polynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut new_coefficients = vec![0.0; max(self.coefficients.len(), rhs.coefficients.len())];
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            new_coefficients[i] += coefficient;
        }
        for (i, coefficient) in rhs.coefficients.iter().enumerate() {
            new_coefficients[i] += coefficient;
        }
        Self::new(new_coefficients)
    }
}
impl Sub for Polynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}
impl Mul<f64> for Polynomial {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let mut coefficients = self.coefficients;
        for coefficient in &mut coefficients {
            *coefficient *= rhs;
        }
        Self::new(coefficients)
    }
}
impl Div<f64> for Polynomial {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        let mut coefficients = self.coefficients;
        for coefficient in &mut coefficients {
            *coefficient /= rhs;
        }
        Self::new(coefficients)
    }
}
impl Mul for Polynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut coefficients = vec![0.0; self.coefficients.len() + rhs.coefficients.len() - 1];
        for (my_exponent, my_coefficient) in self.coefficients.iter().enumerate() {
            for (rhs_exponent, rhs_coefficient) in rhs.coefficients.iter().enumerate() {
                coefficients[my_exponent + rhs_exponent] += my_coefficient * rhs_coefficient;
            }
        }
        Polynomial::new(coefficients)
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointXYT {
    pub x: f64,
    pub y: f64,
    pub t: f64,
}
impl PointXYT {
    #[inline]
    pub fn new(x: f64, y: f64, t: f64) -> Self {
        Self { x, y, t }
    }
    #[inline]
    pub fn xy(&self) -> PointXY {
        PointXY::new(self.x, self.y)
    }
    #[inline]
    pub fn tx(&self) -> PointXY {
        PointXY::new(self.t, self.x)
    }
    #[inline]
    pub fn ty(&self) -> PointXY {
        PointXY::new(self.t, self.y)
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct XYTCurve {
    x_polynomial: Polynomial,
    y_polynomial: Polynomial,
}
impl XYTCurve {
    pub fn new(points: Vec<PointXYT>) -> Self {
        let mut x_polynomial_points = Vec::with_capacity(points.len());
        for point in &points {
            x_polynomial_points.push(point.tx());
        }
        let mut y_polynomial_points = Vec::with_capacity(points.len());
        for point in points {
            y_polynomial_points.push(point.ty());
        }
        Self {
            x_polynomial: Polynomial::interpolate(x_polynomial_points),
            y_polynomial: Polynomial::interpolate(y_polynomial_points),
        }
    }
    #[inline]
    pub fn evaluate(&self, t: f64) -> PointXYT {
        PointXYT::new(
            self.x_polynomial.evaluate(t).y,
            self.y_polynomial.evaluate(t).y,
            t,
        )
    }
    #[inline]
    pub fn derivative(&self) -> Self {
        Self {
            x_polynomial: self.x_polynomial.derivative(),
            y_polynomial: self.y_polynomial.derivative(),
        }
    }
    pub fn newtons_method_x(&self, x: f64, t_guess: f64, iterations: u16) -> PointXYT {
        let newton_output = self.x_polynomial.newtons_method(x, t_guess, iterations);
        let x = newton_output.y;
        let t = newton_output.x;
        let y = self.y_polynomial.evaluate(t).y;
        PointXYT::new(x, y, t)
    }
    pub fn newtons_method_y(&self, y: f64, t_guess: f64, iterations: u16) -> PointXYT {
        let newton_output = self.y_polynomial.newtons_method(y, t_guess, iterations);
        let y = newton_output.y;
        let t = newton_output.x;
        let x = self.x_polynomial.evaluate(t).y;
        PointXYT::new(x, y, t)
    }
    fn t_to_speed_squared(&self, t: f64) -> f64 {
        let derivative = self.derivative();
        let x = derivative.x_polynomial.evaluate(t).y;
        let y = derivative.y_polynomial.evaluate(t).y;
        x.powi(2) + y.powi(2)
    }
    pub fn t_to_distance(&self, t: f64) -> f64 {
        //Integral of square root
        self.t_to_speed_squared(t).powf(1.5) / 1.5
    }
    pub fn newtons_method_distance_to_t(
        &self,
        distance: f64,
        t_guess: f64,
        iterations: u16,
    ) -> f64 {
        newtons_method(
            |t| self.t_to_distance(t),
            |t| self.t_to_speed_squared(t).sqrt(),
            distance,
            t_guess,
            iterations,
        )
        .y
    }
}
