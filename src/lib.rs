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
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Polynomial {
    coefficients: Vec<f64>,
}
impl Polynomial {
    #[inline]
    pub fn new(coefficients: Vec<f64>) -> Self {
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
        let mut zeroing_polynomials = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let mut new_zeros = Vec::with_capacity(i);
            for j in 0..i {
                new_zeros.push(points[j].x);
            }
            zeroing_polynomials.push(Self::from_zeros(new_zeros));
        }
        let mut coefficients = Vec::with_capacity(points.len());
        for i in 0..points.len() {
            let mut numerator = points[i].y;
            for j in 0..i {
                numerator =
                    numerator - zeroing_polynomials[j].evaluate(points[i].x).y * coefficients[j];
            }
            coefficients.push(numerator / zeroing_polynomials[i].evaluate(points[i].x).y);
        }
        let mut final_polynomial = Self::default();
        for (i, zeroing_polynomial) in zeroing_polynomials.into_iter().enumerate() {
            final_polynomial = final_polynomial + zeroing_polynomial * coefficients[i];
        }
        final_polynomial
    }
    pub fn evaluate(&self, x: f64) -> PointXY {
        let mut y = 0.0;
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            y += coefficient * x.powi(i as i32);
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
