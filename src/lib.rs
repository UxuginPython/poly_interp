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
    //TEMPORARY
    pub fn interpolate3(p0: PointXY, p1: PointXY, p2: PointXY) -> Self {
        //y0 = c0
        //y1 = c0 + c1(x1 - x0)
        //y1 - c0 = c1(x1 - x0)
        //y = c0 + c1(x - x0)
        let c0 = p0.y;
        let c1 = (p1.y - c0) / (p1.x - p0.x);
        let c2 = (p2.y - c0 - c1 * (p2.x - p0.x)) / (p2.x - p0.x) / (p2.x - p1.x);
        //y = c0 + c1(x - x0) + c2(x - x0)(x - x1)
        Self::from(c0)
            + Self::from(c1) * Self::from_zeros(vec![p0.x])
            + Self::from(c2) * Self::from_zeros(vec![p0.x, p1.x])
    }
    //TEMPORARY
    pub fn interpolate3_2(p0: PointXY, p1: PointXY, p2: PointXY) -> Self {
        let c0 = p0.y;
        let c1 = (p1.y - c0) / Self::from_zeros(vec![p0.x]).evaluate(p1.x).y;
        let c2 = (p2.y - c0 - c1 * Self::from_zeros(vec![p0.x]).evaluate(p2.x).y)
            / Self::from_zeros(vec![p0.x, p1.x]).evaluate(p2.x).y;
        Self::from(c0)
            + Self::from(c1) * Self::from_zeros(vec![p0.x])
            + Self::from(c2) * Self::from_zeros(vec![p0.x, p1.x])
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
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NewtonPolynomial3 {
    point0: PointXY,
    point1: PointXY,
    point2: PointXY,
}
impl NewtonPolynomial3 {
    pub fn new(point0: PointXY, point1: PointXY, point2: PointXY) -> Self {
        Self {
            point0,
            point1,
            point2,
        }
    }
    pub fn get(&self, x: f64) -> PointXY {
        /*y0=c0
        y1=c0+c1(x1-x0)
        y2=c0+c1(x2-x0)+c2(x2-x0)(x2-x1)*/
        let c0 = self.point0.y;
        let c1 = (self.point1.y - c0) / (self.point1.x - self.point0.x);
        let c2 = (self.point2.y - c0 - c1 * (self.point2.x - self.point0.x))
            / ((self.point2.x - self.point0.x) * (self.point2.x - self.point1.x));
        let y = c0 + c1 * (x - self.point0.x) + c2 * (x - self.point0.x) * (x - self.point1.x);
        PointXY::new(x, y)
    }
}
