// SPDX-License-Identifier: BSD-3-Clause
// Copyright 2025 UxuginPython
//!A simple but powerful polynomial library focused on interpolation between points. Works with
//!both standard polynomials and parametric curves. Mostly no_std although a dynamic allocator is
//!required.
#![warn(missing_docs)]
#![cfg_attr(not(feature = "std"), no_std)]
extern crate alloc;
use alloc::{vec, vec::Vec};
use core::cmp::max;
use core::ops::*;
///A point with x and y coordinates.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointXY {
    ///The x coordinate.
    pub x: f64,
    ///The y coordinate.
    pub y: f64,
}
impl PointXY {
    ///Constructor for `PointXY`.
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}
///Performs an iterative process to approximate the inverse of a function using its derivative at a
///given y coordinate. (Solves f(x)=y for x at a given y.) This requires an initial "guess" at the
///desired x value and a maximum number of iterations to perform when calculating it. The function
///returns early if it detects that x has been calculated perfectly (`function(x)` returns exactly
///`y`). It may return NaN or infinity if the derivative returns 0 at one of the progressive
///estimates of x. The function returns (x, f(x)) and not (x, y) if they are different.
pub fn newtons_method<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(
    function: F,
    derivative: D,
    y: f64,
    x_guess: f64,
    max_iterations: u16,
) -> PointXY {
    let mut x = x_guess;
    for _ in 0..max_iterations {
        let function_evaluation_minus_y = function(x) - y;
        if function_evaluation_minus_y == 0.0 {
            break;
        }
        x -= function_evaluation_minus_y / derivative(x);
    }
    PointXY::new(x, function(x))
}
///A polynomial function.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Polynomial {
    coefficients: Vec<f64>,
}
impl Polynomial {
    ///Constructor for `Polynomial`. The index of each coefficient in the `Vec` is its
    ///corresponding exponent. For example, the polynomial x^4+2x^2+3x+1 would be constructed with
    ///`Polynomial::new(vec![1.0, 3.0, 2.0, 0.0, 1.0])`.
    #[inline]
    pub fn new(coefficients: Vec<f64>) -> Self {
        let mut coefficients = coefficients;
        while coefficients.last() == Some(&0.0) {
            coefficients.pop();
        }
        Self { coefficients }
    }
    ///Construct a `Polynomial` where `polynomial.evaluate(x) == PointXY::new(x, 0.0)` for every
    ///value of `x` in the `zeros` `Vec`.
    pub fn from_zeros(zeros: Vec<f64>) -> Self {
        let mut new_self = Polynomial::new(vec![1.0]);
        for x in zeros {
            new_self = new_self * Polynomial::new(vec![-x, 1.0]);
        }
        new_self
    }
    ///Interpolate the lowest-degree polynomial going through every point in the `points` `Vec`
    ///using an algorithm based on Newton's form.
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
    ///Evaluate the polynomial. Returns a `PointXY` containing (x, f(x)) where x is `evaluate`'s
    ///input and f(x) is the mathematical function the `Polynomial` object represents.
    pub fn evaluate(&self, x: f64) -> PointXY {
        let mut y = 0.0;
        let mut x_power = 1.0;
        for coefficient in &self.coefficients {
            y += coefficient * x_power;
            x_power *= x;
        }
        PointXY::new(x, y)
    }
    ///Solve for the derivative function of the polynomial using the power rule.
    pub fn derivative(&self) -> Polynomial {
        let mut derivative_coefficients = Vec::with_capacity(self.coefficients.len() - 1);
        for (i, coefficient) in self.coefficients.iter().enumerate().skip(1) {
            derivative_coefficients.push(i as f64 * coefficient);
        }
        Self::new(derivative_coefficients)
    }
    ///Calculate the indefinite integral of the polynomial using the power rule given c, the
    ///constant of integration.
    pub fn integral(&self, c: f64) -> Polynomial {
        let mut integral_coefficients = Vec::with_capacity(self.coefficients.len() + 1);
        integral_coefficients.push(c);
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            integral_coefficients.push(coefficient / (i + 1) as f64)
        }
        Polynomial::new(integral_coefficients)
    }
    ///Estimate an x for which `polynomial.evaluate(x)` will return (x, y) for a given y using
    ///Newton's method. This requires an initial "guess" at this x value and a maximum number of
    ///iterations to perform when estimating x. See the documentation for [`newtons_method`], the
    ///function which this uses internally, for more information.
    pub fn newtons_method(&self, y: f64, x_guess: f64, max_iterations: u16) -> PointXY {
        let derivative = self.derivative();
        newtons_method(
            |x| self.evaluate(x).y,
            |x| derivative.evaluate(x).y,
            y,
            x_guess,
            max_iterations,
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
///A point with x, y, and t coordinates. t is usually the independent variable here.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PointXYT {
    ///The x coordinate.
    pub x: f64,
    ///The y coordinate.
    pub y: f64,
    ///The t coordinate. This is usually the independent variable.
    pub t: f64,
}
impl PointXYT {
    ///Constructor for `PointXYT`.
    #[inline]
    pub fn new(x: f64, y: f64, t: f64) -> Self {
        Self { x, y, t }
    }
    ///Get the point (x, y) as a `PointXY`.
    #[inline]
    pub fn xy(&self) -> PointXY {
        PointXY::new(self.x, self.y)
    }
    ///Get the point (t, x) as a `PointXY`.
    #[inline]
    pub fn tx(&self) -> PointXY {
        PointXY::new(self.t, self.x)
    }
    ///Get the point (t, y) as a `PointXY`.
    #[inline]
    pub fn ty(&self) -> PointXY {
        PointXY::new(self.t, self.y)
    }
}
///A smooth curve based on functions x(t) and y(t).
#[derive(Clone, Debug, PartialEq)]
pub struct XYTCurve {
    x_polynomial: Polynomial,
    y_polynomial: Polynomial,
}
impl XYTCurve {
    ///Interpolate an `XYTCurve` that will go through a set of points using lowest-degree
    ///polynomials for the internal x(t) and y(t) functions.
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
    ///Evaluate the curve at a given t. Returns (x(t), y(t), t) as a `PointXYT`.
    #[inline]
    pub fn evaluate(&self, t: f64) -> PointXYT {
        PointXYT::new(
            self.x_polynomial.evaluate(t).y,
            self.y_polynomial.evaluate(t).y,
            t,
        )
    }
    ///Return the derivative of the curve. This can be thought of as the "velocity" at which a
    ///point moves down the curve as t increases at a constant rate.
    #[inline]
    pub fn derivative(&self) -> Self {
        Self {
            x_polynomial: self.x_polynomial.derivative(),
            y_polynomial: self.y_polynomial.derivative(),
        }
    }
    ///Estimate a t at which `curve.evaluate(t)` will return (x, y(t), t) for a given x using
    ///Newton's method. This requires a "guess" at this t value and a maximum number of iterations
    ///to perform when estimating t. See the documentation of [`newtons_method`], the function
    ///which this calls internally, for more information.
    pub fn newtons_method_x(&self, x: f64, t_guess: f64, max_iterations: u16) -> PointXYT {
        let newton_output = self.x_polynomial.newtons_method(x, t_guess, max_iterations);
        let x = newton_output.y;
        let t = newton_output.x;
        let y = self.y_polynomial.evaluate(t).y;
        PointXYT::new(x, y, t)
    }
    ///Estimate a t at which `curve.evaluate(t)` will return (x(t), y, t) for a given y using
    ///Newton's method. This requires a "guess" at this t value and a maximum number of iterations
    ///to perform when estimating t. See the documentation of [`newtons_method`], the function
    ///which this calls internally, for more information.
    pub fn newtons_method_y(&self, y: f64, t_guess: f64, max_iterations: u16) -> PointXYT {
        let newton_output = self.y_polynomial.newtons_method(y, t_guess, max_iterations);
        let y = newton_output.y;
        let t = newton_output.x;
        let x = self.x_polynomial.evaluate(t).y;
        PointXYT::new(x, y, t)
    }
    //This would work in no_std if the last line was changed to `x * x + y * y`; it's just that
    //it's a private method and it's not used outside of std-only functions.
    #[cfg(feature = "std")]
    fn t_to_speed_squared(&self, t: f64) -> f64 {
        let derivative = self.derivative();
        let x = derivative.x_polynomial.evaluate(t).y;
        let y = derivative.y_polynomial.evaluate(t).y;
        x.powi(2) + y.powi(2)
    }
    ///Calculates the distance along the curve between 0 and t.
    #[cfg(feature = "std")]
    pub fn t_to_distance(&self, t: f64) -> f64 {
        //Integral of square root
        self.t_to_speed_squared(t).powf(1.5) / 1.5
    }
    ///Uses Newton's method to calculate t at a given distance along the curve. This requires an
    ///initial "guess" at this t and a maximum number of iterations to perform when estimating it.
    ///See the documentation of [`newtons_method`], the function which this calls internally, for
    ///more information.
    #[cfg(feature = "std")]
    pub fn newtons_method_distance_to_t(
        &self,
        distance: f64,
        t_guess: f64,
        max_iterations: u16,
    ) -> f64 {
        newtons_method(
            |t| self.t_to_distance(t),
            |t| self.t_to_speed_squared(t).sqrt(),
            distance,
            t_guess,
            max_iterations,
        )
        .y
    }
}
