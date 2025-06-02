# poly_interp

**A simple but powerful polynomial library focused on interpolation between points.** Works with both standard polynomials and parametric curves. Mostly no_std although a dynamic allocator is required. Based on the f64 number type. BSD licensed.

## Features
- Construct polynomials from coefficients or zeros.
- Interpolate polynomials of least degree between points.
- Solve for derivatives and integrals of polynomials.
- Add, subtract, and multiply polynomials.
- Use Newton's method of iteratively finding roots of functions generically or with the included types.
- Construct polynomial parametric curves from (x, y, t) points that they pass through.
  - Convert between *t* and distance along the curve (std-only).
