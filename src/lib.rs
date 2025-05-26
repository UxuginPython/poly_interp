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
