use std::ops::Add;
use std::ops::Sub;

pub type H3Resolution = u8;
pub type BaseCell = u8;

#[derive(Debug, PartialEq, PartialOrd)]
// Base cell at a given ijk
// and required rotations into its system
pub struct BaseCellOrientation {
    base_cell: BaseCell,
    // number of ccw 60 degree rotations relative to current face
    ccw_rotations: usize
}

#[derive(Debug, PartialEq, PartialOrd)]
pub enum Direction {
    Center = 0,
    KAxes = 1,
    JAxes = 2,
    JKAxes = 3,
    IAxes = 4,
    IKAxes = 5,
    IJAxes = 6,
    Invalid = 7
}

impl Direction {
    pub fn from_int(i: i64) -> Direction {
        match i {
            0 => Direction::Center,
            1 => Direction::KAxes,
            2 => Direction::JAxes,
            3 => Direction::JKAxes,
            4 => Direction::IAxes,
            5 => Direction::IKAxes,
            6 => Direction::IJAxes,
            _ => Direction::Invalid,
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
pub struct CoordIJK {
    pub i: i64,
    pub j: i64,
    pub k: i64
}

impl CoordIJK {
    pub fn new(i: i64, j: i64, k: i64) -> CoordIJK {
        CoordIJK{i: i, j: j, k: k}
    }

    pub fn scale(self, mult: i64) -> CoordIJK {
        Self::new(self.i * mult, self.j * mult, self.k * mult)
    }
}

impl Add for CoordIJK {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.i + other.i,
                  self.j + other.j,
                  self.k + other.k)
    }
}

impl Sub for CoordIJK {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.i - other.i,
                  self.j - other.j,
                  self.k - other.k)
    }
}
