pub type BaseCell = u8;

#[derive(Debug, PartialEq, PartialOrd)]
// Base cell at a given ijk
// and required rotations into its system
pub struct BaseCellOrientation {
    base_cell: BaseCell,
    // number of ccw 60 degree rotations relative to current face
    ccw_rotations: usize
}
