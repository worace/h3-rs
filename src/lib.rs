#![allow(dead_code)]
extern crate geo_types;
use std::fmt;
use std::f64::consts::PI;
use std::cmp::min;

// use geo_types::Coordinate;

const TWO_PI: f64 = PI * 2.0;
mod constants;
mod types;

use types::*;

// [ ] geoToH3(coord, res) -> h3 id
// [X] geoToFace(coord) -> face (numeric id)
// [x] geoToFaceIJK(coord, res) -> faceIJK coord
// [X] geoToHex2d(coord, res) -> (face, vec2D)
// [x] geoToVec3d(coord) -> 3dCoord
// [X] hex2dToCoordIJK(vec2D) -> CoordIJK
// [X] pointSquareDist(vec3d, vec3d)
// Most coord ops in radians

#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
pub struct FaceIJK {
    face: usize,
    coord: CoordIJK
}
impl FaceIJK {
    fn new(face: usize, i: i64, j: i64, k: i64) -> FaceIJK {
        FaceIJK{face: face,
                coord: CoordIJK{i: i, j: j, k: k}}
    }
}

#[derive(Debug)]
pub struct GeoCoord {
    lat: f64,
    lon: f64
}

impl fmt::Display for GeoCoord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "GeoCoord(lat: {}, lon: {})", self.lat, self.lon)
    }
}

impl GeoCoord {
    fn new(lat: f64, lon: f64) -> GeoCoord {
        GeoCoord{lat: lat, lon: lon}
    }
}

#[derive(Debug, PartialEq, PartialOrd)]
pub struct Vec2d {
    x: f64,
    y: f64
}

impl Vec2d {
    fn zero() -> Vec2d { Vec2d{x: 0.0, y: 0.0} }
    fn new(x: f64, y: f64) -> Vec2d { Vec2d{x: x, y: y} }
}

#[derive(Debug, PartialEq, PartialOrd)]
struct FaceCoord2d {
    face: usize,
    coords: Vec2d
}

impl fmt::Display for FaceCoord2d {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "FaceCoord2d(face: {}, coords: ({}, {}))",
               self.face, self.coords.x, self.coords.y)
    }
}

#[derive(Debug, PartialEq, PartialOrd)]
pub struct Vec3d {
    x: f64,
    y: f64,
    z: f64
}

impl Vec3d {
    fn new(x: f64, y: f64, z: f64) -> Vec3d {
        Vec3d{x: x, y: y, z: z}
    }
}

fn square(n: f64) -> f64 {
    n.powi(2)
}

fn square_distance_3d(a: &Vec3d, b: &Vec3d) -> f64 {
    square(a.x - b.x) + square(a.y - b.y) + square(a.z - b.z)
}

fn nearest_face_to_geo(geo: &GeoCoord) -> (usize, f64) {
    let v3d = geo_to_coord_3d(geo);
    let mut face = 0;
    let mut min_dist = square_distance_3d(&v3d, &constants::FACE_CENTERS[0]);
    for i in 1..constants::NUM_ICOSA_FACES {
        let dist = square_distance_3d(&v3d, &constants::FACE_CENTERS[i]);
        if dist < min_dist {
            face = i;
            min_dist = dist;
        }
    }
    (face, min_dist)
}

fn geo_to_coord_3d(geo: &GeoCoord) -> Vec3d {
    let r = geo.lat.cos();

    Vec3d::new(
        geo.lon.cos() * r,
        geo.lon.sin() * r,
        geo.lat.sin()
    )
}

fn geo_azimuth_radians(a: &GeoCoord, b: &GeoCoord) -> f64 {
    let y = b.lat.cos() * (b.lon - a.lon).sin();
    let x = (a.lat.cos() * b.lat.sin()) -
        (a.lat.sin() * b.lat.cos() * (b.lon - a.lon).cos());
    y.atan2(x)
}

/**
 * Normalizes radians to a value between 0.0 and two PI.
 *
 * @param rads The input radians value.
 * @return The normalized radians value.
 */
fn positive_angle_radians(angle: f64) -> f64 {
    if angle < 0.0 {
        angle + TWO_PI
    } else if angle > TWO_PI {
        angle - TWO_PI
    } else {
        angle
    }
}

fn face_theta(geo: &GeoCoord, face: usize) -> f64 {
    let face_axis_i_azimuth = constants::CLASS_II_FACE_IJK_AXES[face][0];
    let face_center = &constants::FACE_CENTER_GEO_COORDS[face];
    let point_face_azimuth = positive_angle_radians(geo_azimuth_radians(face_center, geo));
    positive_angle_radians(face_axis_i_azimuth - point_face_azimuth)
}

/**
 * Returns whether or not a resolution is a Class III grid. Note that odd
 * resolutions are Class III and even resolutions are Class II.
 * @param res The H3 resolution.
 * @return 1 if the resolution is a Class III grid, and 0 if the resolution is
 *         a Class II grid.
 */
fn is_class_iii_resolution(res: H3Resolution) -> bool {
    res % 2 == 1
}

fn is_class_ii_resolution(res: H3Resolution) -> bool {
    res % 2 == 0
}

fn geo_to_hex_2d(geo: &GeoCoord, res: H3Resolution) -> FaceCoord2d {
    let (face, dist) = nearest_face_to_geo(geo);

    let r = (1.0 - dist / 2.0).acos();
    if r < constants::FACE_DISTANCE_EPSILON {
        return FaceCoord2d{face: face, coords: Vec2d::zero()};
    }

    let face_angle = if is_class_iii_resolution(res) {
        positive_angle_radians(
            face_theta(geo, face) - constants::AP7_ROT_RADS
        )
    } else {
        face_theta(geo, face)
    };


    let mut gnomonic_r = r.tan() / constants::RES0_U_GNOMONIC;
    for _ in 0..res {
        gnomonic_r *= constants::SQRT7;
    }

    let coords = Vec2d{x: gnomonic_r * face_angle.cos(),
                       y: gnomonic_r * face_angle.sin()};
    FaceCoord2d{face: face, coords: coords}
}

fn unnormalized_coord_ijk(h2d: &Vec2d) -> CoordIJK {
    let mut coord = CoordIJK{i: 0, j: 0, k: 0};
    let a1 = h2d.x.abs();
    let a2 = h2d.y.abs();

    let x2 = a2 / constants::M_SIN60;
    let x1 = a1 + x2 / 2.0;

    let m1 = x1 as i64;
    let m2 = x2 as i64;

    let r1 = x1 - (m1 as f64);
    let r2 = x2 - (m2 as f64);

    if r1 < 0.5 {
        if r1 < 1.0 / 3.0 {
            if r2 < (1.0 + r1) / 2.0 {
                coord.i = m1;
                coord.j = m2;
            } else {
                coord.i = m1;
                coord.j = m2 + 1;
            }
        } else {
            if r2 < (1.0 - r1) {
                coord.j = m2;
            } else {
                coord.j = m2 + 1;
            }

            if (1.0 - r1) <= r2 && r2 < (2.0 * r1) {
                coord.i = m1 + 1;
            } else {
                coord.i = m1;
            }
        }
    } else {
        if r1 < (2.0 / 3.0) {
            if r2 < (1.0 - r1) {
                coord.j = m2;
            } else {
                coord.j = m2 + 1;
            }

            if (2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1) {
                coord.i = m1;
            } else {
                coord.i = m1 + 1;
            }
        } else {
            if r2 < (r1 / 2.0) {
                coord.i = m1 + 1;
                coord.j = m2;
            } else {
                coord.i = m1 + 1;
                coord.j = m2 + 1;
            }
        }
    }

    if h2d.x < 0.0 {
        if coord.j % 2 == 0 {
            let axis_i = coord.j / 2;
            let diff = coord.i - axis_i;
            coord.i = coord.i - 2 * diff;
        } else {
            let axis_i = (coord.j + 1) / 2;
            let diff = coord.i - axis_i;
            coord.i = coord.i - (2 * diff + 1);
        }
    }

    if h2d.y < 0.0 {
        coord.i = coord.i - (2 * coord.j + 1) / 2;
        coord.j = coord.j * -1;
    }
    coord
}

fn normalize_ijk_coord(raw: CoordIJK) -> CoordIJK {
    let mut norm = CoordIJK{i: raw.i, j: raw.j, k: raw.k};
    if norm.i < 0 {
        norm.j -= norm.i;
        norm.k -= norm.i;
        norm.i = 0
    }

    if norm.j < 0 {
        norm.i -= norm.j;
        norm.k -= norm.j;
        norm.j = 0;
    }

    if norm.k < 0 {
        norm.i -= norm.k;
        norm.j -= norm.k;
        norm.k = 0;
    }

    let min = min(norm.i, min(norm.j, norm.k));
    if min > 0 {
        norm.i -= min;
        norm.j -= min;
        norm.k -= min;
    }
    norm
}

fn hex_2d_to_coord_ijk(h2d: &Vec2d) -> CoordIJK {
    normalize_ijk_coord(unnormalized_coord_ijk(h2d))
}


/**
 * Encodes a coordinate on the sphere to the FaceIJK address of the containing
 * cell at the specified resolution.
 *
 * geo: Lat/Lon coord in radians
 * res: Desired H3 Resolution
 */
fn geo_to_face_ijk(geo: &GeoCoord, res: H3Resolution) -> FaceIJK {
    let face_2d = geo_to_hex_2d(geo, res);
    let ijk = hex_2d_to_coord_ijk(&face_2d.coords);
    FaceIJK{face: face_2d.face, coord: ijk}
}

// H3 Index
// https://uber.github.io/h3/#/documentation/core-library/h3-index-representations
// 64 bits
// 0-3   - Mode setting -- 1000 for hexagon mode
// 4-6   - Reserved -- reserved / empty
// 7-10  - Resolution (4 bits, 0 - 15)
// 11-17 - Base Cell (0 - 121)
// 18-63 - Cell encoding - 45 bits
//         3 bits per level from 1 - 15

type H3Index = u64;
fn face_ijk_to_h3(fijk: &FaceIJK, res: H3Resolution) -> H3Index {
    let mut h3 = constants::H3_INIT;
    h3 = set_h3_mode(h3, constants::H3_HEXAGON_MODE);
    h3 = set_h3_resolution(h3, res);
    // * [X] Set Base Cell
    if res == 0 {
        let base_cell = face_ijk_to_base_cell(fijk);
        h3 = set_h3_base_cell(h3, base_cell);
        return h3;
    }
    h3
}

fn unit_ijk_to_direction(ijk: &CoordIJK) -> Direction {
    let ijk_norm = normalize_ijk_coord(ijk.clone());
    let mut dir = Direction::Invalid;
    for (i, unit_ijk) in constants::UNIT_VECS.iter().enumerate() {
        if ijk_norm == *unit_ijk {
            dir = Direction::from_int(i as i64);
            break;
        }
    }
    return dir
}

fn set_h3_index_digit(h3: H3Index, res: H3Resolution, dir: Direction) -> H3Index {
    h3
}

fn set_h3_index_digits(fijk: &FaceIJK, res: H3Resolution, h3: H3Index) -> H3Index {
    let mut fijk_bc = fijk.clone();
    let mut ijk = fijk.coord;
    for r in (0..res).rev() {
        println!("{}", r);
        let last_ijk = ijk;
        let last_center: CoordIJK;

        if is_class_iii_resolution(r + 1) {
            ijk = up_ap7(ijk);
            last_center = down_ap7(ijk);
        } else {
            ijk = up_ap7_rot(ijk);
            last_center = down_ap7_rot(ijk);
        }

        let diff = normalize_ijk_coord(last_ijk - last_center);
    }
    h3
}

fn set_h3_mode(index: H3Index, mode: u64) -> H3Index {
    (index & constants::H3_MODE_MASK_NEGATIVE) | (mode << constants::H3_MODE_OFFSET)
}

fn set_h3_resolution(index: H3Index, res: H3Resolution) -> H3Index {
    (index & constants::H3_RES_MASK_NEGATIVE) | ((res as u64) << constants::H3_RES_OFFSET)
}

fn set_h3_base_cell(index: H3Index, cell: BaseCell) -> H3Index {
    (index & constants::H3_BC_MASK_NEGATIVE) | ((cell as u64) << constants::H3_BC_OFFSET)
}

fn face_ijk_to_base_cell(fijk: &FaceIJK) -> BaseCell {
    constants::FACE_IJK_BASE_CELLS
        [fijk.face]
        [fijk.coord.i as usize]
        [fijk.coord.j as usize]
        [fijk.coord.k as usize]
        [0]
}

fn up_ap7(coord: CoordIJK) -> CoordIJK {
    let mut out = coord.clone();
    let i = coord.i - coord.k;
    let j = coord.j - coord.k;

    out.i = (((3 * i - j) as f64) / 7.0f64).round() as i64;
    out.j = (((i + 2 * j) as f64) / 7.0f64).round() as i64;
    out.k = 0;
    normalize_ijk_coord(out)
}

fn up_ap7_rot(coord: CoordIJK) -> CoordIJK {
    let mut out = coord.clone();
    let i = coord.i - coord.k;
    let j = coord.j - coord.k;

    out.i = (((2 * i + j) as f64) / 7.0f64).round() as i64;
    out.j = (((3 * j - i) as f64) / 7.0f64).round() as i64;
    out.k = 0;
    normalize_ijk_coord(out)
}

fn down_ap7(coord: CoordIJK) -> CoordIJK {
    let ivec = CoordIJK::new(3,0,1).scale(coord.i);
    let jvec = CoordIJK::new(1,3,0).scale(coord.j);
    let kvec = CoordIJK::new(0,1,3).scale(coord.k);

    let out = ivec + jvec + kvec;
    normalize_ijk_coord(out)
}

fn down_ap7_rot(coord: CoordIJK) -> CoordIJK {
    let ivec = CoordIJK::new(3,1,0).scale(coord.i);
    let jvec = CoordIJK::new(0,3,1).scale(coord.j);
    let kvec = CoordIJK::new(1,0,3).scale(coord.k);

    let out = ivec + jvec + kvec;
    normalize_ijk_coord(out)
}

#[cfg(test)]
mod tests {
    use *;
    use geo_to_coord_3d;
    use geo_to_hex_2d;
    use nearest_face_to_geo;
    use square_distance_3d;
    use Vec3d;
    use GeoCoord;
    use std::fs::File;
    use std::io::prelude::*;

    const DEFAULT_FLOAT_EPSILON: f64 = 0.0001;
    fn assert_eq_floats(a: f64, b: f64) -> () {
        assert_eq_floats_with_eps(a, b, DEFAULT_FLOAT_EPSILON);
    }
    fn assert_eq_floats_with_eps(a: f64, b: f64, epsilon: f64) -> () {
        assert!((a - b).abs() < epsilon,
                format!("Expected {} vs {} within {}", a, b, epsilon));
    }

    fn assert_eq_v3d(a: Vec3d, b: Vec3d) -> () {
        assert_eq_floats(a.x, b.x);
        assert_eq_floats(a.y, b.y);
        assert_eq_floats(a.z, b.z);
    }

    #[test]
    fn test_geo_to_coord_3d() {
        let pairs = vec![
            (GeoCoord::new(-0.094485, 1.552838),
             Vec3d::new(0.017878, 0.995379, -0.094345)),
            (GeoCoord::new(-0.056663, 1.634397),
             Vec3d::new(-0.063456, 0.996376, -0.056633)),
            (GeoCoord::new(-0.219365, 1.115781),
             Vec3d::new(0.428944, 0.876729, -0.217610))
        ];
        for (coord, v3d) in pairs {
            assert_eq_v3d(v3d, geo_to_coord_3d(&coord));
        }
    }

    fn test_tsv(name: &str) -> Vec<Vec<String>> {
        let mut f = File::open(format!("./tests/resources/{}", name)).expect("file not found");
        let mut contents = String::new();
        f.read_to_string(&mut contents).expect("Failed to read test file");

        let lines = contents.split("\n");

        let tsv = lines
            .filter(|l| !l.is_empty())
            .map(|l| l.split("\t").map(|s| s.to_owned()).collect());
        tsv.collect()
    }

    fn geo_to_face_test_cases() -> Vec<(GeoCoord, usize)> {
        test_tsv("geo_coord_face_test_cases.tsv").iter().map(|cols| {
            let lat_rads: f64 = cols[0].parse().unwrap();
            let lon_rads: f64 = cols[1].parse().unwrap();
            let face: usize = cols[2].parse().unwrap();
            (GeoCoord::new(lat_rads, lon_rads), face)
        }).collect()
    }

    #[test]
    fn test_geo_to_face() {
        for pair in geo_to_face_test_cases() {
            assert_eq!(
                nearest_face_to_geo(&pair.0).0,
                pair.1,
                "Expected GeoCoord {} to be nearest face {}",
                pair.0,
                pair.1
            );
        }
    }

    #[test]
    fn test_face_2d_coords_for_face_centers() {
        test_tsv("face_2d_zero_coord_examples.tsv").iter().for_each(|cols| {
            let lat_rads: f64 = cols[0].parse().unwrap();
            let lon_rads: f64 = cols[1].parse().unwrap();
            let face: usize = cols[2].parse().unwrap();
            let geo = GeoCoord::new(lat_rads, lon_rads);

            let exp_face_coord = FaceCoord2d{face: face, coords: Vec2d::zero()};
            assert_eq!(
                geo_to_hex_2d(&geo, 0),
                exp_face_coord,
                "Expected Geo Coord {} to map to FaceCoord2d {}",
                geo,
                exp_face_coord
            )
        });
    }

    #[test]
    fn test_square_distance_3d() {
        let v1 = Vec3d::new(0.0, 0.0, 0.0);
        let v2 = Vec3d::new(1.0, 0.0, 0.0);
        let v3 = Vec3d::new(0.0, 1.0, 1.0);
        let v4 = Vec3d::new(1.0, 1.0, 1.0);
        let v5 = Vec3d::new(1.0, 1.0, 2.0);

        assert_eq_floats(square_distance_3d(&v1, &v1), 0.0);
        assert_eq_floats(square_distance_3d(&v1, &v2), 1.0);
        assert_eq_floats(square_distance_3d(&v1, &v3), 2.0);
        assert_eq_floats(square_distance_3d(&v1, &v4), 3.0);
        assert_eq_floats(square_distance_3d(&v1, &v5), 6.0);
    }

    #[test]
    fn test_geo_azimuth_radians() {
        let cases = vec![
            (GeoCoord::new(1.054751, -1.347517) , GeoCoord::new(0.757873, -1.700381) , -2.361377),
            (GeoCoord::new(-0.491715, -2.739604), GeoCoord::new(-0.679978, -2.580232), 2.569891),
            (GeoCoord::new(-0.491715, -2.739604), GeoCoord::new(-0.211030, -2.975053), -0.712065),
            (GeoCoord::new(0.605929, 2.953923)  , GeoCoord::new(0.919620, 2.719804)  , -0.416341),
            (GeoCoord::new(0.079066, 2.408163)  , GeoCoord::new(0.101005, 2.400045)  , -0.352712),
            (GeoCoord::new(-0.605929, -0.187669), GeoCoord::new(-0.189070, -0.236293), -0.117554),
            (GeoCoord::new(-0.172745, -1.463446), GeoCoord::new(-0.091844, -1.669910), -1.209337),
            (GeoCoord::new(-0.172745, -1.463446), GeoCoord::new(-0.039764, -1.769705), -1.178624),
            (GeoCoord::new(-0.803583, -1.893195), GeoCoord::new(-0.567315, -1.768839), 0.427850),
            (GeoCoord::new(-1.307748, -0.604648), GeoCoord::new(-1.186825, -0.604890), -0.000753),
            (GeoCoord::new(-0.600192, 2.690989) , GeoCoord::new(-0.704096, 2.934774) , 2.135086),
            (GeoCoord::new(1.307748, 2.536945)  , GeoCoord::new(1.119119, -2.838903) , 1.645120),
            (GeoCoord::new(0.803583, 1.248397)  , GeoCoord::new(0.921756, 0.832705)  , -1.005149),
            (GeoCoord::new(-1.307748, -0.604648), GeoCoord::new(-1.011638, 0.063376) , 1.065688),
            (GeoCoord::new(-0.491715, -2.739604), GeoCoord::new(-0.344152, -2.672664), 0.407104),
            (GeoCoord::new(-1.054751, 1.794075) , GeoCoord::new(-1.223819, 2.037650) , 2.707681),
            (GeoCoord::new(0.427371, -1.888876) , GeoCoord::new(0.606291, -1.851838) , 0.169130),
            (GeoCoord::new(1.054751, -1.347517) , GeoCoord::new(0.757870, -1.700376) , -2.361388),
            (GeoCoord::new(0.427371, -1.888876) , GeoCoord::new(0.019223, -1.958678) , -2.967246),
            (GeoCoord::new(-1.054751, 1.794075) , GeoCoord::new(-1.158853, 1.237223) , -2.207662)
        ];

        for (a, b, azimuth) in cases {
            assert_eq_floats(azimuth, geo_azimuth_radians(&a, &b));
        }
    }

    #[test]
    fn test_positive_angle_rads() {
        let cases = vec![
            (-1.387236, 4.895949),
            (-0.895065, 5.388120),
            (-0.642279, 5.640906),
            (-1.278988, 5.004198),
            (-1.142544, 5.140641),
            (-2.778123, 3.505063),
            (-0.648280, 5.634906),
            (-2.500741, 3.782444),
            (9.004937, 2.7217516),
            (-1.314182, 4.969003),
            (-1.097476, 5.185710),
            (-2.181562, 4.101623),
            (7.209724, 0.9265386),
            (-2.567716, 3.715469),
            (-0.522846, 5.760339),
            (-0.748559, 5.534626),
            (-0.231651, 6.051535),
            (-2.598210, 3.684975),
            (-0.999930, 5.283255),
            (-1.474177, 4.809009)
        ];

        for (a,b) in cases {
            assert_eq_floats(b, positive_angle_radians(a));
        }
    }

    #[test]
    fn test_face_theta() {
        let cases = vec![
            (GeoCoord::new(1.258573, 2.039761), 1, 1.138107),
            (GeoCoord::new(-0.677089, -0.758010), 13, 2.053503),
            (GeoCoord::new(-0.938231, 3.000496), 15, 0.047635),
            (GeoCoord::new(1.157020, -2.137218), 2, 1.726973),
            (GeoCoord::new(0.375171, -1.478357), 7, 1.907947),
            (GeoCoord::new(0.666978, 2.478336), 6, 4.255309),
            (GeoCoord::new(0.499680, -0.299140), 3, 4.522263),
            (GeoCoord::new(0.770502, -1.553657), 2, 3.419251),
            (GeoCoord::new(0.038155, 1.586939), 5, 5.235988),
            (GeoCoord::new(-0.439301, 0.681398), 9, 0.505100),
            (GeoCoord::new(1.009809, 0.577651), 0, 0.191516),
            (GeoCoord::new(0.713048, -0.280992), 3, 5.892041),
            (GeoCoord::new(-0.609226, -0.388859), 13, 1.806873),
            (GeoCoord::new(-0.677665, 2.392681), 15, 4.682011),
            (GeoCoord::new(-0.108839, 2.286307), 10, 2.213284),
            (GeoCoord::new(0.204382, -0.712922), 8, 3.422628),
            (GeoCoord::new(-0.629125, 0.181376), 13, 4.687292),
            (GeoCoord::new(0.340684, -2.290744), 11, 5.217366),
            (GeoCoord::new(1.304647, -2.851584), 1, 4.610854),
            (GeoCoord::new(-0.161922, -2.342135), 11, 3.932146),
            (GeoCoord::new(-0.747909, 3.133747), 15, 0.603066),
            (GeoCoord::new(-0.463428, -1.088851), 12, 4.423758),
            (GeoCoord::new(-0.486034, 2.020983), 15, 4.271083),
            (GeoCoord::new(0.034264, -2.482857), 11, 3.946226),
            (GeoCoord::new(1.251998, 5.660560), 2, 0.114276),
            (GeoCoord::new(-0.807800, -2.094861), 17, 5.478552),
            (GeoCoord::new(0.218812, 4.900052), 12, 0.246037),
            (GeoCoord::new(-0.668244, 2.984760), 15, 0.774780),
            (GeoCoord::new(0.722083, -0.884139), 3, 1.529580),
            (GeoCoord::new(-0.830868, 1.500760), 19, 3.141593)
        ];

        for (geo, face, exp_theta) in cases {
            let theta = face_theta(&geo, face);
            assert_eq_floats(theta, exp_theta);
        }
    }

    #[test]
    fn test_checking_resolution_class() {
        assert!(is_class_iii_resolution(1));
        assert!(is_class_iii_resolution(5));
        assert!(!is_class_iii_resolution(2));
        assert!(!is_class_iii_resolution(4));
        assert!(is_class_ii_resolution(0));
        assert!(is_class_ii_resolution(2));
        assert!(!is_class_ii_resolution(7));
        assert!(!is_class_ii_resolution(9));
    }

    fn assert_eq_face_coords(a: FaceCoord2d, b: FaceCoord2d) -> () {
        assert_eq!(a.face, b.face, "Expected face {} to match {}", a.face, b.face);
        assert_eq_floats_with_eps(a.coords.x, b.coords.x, 0.00001);
        assert_eq_floats_with_eps(a.coords.y, b.coords.y, 0.00001);
    }

    #[test]
    fn test_geo_to_hex2d_results() {
        for cols in test_tsv("hex_2d_non_zero_cases.tsv") {
            let lat_rads: f64 = cols[0].parse().unwrap();
            let lon_rads: f64 = cols[1].parse().unwrap();
            let res: H3Resolution = cols[2].parse().unwrap();
            let face: usize = cols[3].parse().unwrap();
            let x: f64 = cols[4].parse().unwrap();
            let y: f64 = cols[5].parse().unwrap();

            let geo = GeoCoord::new(lat_rads, lon_rads);

            assert_eq_face_coords(
                geo_to_hex_2d(&geo, res),
                FaceCoord2d{
                    face: face,
                    coords: Vec2d::new(x, y)
                }
            );

        }
    }

    #[test]
    fn test_hex2d_unnorm_ijk_conversion() {
        for cols in test_tsv("unnorm_ijk_cases.tsv") {
            let h2d_x: f64 = cols[0].parse().unwrap();
            let h2d_y: f64 = cols[1].parse().unwrap();
            let exp_i: i64 = cols[2].parse().unwrap();
            let exp_j: i64 = cols[3].parse().unwrap();
            let exp_k: i64 = cols[4].parse().unwrap();

            let h2d = Vec2d::new(h2d_x, h2d_y);
            let exp_ijk = CoordIJK::new(exp_i, exp_j, exp_k);
            let ijk = unnormalized_coord_ijk(&h2d);

            assert_eq!(exp_ijk, ijk);
        }
    }

    #[test]
    fn test_normalize_ijk_coords() {
        for cols in test_tsv("ijk_normalize_cases.tsv") {
            let raw_i: i64 = cols[0].parse().unwrap();
            let raw_j: i64 = cols[1].parse().unwrap();
            let raw_k: i64 = cols[2].parse().unwrap();
            let norm_i: i64 = cols[3].parse().unwrap();
            let norm_j: i64 = cols[4].parse().unwrap();
            let norm_k: i64 = cols[5].parse().unwrap();

            let raw_ijk = CoordIJK::new(raw_i, raw_j, raw_k);
            let exp_ijk = CoordIJK::new(norm_i, norm_j, norm_k);

            assert_eq!(exp_ijk, normalize_ijk_coord(raw_ijk));
        }
    }

    #[test]
    fn test_hex2d_to_coord_ijk() {
        for cols in test_tsv("hex2d_coordijk_cases.tsv") {
            let hex_x: f64 = cols[0].parse().unwrap();
            let hex_y: f64 = cols[1].parse().unwrap();
            let coord_i: i64 = cols[2].parse().unwrap();
            let coord_j: i64 = cols[3].parse().unwrap();
            let coord_k: i64 = cols[4].parse().unwrap();

            let h2d = Vec2d::new(hex_x, hex_y);
            let ijk = CoordIJK::new(coord_i, coord_j, coord_k);

            assert_eq!(ijk, hex_2d_to_coord_ijk(&h2d));
        }
    }

    #[test]
    fn test_geo_to_face_ijk() {
        // TODO -- why is this one bad?
        // 1.125543981167061113879412914685	-3.089471475623864371584659238579	4	3	92	6	0
        for cols in test_tsv("geo_to_face_ijk_cases.tsv") {
            let geo_lat: f64 = cols[0].parse().unwrap();
            let geo_lon: f64 = cols[1].parse().unwrap();
            let res: H3Resolution = cols[2].parse().unwrap();
            let face: usize = cols[3].parse().unwrap();
            let i: i64 = cols[4].parse().unwrap();
            let j: i64 = cols[5].parse().unwrap();
            let k: i64 = cols[6].parse().unwrap();

            let geo = GeoCoord::new(geo_lat, geo_lon);
            let face_ijk = FaceIJK::new(face, i, j, k);

            assert_eq!(face_ijk, geo_to_face_ijk(&geo, res));
        }
    }

    #[test]
    fn test_h3_set_mode() {
        assert_eq!(576495936675512319, set_h3_mode(constants::H3_INIT, constants::H3_HEXAGON_MODE));
    }

    #[test]
    fn test_h3_set_resolution() {
        assert_eq!(612524733694476287, set_h3_resolution(576495936675512319, 8));
        assert_eq!(617028333321846783, set_h3_resolution(576495936675512319, 9));
        assert_eq!(621531932949217279, set_h3_resolution(576495936675512319, 10));
        assert_eq!(630539132203958271, set_h3_resolution(576495936675512319, 12));
        assert_eq!(635042731831328767, set_h3_resolution(576495936675512319, 13));
        assert_eq!(639546331458699263, set_h3_resolution(576495936675512319, 14));
    }

    #[test]
    #[ignore]
    fn test_face_ijk_to_h3_invalid_index() {
    }

    #[test]
    fn test_h3_set_base_cell() {
        let start = 576495936675512319;
        assert_eq!(576495936675512319, set_h3_base_cell(start, 0));
        assert_eq!(576566305419689983, set_h3_base_cell(start, 2));
        assert_eq!(576812596024311807, set_h3_base_cell(start, 9));
        assert_eq!(577058886628933631, set_h3_base_cell(start, 16));
        assert_eq!(577269992861466623, set_h3_base_cell(start, 22));
    }

    #[test]
    fn test_face_ijk_to_base_cell() {
        for cols in test_tsv("face_ijk_to_base_cell_cases.tsv") {
            let face: usize = cols[0].parse().unwrap();
            let i: i64 = cols[1].parse().unwrap();
            let j: i64 = cols[2].parse().unwrap();
            let k: i64 = cols[3].parse().unwrap();
            let base_cell: u8 = cols[4].parse().unwrap();

            let face_ijk = FaceIJK::new(face, i, j, k);
            assert_eq!(base_cell, face_ijk_to_base_cell(&face_ijk));
        }
    }

    #[test]
    fn test_face_ijk_to_h3_res_0() {
        for cols in test_tsv("face_ijk_to_h3_res_0_cases.tsv") {
            let face: usize = cols[0].parse().unwrap();
            let i: i64 = cols[1].parse().unwrap();
            let j: i64 = cols[2].parse().unwrap();
            let k: i64 = cols[3].parse().unwrap();
            let h3: H3Index = cols[4].parse().unwrap();

            let face_ijk = FaceIJK::new(face, i, j, k);
            assert_eq!(h3, face_ijk_to_h3(&face_ijk, 0));
        }
    }

    #[test]
    fn test_up_ap7() {
        for cols in test_tsv("up_ap7_cases.tsv") {
            let i0: i64 = cols[0].parse().unwrap();
            let j0: i64 = cols[1].parse().unwrap();
            let k0: i64 = cols[2].parse().unwrap();
            let i1: i64 = cols[3].parse().unwrap();
            let j1: i64 = cols[4].parse().unwrap();
            let k1: i64 = cols[5].parse().unwrap();

            let c0 = CoordIJK::new(i0, j0, k0);
            let c1 = CoordIJK::new(i1, j1, k1);
            assert_eq!(c1, up_ap7(c0));
        }
    }

    #[test]
    fn test_down_ap7() {
        for cols in test_tsv("down_ap7_cases.tsv") {
            let i0: i64 = cols[0].parse().unwrap();
            let j0: i64 = cols[1].parse().unwrap();
            let k0: i64 = cols[2].parse().unwrap();
            let i1: i64 = cols[3].parse().unwrap();
            let j1: i64 = cols[4].parse().unwrap();
            let k1: i64 = cols[5].parse().unwrap();

            let c0 = CoordIJK::new(i0, j0, k0);
            let c1 = CoordIJK::new(i1, j1, k1);
            assert_eq!(c1, down_ap7(c0));
        }
    }

    #[test]
    fn test_up_ap7_rot() {
        for cols in test_tsv("up_ap7_rot_cases.tsv") {
            let i0: i64 = cols[0].parse().unwrap();
            let j0: i64 = cols[1].parse().unwrap();
            let k0: i64 = cols[2].parse().unwrap();
            let i1: i64 = cols[3].parse().unwrap();
            let j1: i64 = cols[4].parse().unwrap();
            let k1: i64 = cols[5].parse().unwrap();

            let c0 = CoordIJK::new(i0, j0, k0);
            let c1 = CoordIJK::new(i1, j1, k1);
            assert_eq!(c1, up_ap7_rot(c0));
        }
    }

    #[test]
    fn test_down_ap7_rot() {
        for cols in test_tsv("down_ap7_rot_cases.tsv") {
            let i0: i64 = cols[0].parse().unwrap();
            let j0: i64 = cols[1].parse().unwrap();
            let k0: i64 = cols[2].parse().unwrap();
            let i1: i64 = cols[3].parse().unwrap();
            let j1: i64 = cols[4].parse().unwrap();
            let k1: i64 = cols[5].parse().unwrap();

            let c0 = CoordIJK::new(i0, j0, k0);
            let c1 = CoordIJK::new(i1, j1, k1);
            assert_eq!(c1, down_ap7_rot(c0));
        }
    }

    #[test]
    #[ignore]
    fn test_set_h3_index_digits() {
        for cols in test_tsv("set_h3_index_digits_cases.tsv") {
            let face: usize = cols[0].parse().unwrap();
            let i: i64 = cols[1].parse().unwrap();
            let j: i64 = cols[2].parse().unwrap();
            let k: i64 = cols[3].parse().unwrap();
            let res: H3Resolution = cols[4].parse().unwrap();
            let h3_input: H3Index = cols[5].parse().unwrap();
            let h3_output: H3Index = cols[6].parse().unwrap();

            let fijk = FaceIJK::new(face, i, j, k);
            assert_eq!(h3_output, set_h3_index_digits(&fijk, res, h3_input));
        }
    }

    #[test]
    fn test_unit_ijk_to_digit() {
        for cols in test_tsv("unit_ijk_to_direction_cases.tsv") {
            let i: i64 = cols[0].parse().unwrap();
            let j: i64 = cols[1].parse().unwrap();
            let k: i64 = cols[2].parse().unwrap();
            let dir_i: i64 = cols[3].parse().unwrap();

            let ijk = CoordIJK::new(i,j,k);
            let dir = Direction::from_int(dir_i);

            assert_eq!(dir, unit_ijk_to_direction(&ijk));
        }
    }
}
