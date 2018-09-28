#![allow(dead_code)]
extern crate geo_types;
use std::fmt;
use std::f64::consts::PI;

// use geo_types::Coordinate;

const TWO_PI: f64 = PI * 2.0;
mod constants;

// [ ] geoToH3(coord, res) -> h3 id
// [ ] geoToFace(coord) -> face (numeric id)
// [ ] geoToFaceIJK(coord, res) -> faceIJK coord
// [ ] geoToHex2d(coord, res) -> (face, vec2D)
// [x] geoToVec3d(coord) -> 3dCoord
// [ ] hex2dToCoordIJK(vec2D, )
// [ ] pointSquareDist(vec3d, vec3d)
// Most coord ops in radians

// FaceIJK:
// typedef struct {
//     int face;        ///< face number
//     CoordIJK coord;  ///< ijk coordinates on that face
// } FaceIJK;

// typedef struct {
//     int i;  ///< i component
//     int j;  ///< j component
//     int k;  ///< k component
// } CoordIJK;

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

fn geo_to_hex_2d(geo: &GeoCoord) -> FaceCoord2d {
    let (face, dist) = nearest_face_to_geo(geo);

    let r = (1.0 - dist / 2.0).acos();
    if r < constants::FACE_DISTANCE_EPSILON {
        return FaceCoord2d{face: face, coords: Vec2d::zero()};
    }

    let face_center_geo = &constants::FACE_CENTER_GEO_COORDS[face];
    FaceCoord2d{face: face, coords: Vec2d::zero()}
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

    fn assert_eq_floats(a: f64, b: f64) -> () {
        let epsilon = 0.0001;
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
                geo_to_hex_2d(&geo),
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
}
