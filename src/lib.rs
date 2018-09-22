#![allow(dead_code)]
extern crate geo_types;
use std::fmt;

// use geo_types::Coordinate;

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
struct GeoCoord {
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

fn geo_to_face_2d(geo: &GeoCoord) -> FaceCoord2d {
    let (face, dist) = nearest_face_to_geo(geo);

    let r = (1.0 - dist / 2.0).acos();
    if r < constants::FACE_DISTANCE_EPSILON {
        FaceCoord2d{face: face, coords: Vec2d::zero()}
    } else {
        FaceCoord2d{face: face, coords: Vec2d::zero()}
    }
}

#[cfg(test)]
mod tests {
    use *;
    use geo_to_coord_3d;
    use geo_to_face_2d;
    use nearest_face_to_geo;
    use square_distance_3d;
    use Vec3d;
    use GeoCoord;
    use std::fs::File;
    use std::io::prelude::*;

    fn assert_eq_floats(a: f64, b: f64) -> () {
        let epsilon = 0.000001;
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
                geo_to_face_2d(&geo),
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
}
