extern crate geo_types;

use geo_types::Coordinate;

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

impl GeoCoord {
    fn new(lat: f64, lon: f64) -> GeoCoord {
        GeoCoord{lat: lat, lon: lon}
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

fn geo_to_face(geo: GeoCoord) -> usize {
    let v3d = geo_to_coord_3d(geo);
    let mut face = 0;
    let mut min_dist = square_distance_3d(&v3d, &constants::FACE_CENTERS[0]);
    for i in (1..constants::NUM_ICOSA_FACES) {
        let dist = square_distance_3d(&v3d, &constants::FACE_CENTERS[i]);
        if (dist < min_dist) {
            face = i;
            min_dist = dist;
        }
    }
    face
}

fn geo_to_coord_3d(geo: GeoCoord) -> Vec3d {
    let r = geo.lat.cos();

    Vec3d::new(
        geo.lon.cos() * r,
        geo.lon.sin() * r,
        geo.lat.sin()
    )
}
// 27: ~~~ Get 3d (x,y,z) coords for lat/lon in rads
// 27: ~~~ (x,y,z) for rads (-0.094485, 1.552838)
// 27: ~~~ (0.017878, 0.995379, -0.094345)
// 27: ~~~ Coord: (-0.056663, 1.634397)~~~ 3 - _geoToVec3d
// 27: ~~~ Get 3d (x,y,z) coords for lat/lon in rads
// 27: ~~~ (x,y,z) for rads (-0.056663, 1.634397)
// 27: ~~~ (-0.063456, 0.996376, -0.056633)
// 27: ~~~ Coord: (-0.219365, 1.115781)~~~ 3 - _geoToVec3d
// 27: ~~~ Get 3d (x,y,z) coords for lat/lon in rads
// 27: ~~~ (x,y,z) for rads (-0.219365, 1.115781)
// 27: ~~~ (0.428944, 0.876729, -0.217610)

#[cfg(test)]
mod tests {
    use geo_types::Coordinate;
    use geo_to_coord_3d;
    use Vec3d;
    use GeoCoord;

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
            assert_eq_v3d(v3d, geo_to_coord_3d(coord));
        }
    }

    #[test]
    fn test_geo_to_face() {
    }

    #[test]
    fn test_point_square_distance() {
    // TEST(_pointSquareDist) {
    //     Vec3d v1 = {0, 0, 0};
    //     Vec3d v2 = {1, 0, 0};
    //     Vec3d v3 = {0, 1, 1};
    //     Vec3d v4 = {1, 1, 1};
    //     Vec3d v5 = {1, 1, 2};

    //     t_assert(fabs(_pointSquareDist(&v1, &v1)) < DBL_EPSILON,
    //              "distance to self is 0");
    //     t_assert(fabs(_pointSquareDist(&v1, &v2) - 1) < DBL_EPSILON,
    //              "distance to <1,0,0> is 1");
    //     t_assert(fabs(_pointSquareDist(&v1, &v3) - 2) < DBL_EPSILON,
    //              "distance to <0,1,1> is 2");
    //     t_assert(fabs(_pointSquareDist(&v1, &v4) - 3) < DBL_EPSILON,
    //              "distance to <1,1,1> is 3");
    //     t_assert(fabs(_pointSquareDist(&v1, &v5) - 6) < DBL_EPSILON,
    //              "distance to <1,1,2> is 6");
    // }
    }
}
