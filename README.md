# h3-rs

This is a very **WIP** effort to port Uber's [H3](https://github.com/uber/h3) spatial grid library from C to Rust.

### Initial Goal: Lat/Lon -> H3 ID Conversion

* [X] geoToFace(coord) -> face (numeric id)
* [x] geoToFaceIJK(coord, res) -> faceIJK coord
* [X] geoToHex2d(coord, res) -> (face, vec2D)
* [x] geoToVec3d(coord) -> 3dCoord
* [X] hex2dToCoordIJK(vec2D) -> CoordIJK
* [X] pointSquareDist(vec3d, vec3d)
* [ ] geoToH3(coord, res) -> h3 id
