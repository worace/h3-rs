use types::CoordIJK;
use types::H3Resolution;
use types::BaseCellData;
use types::FaceIJK;
use Vec3d;
use GeoCoord;

pub const NUM_ICOSA_FACES: usize = 20;

pub static FACE_CENTERS: [Vec3d; NUM_ICOSA_FACES] = [
   Vec3d{x: 0.2199307791404606,  y: 0.6583691780274996,  z: 0.7198475378926182},  // face  0
   Vec3d{x: -0.2139234834501421, y: 0.1478171829550703,  z: 0.9656017935214205},  // face  1
   Vec3d{x: 0.1092625278784797,  y: -0.4811951572873210, z: 0.8697775121287253},  // face  2
   Vec3d{x: 0.7428567301586791,  y: -0.3593941678278028, z: 0.5648005936517033},  // face  3
   Vec3d{x: 0.8112534709140969,  y: 0.3448953237639384,  z: 0.4721387736413930},  // face  4
   Vec3d{x: -0.1055498149613921, y: 0.9794457296411413,  z: 0.1718874610009365},  // face  5
   Vec3d{x: -0.8075407579970092, y: 0.1533552485898818,  z: 0.5695261994882688},  // face  6
   Vec3d{x: -0.2846148069787907, y: -0.8644080972654206, z: 0.4144792552473539},  // face  7
   Vec3d{x: 0.7405621473854482,  y: -0.6673299564565524, z: -0.0789837646326737}, // face  8
   Vec3d{x: 0.8512303986474293,  y: 0.4722343788582681,  z: -0.2289137388687808}, // face  9
   Vec3d{x: -0.7405621473854481, y: 0.6673299564565524,  z: 0.0789837646326737},  // face 10
   Vec3d{x: -0.8512303986474292, y: -0.4722343788582682, z: 0.2289137388687808},  // face 11
   Vec3d{x: 0.1055498149613919,  y: -0.9794457296411413, z: -0.1718874610009365}, // face 12
   Vec3d{x: 0.8075407579970092,  y: -0.1533552485898819, z: -0.5695261994882688}, // face 13
   Vec3d{x: 0.2846148069787908,  y: 0.8644080972654204,  z: -0.4144792552473539}, // face 14
   Vec3d{x: -0.7428567301586791, y: 0.3593941678278027,  z: -0.5648005936517033}, // face 15
   Vec3d{x: -0.8112534709140971, y: -0.3448953237639382, z: -0.4721387736413930}, // face 16
   Vec3d{x: -0.2199307791404607, y: -0.6583691780274996, z: -0.7198475378926182}, // face 17
   Vec3d{x: 0.2139234834501420,  y: -0.1478171829550704, z: -0.9656017935214205}, // face 18
   Vec3d{x: -0.1092625278784796, y: 0.4811951572873210,  z: -0.8697775121287253}, // face 19
];

pub const FACE_DISTANCE_EPSILON: f64 = 0.0000000000000001;

pub static FACE_CENTER_GEO_COORDS: [GeoCoord; NUM_ICOSA_FACES] = [
    GeoCoord{lat: 0.803582649718989942, lon: 1.248397419617396099},    // face  0
    GeoCoord{lat: 1.307747883455638156, lon: 2.536945009877921159},    // face  1
    GeoCoord{lat: 1.054751253523952054, lon: -1.347517358900396623},   // face  2
    GeoCoord{lat: 0.600191595538186799, lon: -0.450603909469755746},   // face  3
    GeoCoord{lat: 0.491715428198773866, lon: 0.401988202911306943},    // face  4
    GeoCoord{lat: 0.172745327415618701, lon: 1.678146885280433686},    // face  5
    GeoCoord{lat: 0.605929321571350690, lon: 2.953923329812411617},    // face  6
    GeoCoord{lat: 0.427370518328979641, lon: -1.888876200336285401},   // face  7
    GeoCoord{lat: -0.079066118549212831, lon: -0.733429513380867741},  // face  8
    GeoCoord{lat: -0.230961644455383637, lon: 0.506495587332349035},   // face  9
    GeoCoord{lat: 0.079066118549212831, lon: 2.408163140208925497},    // face 10
    GeoCoord{lat: 0.230961644455383637, lon: -2.635097066257444203},   // face 11
    GeoCoord{lat: -0.172745327415618701, lon: -1.463445768309359553},  // face 12
    GeoCoord{lat: -0.605929321571350690, lon: -0.187669323777381622},  // face 13
    GeoCoord{lat: -0.427370518328979641, lon: 1.252716453253507838},   // face 14
    GeoCoord{lat: -0.600191595538186799, lon: 2.690988744120037492},   // face 15
    GeoCoord{lat: -0.491715428198773866, lon: -2.739604450678486295},  // face 16
    GeoCoord{lat: -0.803582649718989942, lon: -1.893195233972397139},  // face 17
    GeoCoord{lat: -1.307747883455638156, lon: -0.604647643711872080},  // face 18
    GeoCoord{lat: -1.054751253523952054, lon: 1.794075294689396615},   // face 19
];

pub static CLASS_II_FACE_IJK_AXES: [[f64; 3]; NUM_ICOSA_FACES] = [
    [5.619958268523939882, 3.525563166130744542, 1.431168063737548730],  // face  0
    [5.760339081714187279, 3.665943979320991689, 1.571548876927796127],  // face  1
    [0.780213654393430055, 4.969003859179821079, 2.874608756786625655],  // face  2
    [0.430469363979999913, 4.619259568766391033, 2.524864466373195467],  // face  3
    [6.130269123335111400, 4.035874020941915804, 1.941478918548720291],  // face  4
    [2.692877706530642877, 0.598482604137447119, 4.787272808923838195],  // face  5
    [2.982963003477243874, 0.888567901084048369, 5.077358105870439581],  // face  6
    [3.532912002790141181, 1.438516900396945656, 5.627307105183336758],  // face  7
    [3.494305004259568154, 1.399909901866372864, 5.588700106652763840],  // face  8
    [3.003214169499538391, 0.908819067106342928, 5.097609271892733906],  // face  9
    [5.930472956509811562, 3.836077854116615875, 1.741682751723420374],  // face 10
    [0.138378484090254847, 4.327168688876645809, 2.232773586483450311],  // face 11
    [0.448714947059150361, 4.637505151845541521, 2.543110049452346120],  // face 12
    [0.158629650112549365, 4.347419854898940135, 2.253024752505744869],  // face 13
    [5.891865957979238535, 3.797470855586042958, 1.703075753192847583],  // face 14
    [2.711123289609793325, 0.616728187216597771, 4.805518392002988683],  // face 15
    [3.294508837434268316, 1.200113735041072948, 5.388903939827463911],  // face 16
    [3.804819692245439833, 1.710424589852244509, 5.899214794638635174],  // face 17
    [3.664438879055192436, 1.570043776661997111, 5.758833981448388027],  // face 18
    [2.361378999196363184, 0.266983896803167583, 4.455774101589558636],  // face 19
];

pub const AP7_ROT_RADS: f64 = 0.333473172251832115336090;
pub const RES0_U_GNOMONIC: f64 = 0.38196601125010500003;
pub const SQRT7: f64 = 2.6457513110645905905016157536392604257102;
pub const M_SQRT3_2: f64 = 0.8660254037844386467637231707529361834714;
pub const M_SIN60: f64 = M_SQRT3_2;
// H3 index with mode 0, res 0, base cell 0, and 7 for all index digits.
pub const H3_INIT: u64 = 35184372088831;
pub const H3_HEXAGON_MODE: u64 = 1;
pub const H3_MODE_OFFSET: u64 = 59;
// 1's in the 4 mode bits, 0's everywhere else.
pub const H3_MODE_MASK: u64 = 15 << H3_MODE_OFFSET;
pub const H3_MODE_MASK_NEGATIVE: u64 = !H3_MODE_MASK;
pub const H3_RES_OFFSET: u64 = 52;
pub const H3_RES_MASK: u64 = 15 << H3_RES_OFFSET;
pub const H3_RES_MASK_NEGATIVE: u64 = !H3_RES_MASK;
pub const H3_BC_OFFSET: u64 = 45;
pub const H3_BC_MASK: u64 = 127 << H3_BC_OFFSET;
pub const H3_BC_MASK_NEGATIVE: u64 = !H3_BC_MASK;

// static const double faceAxesAzRadsCII[NUM_ICOSA_FACES][3] = {
//     {5.619958268523939882, 3.525563166130744542,
//      1.431168063737548730},  // face  0
//     {5.760339081714187279, 3.665943979320991689,
//      1.571548876927796127},  // face  1
//     {0.780213654393430055, 4.969003859179821079,
//      2.874608756786625655},  // face  2
//     {0.430469363979999913, 4.619259568766391033,
//      2.524864466373195467},  // face  3
//     {6.130269123335111400, 4.035874020941915804,
//      1.941478918548720291},  // face  4
//     {2.692877706530642877, 0.598482604137447119,
//      4.787272808923838195},  // face  5
//     {2.982963003477243874, 0.888567901084048369,
//      5.077358105870439581},  // face  6
//     {3.532912002790141181, 1.438516900396945656,
//      5.627307105183336758},  // face  7
//     {3.494305004259568154, 1.399909901866372864,
//      5.588700106652763840},  // face  8
//     {3.003214169499538391, 0.908819067106342928,
//      5.097609271892733906},  // face  9
//     {5.930472956509811562, 3.836077854116615875,
//      1.741682751723420374},  // face 10
//     {0.138378484090254847, 4.327168688876645809,
//      2.232773586483450311},  // face 11
//     {0.448714947059150361, 4.637505151845541521,
//      2.543110049452346120},  // face 12
//     {0.158629650112549365, 4.347419854898940135,
//      2.253024752505744869},  // face 13
//     {5.891865957979238535, 3.797470855586042958,
//      1.703075753192847583},  // face 14
//     {2.711123289609793325, 0.616728187216597771,
//      4.805518392002988683},  // face 15
//     {3.294508837434268316, 1.200113735041072948,
//      5.388903939827463911},  // face 16
//     {3.804819692245439833, 1.710424589852244509,
//      5.899214794638635174},  // face 17
//     {3.664438879055192436, 1.570043776661997111,
//      5.758833981448388027},  // face 18
//     {2.361378999196363184, 0.266983896803167583,
//      4.455774101589558636},  // face 19
// };

/** Resolution 0 base cell lookup table for each face.
 *
 * Given the face number and a resolution 0 ijk+ coordinate in that face's
 * face-centered ijk coordinate system, gives the base cell located at that
 * coordinate and the number of 60 ccw rotations to rotate into that base
 * cell's orientation.
 *
 * Valid lookup coordinates are from (0, 0, 0) to (2, 2, 2).
 */
pub static FACE_IJK_BASE_CELLS: [[[[[u8; 2]; 3]; 3]; 3]; NUM_ICOSA_FACES] = [
    [// face 0
        [
            // i 0
            [[16, 0], [18, 0], [24, 0]],  // j 0
            [[33, 0], [30, 0], [32, 3]],  // j 1
            [[49, 1], [48, 3], [50, 3]]   // j 2
        ],
        [
            // i 1
            [[8, 0], [5, 5], [10, 5]],    // j 0
            [[22, 0], [16, 0], [18, 0]],  // j 1
            [[41, 1], [33, 0], [30, 0]]   // j 2
        ],
        [
            // i 2
            [[4, 0], [0, 5], [2, 5]],    // j 0
            [[15, 1], [8, 0], [5, 5]],   // j 1
            [[31, 1], [22, 0], [16, 0]]  // j 2
        ]],
    [// face 1
        [
            // i 0
            [[2, 0], [6, 0], [14, 0]],    // j 0
            [[10, 0], [11, 0], [17, 3]],  // j 1
            [[24, 1], [23, 3], [25, 3]]   // j 2
        ],
        [
            // i 1
            [[0, 0], [1, 5], [9, 5]],    // j 0
            [[5, 0], [2, 0], [6, 0]],    // j 1
            [[18, 1], [10, 0], [11, 0]]  // j 2
        ],
        [
            // i 2
            [[4, 1], [3, 5], [7, 5]],  // j 0
            [[8, 1], [0, 0], [1, 5]],  // j 1
            [[16, 1], [5, 0], [2, 0]]  // j 2
        ]],
    [// face 2
        [
            // i 0
            [[7, 0], [21, 0], [38, 0]],  // j 0
            [[9, 0], [19, 0], [34, 3]],  // j 1
            [[14, 1], [20, 3], [36, 3]]  // j 2
        ],
        [
            // i 1
            [[3, 0], [13, 5], [29, 5]],  // j 0
            [[1, 0], [7, 0], [21, 0]],   // j 1
            [[6, 1], [9, 0], [19, 0]]    // j 2
        ],
        [
            // i 2
            [[4, 2], [12, 5], [26, 5]],  // j 0
            [[0, 1], [3, 0], [13, 5]],   // j 1
            [[2, 1], [1, 0], [7, 0]]     // j 2
        ]],
    [// face 3
        [
            // i 0
            [[26, 0], [42, 0], [58, 0]],  // j 0
            [[29, 0], [43, 0], [62, 3]],  // j 1
            [[38, 1], [47, 3], [64, 3]]   // j 2
        ],
        [
            // i 1
            [[12, 0], [28, 5], [44, 5]],  // j 0
            [[13, 0], [26, 0], [42, 0]],  // j 1
            [[21, 1], [29, 0], [43, 0]]   // j 2
        ],
        [
            // i 2
            [[4, 3], [15, 5], [31, 5]],  // j 0
            [[3, 1], [12, 0], [28, 5]],  // j 1
            [[7, 1], [13, 0], [26, 0]]   // j 2
        ]],
    [// face 4
        [
            // i 0
            [[31, 0], [41, 0], [49, 0]],  // j 0
            [[44, 0], [53, 0], [61, 3]],  // j 1
            [[58, 1], [65, 3], [75, 3]]   // j 2
        ],
        [
            // i 1
            [[15, 0], [22, 5], [33, 5]],  // j 0
            [[28, 0], [31, 0], [41, 0]],  // j 1
            [[42, 1], [44, 0], [53, 0]]   // j 2
        ],
        [
            // i 2
            [[4, 4], [8, 5], [16, 5]],    // j 0
            [[12, 1], [15, 0], [22, 5]],  // j 1
            [[26, 1], [28, 0], [31, 0]]   // j 2
        ]],
    [// face 5
        [
            // i 0
            [[50, 0], [48, 0], [49, 3]],  // j 0
            [[32, 0], [30, 3], [33, 3]],  // j 1
            [[24, 3], [18, 3], [16, 3]]   // j 2
        ],
        [
            // i 1
            [[70, 0], [67, 0], [66, 3]],  // j 0
            [[52, 3], [50, 0], [48, 0]],  // j 1
            [[37, 3], [32, 0], [30, 3]]   // j 2
        ],
        [
            // i 2
            [[83, 0], [87, 3], [85, 3]],  // j 0
            [[74, 3], [70, 0], [67, 0]],  // j 1
            [[57, 1], [52, 3], [50, 0]]   // j 2
        ]],
    [// face 6
        [
            // i 0
            [[25, 0], [23, 0], [24, 3]],  // j 0
            [[17, 0], [11, 3], [10, 3]],  // j 1
            [[14, 3], [6, 3], [2, 3]]     // j 2
        ],
        [
            // i 1
            [[45, 0], [39, 0], [37, 3]],  // j 0
            [[35, 3], [25, 0], [23, 0]],  // j 1
            [[27, 3], [17, 0], [11, 3]]   // j 2
        ],
        [
            // i 2
            [[63, 0], [59, 3], [57, 3]],  // j 0
            [[56, 3], [45, 0], [39, 0]],  // j 1
            [[46, 3], [35, 3], [25, 0]]   // j 2
        ]],
    [// face 7
        [
            // i 0
            [[36, 0], [20, 0], [14, 3]],  // j 0
            [[34, 0], [19, 3], [9, 3]],   // j 1
            [[38, 3], [21, 3], [7, 3]]    // j 2
        ],
        [
            // i 1
            [[55, 0], [40, 0], [27, 3]],  // j 0
            [[54, 3], [36, 0], [20, 0]],  // j 1
            [[51, 3], [34, 0], [19, 3]]   // j 2
        ],
        [
            // i 2
            [[72, 0], [60, 3], [46, 3]],  // j 0
            [[73, 3], [55, 0], [40, 0]],  // j 1
            [[71, 3], [54, 3], [36, 0]]   // j 2
        ]],
    [// face 8
        [
            // i 0
            [[64, 0], [47, 0], [38, 3]],  // j 0
            [[62, 0], [43, 3], [29, 3]],  // j 1
            [[58, 3], [42, 3], [26, 3]]   // j 2
        ],
        [
            // i 1
            [[84, 0], [69, 0], [51, 3]],  // j 0
            [[82, 3], [64, 0], [47, 0]],  // j 1
            [[76, 3], [62, 0], [43, 3]]   // j 2
        ],
        [
            // i 2
            [[97, 0], [89, 3], [71, 3]],  // j 0
            [[98, 3], [84, 0], [69, 0]],  // j 1
            [[96, 3], [82, 3], [64, 0]]   // j 2
        ]],
    [// face 9
        [
            // i 0
            [[75, 0], [65, 0], [58, 3]],  // j 0
            [[61, 0], [53, 3], [44, 3]],  // j 1
            [[49, 3], [41, 3], [31, 3]]   // j 2
        ],
        [
            // i 1
            [[94, 0], [86, 0], [76, 3]],  // j 0
            [[81, 3], [75, 0], [65, 0]],  // j 1
            [[66, 3], [61, 0], [53, 3]]   // j 2
        ],
        [
            // i 2
            [[107, 0], [104, 3], [96, 3]],  // j 0
            [[101, 3], [94, 0], [86, 0]],   // j 1
            [[85, 3], [81, 3], [75, 0]]     // j 2
        ]],
    [// face 10
        [
            // i 0
            [[57, 0], [59, 0], [63, 3]],  // j 0
            [[74, 0], [78, 3], [79, 3]],  // j 1
            [[83, 3], [92, 3], [95, 3]]   // j 2
        ],
        [
            // i 1
            [[37, 0], [39, 3], [45, 3]],  // j 0
            [[52, 0], [57, 0], [59, 0]],  // j 1
            [[70, 3], [74, 0], [78, 3]]   // j 2
        ],
        [
            // i 2
            [[24, 0], [23, 3], [25, 3]],  // j 0
            [[32, 3], [37, 0], [39, 3]],  // j 1
            [[50, 3], [52, 0], [57, 0]]   // j 2
        ]],
    [// face 11
        [
            // i 0
            [[46, 0], [60, 0], [72, 3]],  // j 0
            [[56, 0], [68, 3], [80, 3]],  // j 1
            [[63, 3], [77, 3], [90, 3]]   // j 2
        ],
        [
            // i 1
            [[27, 0], [40, 3], [55, 3]],  // j 0
            [[35, 0], [46, 0], [60, 0]],  // j 1
            [[45, 3], [56, 0], [68, 3]]   // j 2
        ],
        [
            // i 2
            [[14, 0], [20, 3], [36, 3]],  // j 0
            [[17, 3], [27, 0], [40, 3]],  // j 1
            [[25, 3], [35, 0], [46, 0]]   // j 2
        ]],
    [// face 12
        [
            // i 0
            [[71, 0], [89, 0], [97, 3]],   // j 0
            [[73, 0], [91, 3], [103, 3]],  // j 1
            [[72, 3], [88, 3], [105, 3]]   // j 2
        ],
        [
            // i 1
            [[51, 0], [69, 3], [84, 3]],  // j 0
            [[54, 0], [71, 0], [89, 0]],  // j 1
            [[55, 3], [73, 0], [91, 3]]   // j 2
        ],
        [
            // i 2
            [[38, 0], [47, 3], [64, 3]],  // j 0
            [[34, 3], [51, 0], [69, 3]],  // j 1
            [[36, 3], [54, 0], [71, 0]]   // j 2
        ]],
    [// face 13
        [
            // i 0
            [[96, 0], [104, 0], [107, 3]],  // j 0
            [[98, 0], [110, 3], [115, 3]],  // j 1
            [[97, 3], [111, 3], [119, 3]]   // j 2
        ],
        [
            // i 1
            [[76, 0], [86, 3], [94, 3]],   // j 0
            [[82, 0], [96, 0], [104, 0]],  // j 1
            [[84, 3], [98, 0], [110, 3]]   // j 2
        ],
        [
            // i 2
            [[58, 0], [65, 3], [75, 3]],  // j 0
            [[62, 3], [76, 0], [86, 3]],  // j 1
            [[64, 3], [82, 0], [96, 0]]   // j 2
        ]],
    [// face 14
        [
            // i 0
            [[85, 0], [87, 0], [83, 3]],     // j 0
            [[101, 0], [102, 3], [100, 3]],  // j 1
            [[107, 3], [112, 3], [114, 3]]   // j 2
        ],
        [
            // i 1
            [[66, 0], [67, 3], [70, 3]],   // j 0
            [[81, 0], [85, 0], [87, 0]],   // j 1
            [[94, 3], [101, 0], [102, 3]]  // j 2
        ],
        [
            // i 2
            [[49, 0], [48, 3], [50, 3]],  // j 0
            [[61, 3], [66, 0], [67, 3]],  // j 1
            [[75, 3], [81, 0], [85, 0]]   // j 2
        ]],
    [// face 15
        [
            // i 0
            [[95, 0], [92, 0], [83, 0]],  // j 0
            [[79, 0], [78, 0], [74, 3]],  // j 1
            [[63, 1], [59, 3], [57, 3]]   // j 2
        ],
        [
            // i 1
            [[109, 0], [108, 0], [100, 5]],  // j 0
            [[93, 1], [95, 0], [92, 0]],     // j 1
            [[77, 1], [79, 0], [78, 0]]      // j 2
        ],
        [
            // i 2
            [[117, 4], [118, 5], [114, 5]],  // j 0
            [[106, 1], [109, 0], [108, 0]],  // j 1
            [[90, 1], [93, 1], [95, 0]]      // j 2
        ]],
    [// face 16
        [
            // i 0
            [[90, 0], [77, 0], [63, 0]],  // j 0
            [[80, 0], [68, 0], [56, 3]],  // j 1
            [[72, 1], [60, 3], [46, 3]]   // j 2
        ],
        [
            // i 1
            [[106, 0], [93, 0], [79, 5]],  // j 0
            [[99, 1], [90, 0], [77, 0]],   // j 1
            [[88, 1], [80, 0], [68, 0]]    // j 2
        ],
        [
            // i 2
            [[117, 3], [109, 5], [95, 5]],  // j 0
            [[113, 1], [106, 0], [93, 0]],  // j 1
            [[105, 1], [99, 1], [90, 0]]    // j 2
        ]],
    [// face 17
        [
            // i 0
            [[105, 0], [88, 0], [72, 0]],  // j 0
            [[103, 0], [91, 0], [73, 3]],  // j 1
            [[97, 1], [89, 3], [71, 3]]    // j 2
        ],
        [
            // i 1
            [[113, 0], [99, 0], [80, 5]],   // j 0
            [[116, 1], [105, 0], [88, 0]],  // j 1
            [[111, 1], [103, 0], [91, 0]]   // j 2
        ],
        [
            // i 2
            [[117, 2], [106, 5], [90, 5]],  // j 0
            [[121, 1], [113, 0], [99, 0]],  // j 1
            [[119, 1], [116, 1], [105, 0]]  // j 2
        ]],
    [// face 18
        [
            // i 0
            [[119, 0], [111, 0], [97, 0]],  // j 0
            [[115, 0], [110, 0], [98, 3]],  // j 1
            [[107, 1], [104, 3], [96, 3]]   // j 2
        ],
        [
            // i 1
            [[121, 0], [116, 0], [103, 5]],  // j 0
            [[120, 1], [119, 0], [111, 0]],  // j 1
            [[112, 1], [115, 0], [110, 0]]   // j 2
        ],
        [
            // i 2
            [[117, 1], [113, 5], [105, 5]],  // j 0
            [[118, 1], [121, 0], [116, 0]],  // j 1
            [[114, 1], [120, 1], [119, 0]]   // j 2
        ]],
    [// face 19
        [
            // i 0
            [[114, 0], [112, 0], [107, 0]],  // j 0
            [[100, 0], [102, 0], [101, 3]],  // j 1
            [[83, 1], [87, 3], [85, 3]]      // j 2
        ],
        [
            // i 1
            [[118, 0], [120, 0], [115, 5]],  // j 0
            [[108, 1], [114, 0], [112, 0]],  // j 1
            [[92, 1], [100, 0], [102, 0]]    // j 2
        ],
        [
            // i 2
            [[117, 0], [121, 5], [119, 5]],  // j 0
            [[109, 1], [118, 0], [120, 0]],  // j 1
            [[95, 1], [108, 1], [114, 0]]    // j 2
        ]
    ]
];

/** CoordIJK unit vectors corresponding to the 7 H3 digits. */
pub static UNIT_VECS: [CoordIJK; 7] = [
   CoordIJK{i: 0, j: 0, k: 0},  // direction 0
   CoordIJK{i: 0, j: 0, k: 1},  // direction 1
   CoordIJK{i: 0, j: 1, k: 0},  // direction 2
   CoordIJK{i: 0, j: 1, k: 1},  // direction 3
   CoordIJK{i: 1, j: 0, k: 0},  // direction 4
   CoordIJK{i: 1, j: 0, k: 1},  // direction 5
   CoordIJK{i: 1, j: 1, k: 0}   // direction 6
];

pub const H3_DIGIT_MASK: u64 = 7;
pub const H3_DIGIT_MASK_NEGATIVE: u64 = !H3_DIGIT_MASK;
pub const MAX_H3_RES: u64 = 15;
pub const H3_PER_DIGIT_OFFSET: u64 = 3;


pub const NUM_BASE_CELLS: usize = 122;

/** Resolution 0 base cell data table.
For each base cell, gives:
 * The "home" face and ijk+ coordinates on that face
 * Whether or not the base cell is a pentagon.
 * Additionally, if the base cell is a pentagon,
   the two cw offset rotation adjacent faces are given if available.
 */
pub static BASE_CELL_DATA: [BaseCellData; NUM_BASE_CELLS] = [
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 0
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 1
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 2
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 3
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: None},            // base cell 4
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 5
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 6
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 7
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 8
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 9
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 10
   BaseCellData{home_fijk: FaceIJK{face: 1,  coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 11
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 12
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 13
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((2, 6))},    // base cell 14
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 15
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 16
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 17
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 18
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 19
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 20
   BaseCellData{home_fijk: FaceIJK{face: 2,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 21
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 22
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 23
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((1, 5))},    // base cell 24
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 25
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 26
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 27
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 28
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 29
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 30
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 31
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 32
   BaseCellData{home_fijk: FaceIJK{face: 0,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 33
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 34
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 35
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 36
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 37
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((3, 7))},    // base cell 38
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 39
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 40
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 41
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 42
   BaseCellData{home_fijk: FaceIJK{face: 3,  coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 43
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 44
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 45
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 46
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 47
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 48
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((0, 9))},    // base cell 49
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 50
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 51
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 52
   BaseCellData{home_fijk: FaceIJK{face: 4,  coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 53
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 54
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 55
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 56
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 57
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((4, 8))},    // base cell 58
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 59
   BaseCellData{home_fijk: FaceIJK{face: 11, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 60
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 61
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 62
   BaseCellData{home_fijk: FaceIJK{face: 6,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((11, 15))},  // base cell 63
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 64
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 65
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 66
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 67
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 68
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 69
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 70
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 71
   BaseCellData{home_fijk: FaceIJK{face: 7,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((12, 16))},  // base cell 72
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 73
   BaseCellData{home_fijk: FaceIJK{face: 10, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 74
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 75
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 76
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 77
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 78
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 79
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 80
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 81
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 1, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 82
   BaseCellData{home_fijk: FaceIJK{face: 5,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((10, 19))},  // base cell 83
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 84
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 85
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 86
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 87
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 88
   BaseCellData{home_fijk: FaceIJK{face: 12, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 89
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 90
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 91
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 92
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 93
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 94
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 95
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 96
   BaseCellData{home_fijk: FaceIJK{face: 8,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((13, 17))},  // base cell 97
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 98
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 99
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 100
   BaseCellData{home_fijk: FaceIJK{face: 14, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 101
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 102
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 103
   BaseCellData{home_fijk: FaceIJK{face: 13, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 104
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 105
   BaseCellData{home_fijk: FaceIJK{face: 16, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 106
   BaseCellData{home_fijk: FaceIJK{face: 9,  coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: Some((14, 18))},  // base cell 107
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 108
   BaseCellData{home_fijk: FaceIJK{face: 15, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 109
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 0, j: 1, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 110
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 111
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 0, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 112
   BaseCellData{home_fijk: FaceIJK{face: 17, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 113
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 114
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 0, j: 1, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 115
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 116
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 2, j: 0, k: 0}}, is_pentagon: true , pentagon_cw_offset_faces: None},  // base cell 117
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 118
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 0, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 119
   BaseCellData{home_fijk: FaceIJK{face: 19, coord: CoordIJK{i: 1, j: 0, k: 1}}, is_pentagon: false, pentagon_cw_offset_faces: None},            // base cell 120
   BaseCellData{home_fijk: FaceIJK{face: 18, coord: CoordIJK{i: 1, j: 0, k: 0}}, is_pentagon: false, pentagon_cw_offset_faces: None}            // base cell 121
];
