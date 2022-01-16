const { inv } = require("mathjs");

function SumVecSpecial(a, b) {
  /**
   * Sum vectors where the first vector is a column vector of type a = [[1],[2]]
   * @param {Array} a first array of arrays of length 1
   * @param {Array} b second array of numbers
   * @return {Array}  array where each element is calculated as a + b
   */
  const len = a.length;
  let out = Array(len);
  for (let col = 0; col < len; col++) {
    out[col] = a[col][0] + b[col];
  }
  return out;
}

function diffVec(a, b) {
  /**
   * Element-wise subtracts vector b from vector a
   * @param {Array} a first array of numbers.
   * @param {Array} b second array of numbers.
   * @return {Array}  array where each element is calculated as a minus b
   */
  const len = a.length;
  let out = Array(len);
  for (let col = 0; col < len; col++) {
    out[col] = a[col] - b[col];
  }
  return out;
}

function multiply(A, B) {
  /**
   * Matrix multiplication of matrices A and B so A*B
   * @param {Array} A Array of arrays where the len of outer array is the number of rows and the inner aray is the number of columnet. Ex. A=[[1,2],[3,4]]
   * @param {Array} B Array of arrays where the len of outer array is the number of rows and the inner aray is the number of columnet. Ex. B=[[5,6],[7,8]]
   * @return {Array}  Array of arrays calculated as a matrix product of A and B; A*B; Ex. out = [[19, 22], [43, 50]]
   */
  const aNumRows = A.length,
    aNumCols = A[0].length,
    bNumCols = B[0].length,
    out = new Array(aNumRows); // initialize array of rows
  for (let row = 0; row < aNumRows; ++row) {
    out[row] = new Array(bNumCols); // initialize the current row
    for (let col = 0; col < bNumCols; ++col) {
      out[row][col] = 0; // initialize the current cell
      for (let i = 0; i < aNumCols; ++i) {
        out[row][col] += A[row][i] * B[i][col];
      }
    }
  }
  return out;
}

function multiplyMatVec(A, b) {
  /**
   * Matrix multiplication of a matrix A and a vector b so A*b
   * @param {Array} A Array of arrays where the len of outer array is the number of rows and the inner aray is the number of columnet. Ex. A=[[1,2],[3,4]]
   * @param {Array} b Array of numbers Ex. B=[5,6]
   * @return {Array}  Array of arrays calculated as a matrix product of A and b; A*b; Ex. out = [[17],[39]]
   */
  let aNumRows = A.length,
    aNumCols = A[0].length,
    out = new Array(aNumRows); // initialize array of rows
  for (let row = 0; row < aNumRows; ++row) {
    out[row] = new Array(1); // initialize the current row
    out[row][0] = 0; // initialize the current cell
    for (let i = 0; i < aNumCols; ++i) {
      out[row][0] += A[row][i] * b[i];
    }
  }
  return out;
}

function constructIdentity(dim) {
  /**
   * Construct an identity matrix where every element is equal to 0 except the diagonal elements that have 1
   * @param {Int} dimension of the identity matrix. Ex dim = 2
   * @return {Array} returns an array of arrays representing the identity matrix out = [[1,0],[0,1]]
   */
  const out = [];
  for (let row = 0; row < dim; row++) {
    if (!out[row]) {
      out[row] = [];
    }
    for (let col = 0; col < dim; col++) {
      if (row === col) {
        out[row][col] = 1;
      } else {
        out[row][col] = 0;
      }
    }
  }
  return out;
}

function constructZeros(nrow, ncol) {
  /**
   * Construct a matrix with nrow rows and ncol columns where every element is equal to 0
   * @param {Int} nrow number of rows of the identity matrix. Ex nrow = 2
   * @param {Int} ncol number of columns of the identity matrix. Ex nrow = 3
   * @return {Array} returns an array of arrays representing the zero matrix out = [[0,0,0],[0,0,0]]
   */

  let out = [];
  for (let row = 0; row < nrow; row++) {
    if (!out[row]) {
      out[row] = [];
    }
    for (let col = 0; col < ncol; col++) {
      out[row][col] = 0;
    }
  }
  return out;
}

function SWHeart(u, v, alpha) {
  /**
   * Calculate the matrix H (Heart of the Wilson function)
   * @param {Array} u  n_1 x 1 array arrays of maturities. Ex. u = [[1], [3]]
   * @param {Array} v  n_2 x 1 array of arrays of maturities. Ex. v = [[1], [2], [3], [5]]
   * @param {Float} alpha  floating number representing the convergence speed parameter alpha. Ex. alpha = 0.05
   * @returns {Array}  n_1 x n_2 array of arrays representing the Heart of the Wilson function for selected maturities and parameter alpha. H is calculated as in the paragraph 132 of the EIOPA documentation.
   * For more information see https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/12092019-technical_documentation.pdf
   */

  const nrow = u.length,
    ncol = v.length;
  let H = new Array(nrow);

  // For each element, calculate heart of W
  for (let row = 0; row < nrow; row++) {
    H[row] = new Array(ncol);
    for (let col = 0; col < ncol; col++) {
      H[row][col] =
        0.5 *
        (alpha * (u[row][0] + v[col][0]) +
          Math.exp(-alpha * (u[row][0] + v[col][0])) -
          alpha * Math.abs(u[row][0] - v[col][0]) -
          Math.exp(-alpha * Math.abs(u[row][0] - v[col][0])));
    }
  }
  return H;
}

function SWCalibrate(r_Obs, M_Obs, ufr, alpha) {
  /**
   * Calculate the calibration vector "b" using a Smith-Wilson algorithm
   * @param {Array} r_Obs n x 1 array of arrays representing rates, for which you wish to calibrate the algorithm. Each rate belongs to an observable zero coupon bond with a known maturity. Ex. r_Obs = [[0.0024], [0.0034]]
   * @param {Array} M_Obs m x 1 array of arrays representing maturities of bonds, that have rates provided in input (r_Obs). Ex. M_Obs=[[1], [3]]
   * @param {Float} ufr floating number, representing the ultimate forward rate. Ex. ufr = 0.042
   * @param {Float} alpha floating number representing the convergence speed parameter alpha. Ex. alpha = 0.05
   * @return {Array} n x 1 array of arrays representing the calibration vector "b" needed to interpolate and extrapolate b =[[14], [-21]]
   *  For more information see https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/12092019-technical_documentation.pdf
   */

  const len = r_Obs.length,
    C = constructIdentity(len);

  let b = [],
    p = Array(len),
    d = Array(len),
    q = Array(len),
    Q = constructZeros(len, len);

  for (let col = 0; col < len; col++) {
    p[col] = Math.pow(1 + r_Obs[col][0], -M_Obs[col][0]); //    p = (1+r).^(-M);
    d[col] = Math.exp(-Math.log(1 + ufr) * M_Obs[col][0]); //    d = exp(-log(1+ufr) .* M);
    Q[col][col] = C[col][col] * d[col]; //    Q = diag(d) * C;
  }

  q = multiplyMatVec(C, d); //    q = C'*d;
  let H = SWHeart(M_Obs, M_Obs, alpha); // Calculate SWHeart
  return multiplyMatVec(inv(multiply(multiply(Q, H), Q)), diffVec(p, q)); //    b = (Q' * H * Q)\(p-q);
}

function SWExtrapolate(M_Tar, M_Obs, b, ufr, alpha) {
  /**
   * Interpolate or/and extrapolate rates for targeted maturities using a Smith-Wilson algorithm.
   * @param {Array} M_Target k x 1 array of arrays. Each element represents a bond maturity of interest. Ex. M_Tar = [[1], [2], [3], [5]]
   * @param {Array} M_Obs n x 1 array of arrays. Observed bond maturities used for calibrating the calibration vector b. Ex. M_Obs = [[1], [3]]
   * @param {Array}  b  n x 1 array of arrays. Calibration vector "b" calculated on observed bonds.
   * @param {Float} ufr floating number, representing the ultimate forward rate. Ex. ufr = 0.042
   * @param {Float} alpha floating number representing the convergence speed parameter alpha. Ex. alpha = 0.05
   * @return {Array} k x 1 array of arrays. Represents the targeted rates for a zero-coupon bond. Each rate belongs to a targeted zero-coupon bond with a maturity from T_Target. Ex. r = [[0.0024], [0.0029], [0.0034], [0.0039]]
   * For more information see https://www.eiopa.europa.eu/sites/default/files/risk_free_interest_rate/12092019-technical_documentation.pdf
   */

  const obsLen = M_Obs.length,
    tarLen = M_Tar.length,
    C = constructIdentity(obsLen);
  let d = Array(obsLen),
    Q = constructZeros(obsLen, obsLen),
    expom = Array(tarLen),
    dDiag = constructZeros(tarLen, tarLen),
    r = Array(tarLen);

  for (let col = 0; col < obsLen; col++) {
    d[col] = Math.exp(-Math.log(1 + ufr) * M_Obs[col][0]);
    Q[col][col] = C[col][col] * d[col];
  }
  let H = SWHeart(M_Tar, M_Obs, alpha);

  for (let col = 0; col < tarLen; col++) {
    dDiag[col][col] = Math.exp(-Math.log(1 + ufr) * M_Tar[col][0]);
    expom[col] = Math.exp(-Math.log(1 + ufr) * M_Tar[col][0]);
  }

  let p = SumVecSpecial(
    multiplyMatVec(multiply(multiply(dDiag, H), Q), b),
    expom
  );
  for (let col = 0; col < tarLen; col++) {
    r[col] = [Math.pow(p[col], -1 / M_Tar[col]) - 1];
  }
  return r;
}

var M_Obs = [
  [1],
  [2],
  [3],
  [4],
  [5],
  [6],
  [7],
  [8],
  [9],
  [10],
  [11],
  [12],
  [13],
  [14],
  [15],
  [16],
  [17],
  [18],
  [19],
  [20],
];
var r_Obs = [
  [0.0131074591432979],
  [0.0222629098372424],
  [0.0273403667327403],
  [0.0317884414257146],
  [0.0327205345299401],
  [0.0332867589595655],
  [0.0336112121443886],
  [0.0341947663149128],
  [0.0345165922380981],
  [0.0346854377006694],
  [0.035717334079127],
  [0.0368501673784445],
  [0.0376263620230677],
  [0.0385237084707761],
  [0.0395043823351044],
  [0.0401574909803133],
  [0.0405715278625131],
  [0.0415574765441695],
  [0.0415582458410996],
  [0.042551132694631],
];
var ufr = 0.042;
var alpha = 0.142068;

var M_Tar = [
  [1],
  [2],
  [3],
  [4],
  [5],
  [6],
  [7],
  [8],
  [9],
  [10],
  [11],
  [12],
  [13],
  [14],
  [15],
  [16],
  [17],
  [18],
  [19],
  [20],
  [21],
  [22],
  [23],
  [24],
  [25],
  [26],
  [27],
  [28],
  [29],
  [30],
  [31],
  [32],
  [33],
  [34],
  [35],
  [36],
  [37],
  [38],
  [39],
  [40],
  [41],
  [42],
  [43],
  [44],
  [45],
  [46],
  [47],
  [48],
  [49],
  [50],
  [51],
  [52],
  [53],
  [54],
  [55],
  [56],
  [57],
  [58],
  [59],
  [60],
  [61],
  [62],
  [63],
  [64],
  [65],
];

// Example of use ////////////////////////
var b = SWCalibrate(r_Obs, M_Obs, ufr, alpha);

console.table(b);
var r = SWExtrapolate(M_Tar, M_Obs, b, ufr, alpha);
console.table(r);
