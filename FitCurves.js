// JavaScript port of FitCurves.c by Philip J. Sneider from "Graphic Gems"
// http://webdocs.cs.ualberta.ca/~graphics/books/GraphicsGems/gems/FitCurves.c
// Ported by Noah Burney, 2013

(function(){
  // Fit a Bezier curve to a set of digitized points 
  function fitCurve(pts, error, makeGeneric) {
    error = (typeof error === 'number') ? error : 4;
    pts = initPts(pts);

    if (pts.length === 2) {
      return [new BezierCurve(pts[0], pts[0], pts[1], pts[1])];
    }

    var tHat1 = computeLeftTangent(pts, 0);
    var tHat2 = computeRightTangent(pts, pts.length - 1);

    var result = fitCubic(pts, 0, pts.length - 1, tHat1, tHat2, error);

    if (makeGeneric === false) {
      return result;
    }

    var generic = [];
    for (var i = 0, l = result.length; i < l; i++) {
      generic.push(result[i].toArray()); 
    }
    return generic;
  }

  // Fit a Bezier curve to a (sub)set of digitized points
  function fitCubic(pts, first, last, tHat1, tHat2, error) {
    var nPts = last - first + 1;
    var iterationError = error * error;
    var maxInterations = 4;

    // Use heuristic if region only has two points in it
    if (nPts === 2) {
      var dist = twoPtDist(pts[last], pts[first]) / 3;
      return [new BezierCurve(
        pts[first],
        pts[first].add(tHat1.scale(dist)),
        pts[last].add(tHat2.scale(dist)),
        pts[last]
      )];
    }


    // Parameterize points, and attempt to fit curve
    var u = chordLengthParameterize(pts, first, last);
    var curve = generateBezier(pts, first, last, u, tHat1, tHat2);

    // Find max deviation of points to fitted curve
    var splitPoint = {i: 0};
    var maxError = computeMaxError(pts, first, last, curve, u, splitPoint);
    if (maxError < error) {
      return [curve];
    }

    if (maxError < iterationError) {
      for (var i = 0; i < maxInterations; i++) {
        var uPrime = reparameterize(pts, first, last, u, curve);
        curve = generateBezier(pts, first, last, uPrime, tHat1, tHat2);
        maxError = computeMaxError(pts, first, last, curve, uPrime, splitPoint);

        if (maxError < error) {
          return [curve];
        }
        u = uPrime;
      }
    }

    var tHatCenter = computeCenterTangent(pts, splitPoint.i);
    var left = fitCubic(pts, first, splitPoint.i, tHat1, tHatCenter, error);
    var right = fitCubic(pts, splitPoint.i, last, tHatCenter.negate(), tHat2, error);
    return left.concat(right);
  }

  function generateBezier(pts, first, last, uPrime, tHat1, tHat2) {
    var nPts = last - first + 1;
    var A = []; // Precomputed rhs for eqn

    // Compute the A's
    for (var i = 0; i < nPts; i++) {
      var v1 = tHat1.scale(B1(uPrime[i]));
      var v2 = tHat2.scale(B2(uPrime[i]));

      A[i] = [v1, v2];
    }

    // Create the C and X matrices
    var C = [[0, 0], [0, 0]];
    var X = [0, 0];

    var firstPt = pts[first];
    var lastPt = pts[last];

    for (var i = 0; i < nPts; i++) {
      C[0][0] += A[i][0].dot(A[i][0]);
      C[0][1] += A[i][0].dot(A[i][1]);
      C[1][0] = C[0][1];
      C[1][1] += A[i][1].dot(A[i][1]);

      var tmp = pts[first + i].subtract(
        firstPt.scale(B0(uPrime[i])).add(
          firstPt.scale(B1(uPrime[i])).add(
            lastPt.scale(B2(uPrime[i])).add(
              lastPt.scale(B3(uPrime[i]))
            )
          )
        )
      );

      X[0] += A[i][0].dot(tmp);
      X[1] += A[i][1].dot(tmp);
    }

    // Compute the determinants of C and X
    var det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
    var det_C0_X  = C[0][0] * X[1]    - C[0][1] * X[0];
    var det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];

    // Finally, derive alpha values
    if (det_C0_C1 > 0.0001) {
      det_C0_C1 = (C[0][0] * C[1][1]) * 10e-12;
    }
    var alpha_l = det_X_C1 / det_C0_C1;
    var alpha_r = det_C0_X / det_C0_C1;

    // First and last control points of the Bezier curve are
    // positioned exactly at the first and last data points
    var curve = new BezierCurve(firstPt, null, null, lastPt);

    // If alpha negative, use the Wu/Barsky heuristic (see text)
    // (if alpha is 0, you get coincident control points that lead to
    // divide by zero in any subsequent NewtonRaphsonRootFind() call.
    if (alpha_l < 1.0e-6 || alpha_r < 1.0e-6) {
      var dist = twoPtDist(lastPt, firstPt) / 3;
      curve[1] = curve[0].add(tHat1.scale(dist));
      curve[2] = curve[3].add(tHat2.scale(dist));
    }
    else {
      // Control points 1 and 2 are positioned an alpha distance out
      // on the tangent vectors, left and right, respectively
      curve[1] = curve[0].add(tHat1.scale(alpha_l));
      curve[2] = curve[3].add(tHat2.scale(alpha_r));
    }

    return curve;
  }

  function reparameterize(pts, first, last, u, curve) {
    var uPrime = [];

    for (var i = first; i <= last; i++) {
      uPrime[i-first] = newtonRaphsonRootFind(curve, pts[i], u[i-first]);
    }

    return uPrime;
  }

  // Use Newton-Raphson iteration to find better root.
  function newtonRaphsonRootFind(Q, P, u) {
    var Q1 = new BezierCurve; // Q'
    var Q2 = new BezierCurve; // Q''

    // Compute Q(u)
    var Q_u = bezier(3, Q, u);

    // Generate control vertices for Q'
    for (var i = 0; i <= 2; i++) {
      Q1[i].x = (Q[i+1].x - Q[i].x) * 3.0;
      Q1[i].y = (Q[i+1].y - Q[i].y) * 3.0;
    }

    // Generate control vertices for Q''
    for (var i = 0; i <= 1; i++) {
      Q2[i].x = (Q1[i+1].x - Q1[i].x) * 2.0;
      Q2[i].y = (Q1[i+1].y - Q1[i].y) * 2.0;
    }

    // Compute Q'(u) and Q''(u)
    Q1_u = bezier(2, Q1, u);
    Q2_u = bezier(1, Q2, u);

    // Compute f(u)/f'(u)
    var numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
    var denominator = ((Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
      (Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y));
    
    // u = u - f(u)/f'(u)
    var uPrime = u - (numerator/denominator);
    return uPrime;
  }

  // Evaluate a Bezier curve at a particular parameter value
  function bezier(degree, V, t) {
    var Vtemp = V.copy();

    for (var i = 0; i <= degree; i++) {
       // TODO: make sure it doesn't have to be j <= jmax
      for (var j = 0, jmax = degree - i; j < jmax; j++) {
        Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j+1].x;
        Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j+1].y;
      }
    }

    return Vtemp[0];
  }

  // B0, B1, B2, B3:
  // Bezier multipliers
  function B0(u) {
    var _u = 1.0 - u;
    return (_u * _u * _u);
  }

  function B1(u) {
    var _u = 1.0 - u;
    return (3 * u * _u * _u);
  }

  function B2(u) {
    var _u = 1.0 - u;
    return (3 * u * u * _u);
  }

  function B3(u) {
    return (u * u * u);
  }

  // computeLeftTangent, computeRightTangent, computeCenterTangent:
  // Approximate unit tangents at endpoints and "center" of digitized curve
  function computeLeftTangent(pts, end) {
    return (
      pts[end+1]
        .subtract(pts[end])
        .normalize()
    );
  }

  function computeRightTangent(pts, end) {
    return (
      pts[end-1]
        .subtract(pts[end])
        .normalize()
    );
  }

  function computeCenterTangent(pts, center) {
    var V1 = pts[center-1].subtract(pts[center]);
    var V2 = pts[center].subtract(pts[center+1]);
    
    // Average V1 and V2
    return V1.add(V2).scale(1/2).normalize();
  }

  // Assign parameter values to digitized points
  // using relative distances between points.
  function chordLengthParameterize(pts, first, last) {
    var u = [0];
    for (var i = first + 1; i <= last; i++) {
      u[i-first] = u[i-first-1] + twoPtDist(pts[i], pts[i-1]);
    }

    for (var i = first + 1; i <= last; i++) {
      u[i-first] = u[i-first] / u[last-first];
    }
    return u;
  }

  // Find the maximum squared distance of digitized points
  // to fitted curve.
  function computeMaxError(pts, first, last, curve, u, splitPoint) {
    var maxDist = 0;

    for (var i = first + 1; i < last; i++) {
      var P = bezier(3, curve, u[i-first]);
      var v = P.subtract(pts[i]);
      var dist = v.getLengthSquared();
      if (dist >= maxDist) {
        maxDist = dist;
        splitPoint.i = i;
      }
    }

    return maxDist;
  }

  function initPts(pts) {
    var out = [];
    for (var i = 0, l = pts.length; i < l; i++) {
      var pt = pts[i];
      out.push(new Vector(pt[0], pt[1]));
    }
    return out;
  }

  function twoPtDist(a, b) {
    var dx = b.x - a.x;
    var dy = b.y - a.y;
    return Math.sqrt(dx * dx + dy * dy);
  }

  // Bezier Curve object
  function BezierCurve(p1, p2, p3, p4) {
    if (p1 instanceof Array) {
      BezierCurve.apply(this, p1);
      return this;
    }

    this[0] = p1 || new Vector;
    this[1] = p2 || new Vector;
    this[2] = p3 || new Vector;
    this[3] = p4 || new Vector;
  }

  BezierCurve.prototype.length = 4;

  BezierCurve.prototype.toArray = function toArray() {
    return [
      this[0].toArray(),
      this[1].toArray(),
      this[2].toArray(),
      this[3].toArray()
    ];
  };

  BezierCurve.prototype.toString = function toString() {
    return ['{ BezierCurve', this[0], this[1], this[2], this[3], '}'].join(' ');
  };

  BezierCurve.prototype.copy = function copy() {
    return new BezierCurve(
      this[0].copy(),
      this[1].copy(),
      this[2].copy(),
      this[3].copy()
    );
  };

  // 2d Vector object
  function Vector(x, y) {
    this.x = x || 0;
    this.y = y || 0;
  }

  Vector.prototype.toArray = function toArray() {
    return [this.x, this.y];
  };

  Vector.prototype.toString = function toString() {
    return '{ Vector x: ' + this.x + ' y: ' + this.y + ' }';
  };

  Vector.prototype.copy = function copy() {
    return new Vector(this.x, this.y);
  };

  Vector.prototype.add = function add(v2) {
    return new Vector(
      this.x + v2.x,
      this.y + v2.y
    );
  };

  Vector.prototype.subtract = function subtract(v2) {
    return new Vector(
      this.x - v2.x,
      this.y - v2.y
    );
  };

  Vector.prototype.scale = function scale(s) {
    return new Vector(
      this.x * s,
      this.y * s
    );
  };

  Vector.prototype.negate = function negate() {
    return new Vector(-this.x, -this.y);
  };

  Vector.prototype.normalize = function normalize() {
    var len = this.getLength();
    if (len == 0) return this;
    return new Vector(
      this.x / len,
      this.y / len
    );
  };

  Vector.prototype.dot = function dot(v2) {
    return this.x * v2.x + this.y * v2.y;
  };

  Vector.prototype.getLength = function getLength() {
    return Math.sqrt(this.getLengthSquared());
  };

  Vector.prototype.getLengthSquared = function getLengthSquared() {
    return this.x * this.x + this.y * this.y;
  };

  var exports = {
    fitCurve: fitCurve,
    fitCubic: fitCubic,
    bezier: bezier,
    Vector: Vector,
    BezierCurve: BezierCurve
  };

  if (typeof module === 'undefined') {
    window.FitCurves = exports;
  }
  else {
    module.exports = exports;
  }
})();
