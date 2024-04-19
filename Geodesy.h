
//  Copyright (c) 2023-2024, Andre Caceres Carrilho
//
//  Geodesy.h common geodesic routines 
//  tags header-only, thread-safe, portable, standard-cpp
//    utm, universal transverse mercator, geodesy, geodesic
//

#ifndef _Geodesy_h
#define _Geodesy_h 1

#if __cplusplus < 201703L
    #error At least C++17 standard is need to build this code
#endif 

#include <cmath>
#include <limits>

namespace Geo {
namespace Math 
{
  template<typename Real>
  inline Real degToRad(const Real& x){
    return x * Real(3.1415926535897932384626) / Real(180);
  }
  template<typename Real>
  inline Real radToDeg(const Real& x){
    return x * Real(180) / Real(3.1415926535897932384626);
  }
  template<typename Real>
  inline Real sqr(const Real& x) {
    return x * x;
  }
}

enum SPHEROID {
    SPHEROID_GRS67,
    SPHEROID_SAD69,
    SPHEROID_WGS72,
    SPHEROID_GRS80,
    SPHEROID_WGS84
};

// Templated because you might want to use long double.
// The term 'Spheroid' is what Open Geospatial Consortium uses.
template<typename Real = double> 
class Spheroid{
private:
  Real majorSemiAxis; // equatorial radius (m)
  Real minorSemiAxis; // polar semi axis (m)
  Real flattening;
  Real eccentricitySqr;
  Real eccentricity;

  // Coefficients for Karney-Kruger UTM equations: 
  // Deakin, R.E.; Hunter, M.N.; Karney, C.F.F. A fresh look at the UTM projection.
  Real rectifyingRadius;
  Real alpha[8], beta[8];

public:
  Spheroid(             
    Real majorSemiAxis_ = Real(6378137), // defaults to WGS-84 
    Real flattening_    = Real(1) / Real(298.257223563)) 
  {
    setFrom(majorSemiAxis_, flattening_);
  }
  Spheroid(SPHEROID ellps) { set(ellps); }
  void set(SPHEROID ellps);
  void setFrom(Real majorSemiAxis_, Real flattening_);

  // Conversion between Geodetic and Geocentric Cartesian coordinates,
  // which is also known as Earth Centered, Earth Fixed (ECEF)
  void geodToEcef(
    const Real& latDeg,
    const Real& lonDeg,
    const Real& hGeo, 
    Real& ecefX, 
    Real& ecefY, 
    Real& ecefZ) const;

  void ecefToGeod(
    const Real& ecefX, 
    const Real& ecefY, 
    const Real& ecefZ, 
    Real& latDeg, 
    Real& lonDeg, 
    Real& hGeo) const;

  // UTM projection (forward and backward)
  // Uses Karney-Kruger equations (less than 5 nm error)
  void geoToUtm(
    const Real& latDeg,
    const Real& lonDeg,
    Real& utmE,
    Real& utmN,
    int& utmZone,
    char& hemisphere) const;

  void utmToGeo(
    const Real& utmE,
    const Real& utmN,
    const int& utmZone,
    const char& hemisphere,
    Real& latDeg,
    Real& lonDeg) const;

  Real utmScaleFactor(const Real& lonDeg) const;

  // Inverse Geodesic problem 
  // Vincenty equations (less than 0.5 mm error) 
  // Note: might fail for near antipodal points
  Real distance(
    const Real& latDeg1, 
    const Real& lonDeg1, 
    const Real& latDeg2, 
    const Real& lonDeg2
  ) const;

private:
  inline Real conformalLatitude(const Real& phiRad) const;
  inline Real conformalLatitude(const Real& u_a, const Real& v_a) const;
  inline Real radiusDenom(const Real& phiRad) const;
  inline Real radiusPrime(const Real& phiRad) const;
  inline Real radiusMeridional(const Real& phiRad) const;
  inline Real utmCentralMeridian(const Real& lonDeg, int &utmZone) const;
  void computeUtmParams();
};

template<typename Real> 
void Spheroid<Real>::set(SPHEROID ellps)
{
  switch(ellps)
  {
    case SPHEROID_GRS67: 
        setFrom(Real(6378160), Real(1) / Real(298.247167427));
        break;
    case SPHEROID_SAD69: 
        setFrom(Real(6378160), Real(1) / Real(298.25));
        break;
    case SPHEROID_WGS72: 
        setFrom(Real(6378135), Real(1) / Real(298.26));
        break;
    case SPHEROID_GRS80:
        setFrom(Real(6378137), Real(1) / Real(298.2572221009));
        break;
    case SPHEROID_WGS84:
        setFrom(Real(6378137), Real(1) / Real(298.257223563));
        break;
  }
}

template<typename Real> 
void Spheroid<Real>::setFrom(Real majorSemiAxis_, Real flattening_){
  majorSemiAxis = majorSemiAxis_;
  flattening    = flattening_;

  eccentricitySqr  = flattening * (Real(2) - flattening);
  eccentricity     = std::sqrt(eccentricitySqr);
  minorSemiAxis    = majorSemiAxis * std::sqrt(Real(1) - eccentricitySqr);
  computeUtmParams();
}

template<typename Real> 
void Spheroid<Real>::geodToEcef(
  const Real& latDeg,
  const Real& lonDeg,
  const Real& hGeo, 
  Real& ecefX, 
  Real& ecefY, 
  Real& ecefZ
) const{
  const Real phiRad = Math::degToRad(latDeg);
  const Real lmbRad = Math::degToRad(lonDeg);
  const Real radiusN = radiusPrime(phiRad);
  ecefX = (radiusN + hGeo) * std::cos(phiRad) * std::cos(lmbRad);
  ecefY = (radiusN + hGeo) * std::cos(phiRad) * std::sin(lmbRad);
  ecefZ = (radiusN * (Real(1) - eccentricitySqr) + hGeo) * std::sin(phiRad);
}

template<typename Real> 
void Spheroid<Real>::ecefToGeod(
  const Real& ecefX, 
  const Real& ecefY, 
  const Real& ecefZ, 
  Real& latDeg, 
  Real& lonDeg, 
  Real& hGeo
) const{
  const Real rho = std::hypot(ecefX, ecefY);
  const Real eps = std::numeric_limits<Real>::epsilon();

  Real radiusN = majorSemiAxis;
  Real phiSin;
  Real zi = ecefZ;
  Real zk = Real(0);

  for (int itrCount = 0; 
    // defining an upper bound for the loop 
    itrCount++ < 25 && std::fabs(zi- zk) > eps; )
  {
    zk = zi;
    phiSin = zi / std::hypot(rho, zi);
    radiusN = majorSemiAxis / 
      std::sqrt(Real(1) - eccentricitySqr * Math::sqr(phiSin));
    zi = ecefZ + radiusN * eccentricitySqr * phiSin;
  }

  Real lonRad = Real(0);
  Real latRad = Real(0);

  if (rho > eps){
    latRad = std::atan(zi / rho);
    lonRad = std::atan2(ecefY, ecefX);
  }
  else{
    if (ecefZ >= Real(0))
      latRad = Real( 0.5) * Real(3.1415926535897932384626);
    else 
      latRad = Real(-0.5) * Real(3.1415926535897932384626);
  }

  latDeg = Math::radToDeg(latRad);
  lonDeg = Math::radToDeg(lonRad);
  hGeo = std::hypot(rho, zi) - radiusN;
}

template<typename Real> 
inline Real Spheroid<Real>::radiusDenom(const Real& phiRad) const{
  return Real(1) - eccentricitySqr * Math::sqr(std::sin(phiRad));
}
template<typename Real> 
inline Real Spheroid<Real>::radiusPrime(const Real& phiRad) const{
  return majorSemiAxis / std::sqrt( radiusDenom(phiRad) );
}
template<typename Real> 
inline Real Spheroid<Real>::radiusMeridional(const Real& phiRad) const{
  return (majorSemiAxis * (Real(1) - eccentricitySqr)) / 
    std::pow(radiusDenom(phiRad) , Real(1.5));
}  

///////////////////////////////////////////////////////////////////////////////////////
// Karney-Kruger UTM equations from: 
// Deakin, R.E.; Hunter, M.N.; Karney, C.F.F. A fresh look at the UTM projection.
// Very accurate, less than 5 nanometer error

template<typename Real> 
inline Real 
Spheroid<Real>::conformalLatitude(const Real& latRad) const
{
// Equations (88) and (89)
// Used on the Forward projection (geoToUtm)
  const Real x = std::tan(latRad);
  const Real y = x * x;
  const Real z = (eccentricity * x) / std::sqrt(Real(1) + y);
  const Real sigma = std::sinh(eccentricity * std::atanh(z));
  const Real r = x * std::sqrt(Real(1) + sigma * sigma);
  const Real s = sigma * std::sqrt(Real(1) + y);
  return std::atan(r - s);
}

template<typename Real> 
inline Real 
Spheroid<Real>::conformalLatitude(const Real& u_a, const Real& v_a) const
{
// Equation (128)
// Used on the Backward projection (utmToGeo)
  const Real x = Math::sqr(std::sinh(v_a))  + Math::sqr(std::cos(u_a));
  const Real y = std::sin(u_a) / std::sqrt(x);
  return std::atan(y);
}

template<typename Real>
inline Real Spheroid<Real>::utmCentralMeridian(
  const Real& lonDeg, int &utmZone) const
{
  utmZone = 1 + static_cast<int>((lonDeg + Real(180)) / Real(6));
  const Real lonMer = static_cast<Real>(((utmZone -1) * 6) - 180 + 3);
  return Math::degToRad(lonMer);
}

template<typename Real>
void Spheroid<Real>::geoToUtm(
  const Real& latDeg,
  const Real& lonDeg,
  Real& utmE,
  Real& utmN,
  int& utmZone,
  char& hemisphere) const
{
  const Real phiRad = Math::degToRad(latDeg);
  const Real lmbRad = Math::degToRad(lonDeg);
  const Real cnfLat = conformalLatitude(phiRad);
  const Real lmb0   = utmCentralMeridian(lonDeg, utmZone);
  const Real omega  = lmbRad - lmb0;

  const Real x = Math::sqr(std::tan(cnfLat)) +  Math::sqr(std::cos(omega)); 

  // Gauss-Schreiber coordinates - Equation (84)
  const Real u = majorSemiAxis * std::atan(std::tan(cnfLat) / std::cos(omega));
  const Real v = majorSemiAxis * std::asinh(std::sin(omega) / std::sqrt(x));
  const Real u_a = u / majorSemiAxis;
  const Real v_a = v / majorSemiAxis;

  // Krueger Equations (108)
  Real sumX(0);
  Real sumY(0);
  Real t(2);

  for (int i = 0; i < 8; 
    i++, t += Real(2))
  {
    sumX += alpha[i] * 
      std::cos(t * u_a) * std::sinh(t * v_a);
    sumY += alpha[i] * 
      std::sin(t * u_a) * std::cosh(t * v_a);
  }

  // TM to UTM
  utmE = Real(0.9996) * rectifyingRadius * (v_a + sumX);
  utmN = Real(0.9996) * rectifyingRadius * (u_a + sumY);

  utmE += Real(5E+05);
  utmN += latDeg < Real(0) ? Real(1E+07) : Real(0);
  hemisphere = latDeg < Real(0) ? 'S' : 'N';
}

template<typename Real> 
void Spheroid<Real>::utmToGeo(
  const Real& utmE,
  const Real& utmN,
  const int& utmZone,
  const char& hemisphere,
  Real& latDeg,
  Real& lonDeg) const
{
  Real tmX(utmE);
  Real tmY(utmN);

  // UTM to TM
  tmX -= Real(5E+05);
  tmY -= hemisphere == 'S' ? Real(1E+07) : Real(0);

  const Real X_A = tmX / (Real(0.9996) * rectifyingRadius);
  const Real Y_A = tmY / (Real(0.9996) * rectifyingRadius);

  // Equation (127)
  Real sumX(0);
  Real sumY(0);
  Real t(2);

  for (int i = 0; i < 8; 
    i++, t += Real(2))
  {
    sumX += beta[i] * 
      std::cos(t * Y_A) * std::sinh(t * X_A);
    sumY += beta[i] * 
      std::sin(t * Y_A) * std::cosh(t * X_A);
  }

  const Real v_a = X_A + sumX;
  const Real u_a = Y_A + sumY;

  // Equation (128)
  const Real cnfLat = conformalLatitude(u_a, v_a);
  const Real omega = std::atan(std::sinh(v_a) / std::cos(u_a));

  // Newton-Raphson solution for latitude
  Real aproxT = std::tan(cnfLat);
  Real initialT = aproxT;

  for (int itrCount = 0; itrCount < 25; itrCount++)
  {
    // helpers
    Real oneT = std::sqrt(Real(1) + Math::sqr(aproxT));
    Real oneEccen = Real(1) - eccentricitySqr;

    // Equation (130)
    Real sigma = std::sinh(eccentricity * 
      std::atanh(eccentricity * aproxT / oneT));

    Real oneSigma = std::sqrt(Real(1) + Math::sqr(sigma));

    // Equation (132)
    Real funcAproxT = aproxT * oneSigma - sigma * oneT - initialT;

    // Equation (133)
    Real derivAproxT = (oneSigma * oneT - sigma * aproxT) * 
      (oneEccen * oneT) / (Real(1) + oneEccen * aproxT * aproxT);
    
    // Equation (131)
    Real newAproxT = aproxT - funcAproxT / derivAproxT;

    // Iterating
    Real deltaT = std::fabs(aproxT - newAproxT);
    aproxT = newAproxT;

    if (deltaT < std::numeric_limits<Real>::epsilon())
      break;
  }

  const Real lambda0 = ((utmZone - 1) * 6) -177;

  latDeg = Math::radToDeg(std::atan(aproxT));
  lonDeg = Math::radToDeg(omega) + lambda0;
}

template<typename Real> 
Real Spheroid<Real>::utmScaleFactor(const Real& lonDeg) const
/// abs(err) < 5.923E-07 compared to the full formulation
{
  int utmZone; 
  Real omegaRad = utmCentralMeridian(lonDeg, utmZone) - Math::degToRad(lonDeg);
  return Real(0.5035348161) * Math::sqr(omegaRad) + Real(0.9996);
}

template<typename Real> 
void Spheroid<Real>::computeUtmParams()
{
  // Most part will be evaluated at compile time.
  // Also, this is only computed once at each instatiation.
    
  // Third flattening (up to the 8th power)
  // This is only needed to compute alpha, beta and the rectifying radius
  // Equation (9)
  Real n1 = flattening / (Real(2) - flattening); 
  Real n2 = n1 * n1;    
  Real n3 = n1 * n2;    
  Real n4 = n2 * n2;
  Real n5 = n2 * n3;    
  Real n6 = n3 * n3;    
  Real n7 = n3 * n4;    
  Real n8 = n4 * n4;
  
  // Equation (41)
  rectifyingRadius = (majorSemiAxis / (Real(1) + n1)) * ( 
    Real(1) +
    Real(1) / Real(    4) * n2 + 
    Real(1) / Real(   64) * n4 + 
    Real(1) / Real(  256) * n6 + 
    Real(1) / Real(16384) * n8 
  );

  // Equation (62)
  // Krueger eq, but extended to order n^8

  alpha[0] = // alpha2
    n1 * (Real(       1) / Real(       2)) - 
    n2 * (Real(       2) / Real(       3)) + 
    n3 * (Real(       5) / Real(      16)) + 
    n4 * (Real(      41) / Real(     180)) -
    n5 * (Real(     127) / Real(     288)) +
    n6 * (Real(    7891) / Real(   37800)) +
    n7 * (Real(   72161) / Real(  387072)) -
    n8 * (Real(18975107) / Real(50803200));

  alpha[1] = // alpha4
    n2 * (Real(       13) / Real(       48)) - 
    n3 * (Real(        3) / Real(        5)) + 
    n4 * (Real(      557) / Real(     1440)) +
    n5 * (Real(      281) / Real(      630)) -
    n6 * (Real(  1983433) / Real(  1935360)) +
    n7 * (Real(    13769) / Real(    28800)) +
    n8 * (Real(148003883) / Real(174182400));

  alpha[2] = // alpha6
    n3 * (Real(      61) / Real(     240)) -
    n4 * (Real(     103) / Real(     140)) +
    n5 * (Real(   15061) / Real(   26880)) +
    n6 * (Real(  167603) / Real(  181440)) -
    n7 * (Real(67102379) / Real(29030400)) +
    n8 * (Real(79682431) / Real(79833600));

  alpha[3] = // alpha8
    n4 * (Real(      49561) / Real(    161280)) -
    n5 * (Real(        179) / Real(       168)) +
    n6 * (Real(    6601661) / Real(   7257600)) +
    n7 * (Real(      97445) / Real(     49896)) -
    n8 * (Real(40176129013) / Real(7664025600));

  alpha[4] = // alpha10
    n5 * (Real(     34729) / Real(    80640)) -
    n6 * (Real(   3418889) / Real(  1995840)) +
    n7 * (Real(  14644087) / Real(  9123840)) +
    n8 * (Real(2605413599) / Real(622702080));

  alpha[5] = // alpha12
    n6 * (Real(   212378941) / Real(  319334400)) -
    n7 * (Real(    30705481) / Real(   10378368)) +
    n8 * (Real(175214326799) / Real(58118860800));

  alpha[6] = // alpha14
    n7 * (Real( 1522256789) / Real(1383782400)) -
    n8 * (Real(16759934899) / Real(3113510400));

  alpha[7] = // alpha16
    n8 * (Real(1424729850961) / Real(743921418240));


  // Equation (64)
  // Krueger eq, but extended to order n^8

  beta[0] = - // beta2
    n1 * (Real(      1) / Real(       2)) +
    n2 * (Real(      2) / Real(       3)) -
    n3 * (Real(     37) / Real(      96)) + 
    n4 * (Real(      1) / Real(     360)) +
    n5 * (Real(     81) / Real(     512)) -
    n6 * (Real(  96199) / Real(  604800)) +
    n7 * (Real(5406467) / Real(38707200)) -
    n8 * (Real(7944359) / Real(67737600));

  beta[1] = - // beta4
    n2 * (Real(       1) / Real(       48)) -
    n3 * (Real(       1) / Real(       15)) +
    n4 * (Real(     437) / Real(     1140)) -
    n5 * (Real(      46) / Real(      105)) +
    n6 * (Real( 1118711) / Real(  3870720)) -
    n7 * (Real(   51841) / Real(  1209600)) -
    n8 * (Real(24749483) / Real(348364800));

  beta[2] = - // beta6
    n3 * (Real(     17) / Real(     480)) +
    n4 * (Real(     37) / Real(     840)) +
    n5 * (Real(    209) / Real(    4480)) -
    n6 * (Real(   5569) / Real(   90720)) -
    n7 * (Real(9261899) / Real(58060800)) +
    n8 * (Real(6457463) / Real(17740800));

  beta[3] = - // beta8
    n4 * (Real(     4397) / Real(    161280)) +
    n5 * (Real(       11) / Real(       504)) +
    n6 * (Real(   830251) / Real(   7257600)) -
    n7 * (Real(   466511) / Real(   2494800)) -
    n8 * (Real(324154477) / Real(7664025600));

  beta[4] = - // beta10
    n5 * (Real(    4583) / Real(   161280)) +
    n6 * (Real(  108847) / Real(  3991680)) +
    n7 * (Real( 8005831) / Real( 63866880)) -
    n8 * (Real(22894433) / Real(124540416));

  beta[5] = - // beta12
    n6 * (Real(  20648693) / Real(  638668800)) +
    n7 * (Real(  16363163) / Real(  518918400)) +
    n8 * (Real(2204645983) / Real(12915302400));
    
  beta[6] = - // beta14
    n7 * (Real(219941297) / Real( 5535129600)) +
    n8 * (Real(497323811) / Real(12454041600));

  beta[7] = - // beta16
    n8 * (Real(191773887257) / Real(3719607091200));
}

///////////////////////////////////////////////////////////////////////////////////////
// Inverse Geodesic problem using Thadeus Vincenty formula (1975-76) 
// This may fail for antipodal points (or take thousands of iterations)

template<typename Real>
Real Spheroid<Real>::distance(
  const Real& latDeg1,
  const Real& lonDeg1,
  const Real& latDeg2,
  const Real& lonDeg2)
const {
// TODO: Implement the Karney solution 
  const Real phi1 = Math::degToRad(latDeg1);
  const Real lmb1 = Math::degToRad(lonDeg1);

  const Real phi2 = Math::degToRad(latDeg2);
  const Real lmb2 = Math::degToRad(lonDeg2);

  // reduced latitudes (latitude on the auxiliary sphere)
  const Real U1 = std::atan((Real(1) - flattening)* std::tan(phi1));
  const Real U2 = std::atan((Real(1) - flattening)* std::tan(phi2));
  const Real L = lmb2 - lmb1;

  // Iterative differencein longitude of the points on the auxiliary sphere
  Real lambdaOld = L;
  Real lambda = L;

  Real sigmaSin, sigmaCos, sigma;
  Real twoSigmaCos, twoSigmaCosSqr;
  Real alphaCosSqr;
  Real h0, h1;

  // Usually takes less than 25 iterations
  for (int itrCount = 0; itrCount < 1000; itrCount++)
  {
    h0 = std::cos(U2) * std::sin(lambda);
    h1 = 
      std::cos(U1) * std::sin(U2) - 
      std::sin(U1) * std::cos(U2) * std::cos(lambda); 

    sigmaSin = std::sqrt(h0 * h0 + h1 * h1);
    sigmaCos = 
      std::sin(U1) * std::sin(U2) + 
      std::cos(U1) * std::cos(U2) * std::cos(lambda);

    // angular separation between points
    sigma = std::atan2(sigmaSin, sigmaCos);

    // forward azimuth of the geodesic at the equator (if it were extended that far)
    Real alphaSin = std::cos(U1) * std::cos(U2) * std::sin(lambda) / sigmaSin;
    alphaCosSqr = Real(1) - Math::sqr(alphaSin);
    twoSigmaCos = sigmaCos - (Real(2) * std::sin(U1) * std::sin(U2) / alphaCosSqr);
    twoSigmaCosSqr = Math::sqr(twoSigmaCos);

    Real C = (flattening / Real(16)) * alphaCosSqr * 
      (Real(4) + flattening * (Real(4) - Real(3) * alphaCosSqr));
    
    // new approximation for lambda
    lambdaOld = lambda;
    lambda = L + (Real(1) - C) * flattening * alphaSin * 
      (sigma + C * sigmaSin * (twoSigmaCos + C * sigmaCos * 
      (Real(2) * twoSigmaCosSqr - Real(1))));

    if (std::abs(lambdaOld - lambda) < std::numeric_limits<Real>::epsilon())
      break;
  }

  h0 = majorSemiAxis / minorSemiAxis;
  Real u2 = alphaCosSqr * (h0 * h0 - Real(1));

  // Simpler formulation for A,B given in Vincenty 1976
  Real oneK1 = std::sqrt(Real(1) + u2);
  Real k1 = (oneK1 - Real(1)) / (oneK1 + Real(1));
  Real A = (Real(1) + (k1 * k1 / Real(4))) / (Real(1) - k1);
  Real B = k1 * (Real(1) - (Real(3) / Real(8)) * k1 * k1);

  Real sigmaDelta = B * sigmaSin * (twoSigmaCos + (B / Real(4)) * 
    ((sigmaCos * (Real(2) * twoSigmaCosSqr - Real(1))) - 
    (B / Real(6) * twoSigmaCos * (Real(4) * Math::sqr(sigmaSin) -Real(3)) * 
    (Real(4) * twoSigmaCosSqr -Real(3)))));

  return minorSemiAxis * A * (sigma - sigmaDelta);
}

} /// !namespace geo
#endif /// !_Geodesy_h

