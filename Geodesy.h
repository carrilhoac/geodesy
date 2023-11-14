
//  Copyright (c) 2023, Andre Caceres Carrilho
//
//  Geodesy.h common geodesic routines 
//  tags header-only, thread-safe, portable, standard-cpp
//    utm, universal transverse mercator, geodesy, geodesic
//

#ifndef _Geodesy_h
#define _Geodesy_h 1

#include <cmath>
#include <limits>

namespace geo {

// Templated because you might want to use long double.
// The term 'Spheroid' is what Open Geospatial Consortium uses.
template<typename FloatType = double> 
class Spheroid{
private:
  FloatType majorSemiAxis; // equatorial radius (m)
  FloatType minorSemiAxis; // polar semi axis (m)
  FloatType flattening;
  FloatType eccentricitySqr;
  FloatType eccentricity;

  // Coefficients for Karney-Kruger UTM equations: 
  // Deakin, R.E.; Hunter, M.N.; Karney, C.F.F. A fresh look at the UTM projection.
  FloatType rectifyingRadius;
  FloatType alpha[8], beta[8];

public:
  Spheroid(             
    FloatType majorSemiAxis_ = FloatType(6378137), // defaults to WGS-84 
    FloatType flattening_    = FloatType(1) / FloatType(298.257223563)) 
  {
    setFrom(majorSemiAxis_, flattening_);
  }

  void setFrom(FloatType majorSemiAxis_, FloatType flattening_);

  void setWgs84(); // GeoJSON only supports EPSG:4326 
  void setGrs80(); // SIRGAS used on EPSG:4674

  // Conversion between Geodetic and Geocentric Cartesian coordinates,
  // which is also known as Earth Centered, Earth Fixed (ECEF)
  void geodToEcef(
    const FloatType& latDeg,
    const FloatType& lonDeg,
    const FloatType& hGeo, 
    FloatType& ecefX, 
    FloatType& ecefY, 
    FloatType& ecefZ) const;

  void ecefToGeod(
    const FloatType& ecefX, 
    const FloatType& ecefY, 
    const FloatType& ecefZ, 
    FloatType& latDeg, 
    FloatType& lonDeg, 
    FloatType& hGeo) const;

  // UTM projection (forward and backward)
  // Uses Karney-Kruger equations (less than 5 nm error)
  void geoToUtm(
    const FloatType& latDeg,
    const FloatType& lonDeg,
    FloatType& utmE,
    FloatType& utmN,
    int& utmZone,
    char& hemisphere) const;

  void utmToGeo(
    const FloatType& utmE,
    const FloatType& utmN,
    const int& utmZone,
    const char& hemisphere,
    FloatType& latDeg,
    FloatType& lonDeg) const;

  // Geodesic Inverse problem 
  // Vincenty equations (less than 0.5 mm error) 
  // Note: might fail for near antipodal points
  FloatType distance(
    const FloatType& latDeg1, 
    const FloatType& lonDeg1, 
    const FloatType& latDeg2, 
    const FloatType& lonDeg2
  ) const;

private:
  inline FloatType conformalLatitude(const FloatType& phiRad) const;
  inline FloatType conformalLatitude(const FloatType& u_a, const FloatType& v_a) const;
  inline FloatType radiusDenom(const FloatType& phiRad) const;
  inline FloatType radiusPrime(const FloatType& phiRad) const;
  inline FloatType radiusMeridional(const FloatType& phiRad) const;
  inline FloatType utmCentralMeridian(const FloatType& lonDeg, int &utmZone) const;
  void computeUtmParams();
};

namespace Math 
{
  template<typename FloatType>
  inline FloatType degToRad(const FloatType& x){
    return x * FloatType(3.1415926535897932384626) / FloatType(180);
  }
  template<typename FloatType>
  inline FloatType radToDeg(const FloatType& x){
    return x * FloatType(180) / FloatType(3.1415926535897932384626);
  }
}

template<typename FloatType> 
void Spheroid<FloatType>::setWgs84(){
  setFrom(FloatType(6378137), FloatType(1) / FloatType(298.257223563));
}

template<typename FloatType> 
void Spheroid<FloatType>::setGrs80(){
  setFrom(FloatType(6378137), FloatType(1) / FloatType(298.2572221009));
}

template<typename FloatType> 
void Spheroid<FloatType>::setFrom(FloatType majorSemiAxis_, FloatType flattening_){
  majorSemiAxis = majorSemiAxis_;
  flattening    = flattening_;

  eccentricitySqr  = flattening * (FloatType(2) - flattening);
  eccentricity     = std::sqrt(eccentricitySqr);
  minorSemiAxis    = majorSemiAxis * std::sqrt(FloatType(1) - eccentricitySqr);
  computeUtmParams();
}

template<typename FloatType> 
void Spheroid<FloatType>::geodToEcef(
  const FloatType& latDeg,
  const FloatType& lonDeg,
  const FloatType& hGeo, 
  FloatType& ecefX, 
  FloatType& ecefY, 
  FloatType& ecefZ
) const{
  const FloatType phiRad = Math::degToRad(latDeg);
  const FloatType lmbRad = Math::degToRad(lonDeg);
  const FloatType radiusN = radiusPrime(phiRad);
  ecefX = (radiusN + hGeo) * std::cos(phiRad) * std::cos(lmbRad);
  ecefY = (radiusN + hGeo) * std::cos(phiRad) * std::sin(lmbRad);
  ecefZ = (radiusN * (FloatType(1) - eccentricitySqr) + hGeo) * std::sin(phiRad);
}

template<typename FloatType> 
void Spheroid<FloatType>::ecefToGeod(
  const FloatType& ecefX, 
  const FloatType& ecefY, 
  const FloatType& ecefZ, 
  FloatType& latDeg, 
  FloatType& lonDeg, 
  FloatType& hGeo
) const{
  const FloatType rho = std::hypot(ecefX, ecefY);
  const FloatType eps = std::numeric_limits<FloatType>::epsilon();

  FloatType radiusN = majorSemiAxis;
  FloatType phiSin;
  FloatType zi = ecefZ;
  FloatType zk = FloatType(0);

  for (int itrCount = 0; 
    // defining an upper bound for the loop 
    itrCount++ < 25 && std::fabs(zi- zk) > eps; )
  {
    zk = zi;
    phiSin = zi / std::hypot(rho, zi);
    radiusN = majorSemiAxis / 
      std::sqrt(FloatType(1) - eccentricitySqr * phiSin * phiSin);
    zi = ecefZ + radiusN * eccentricitySqr * phiSin;
  }

  FloatType lonRad = FloatType(0);
  FloatType latRad = FloatType(0);

  if (rho > eps){
    latRad = std::atan(zi / rho);
    lonRad = std::atan2(ecefY, ecefX);
  }
  else{
    if (ecefZ >= FloatType(0))
      latRad = FloatType( 0.5) * FloatType(3.1415926535897932384626);
    else 
      latRad = FloatType(-0.5) * FloatType(3.1415926535897932384626);
  }

  latDeg = Math::radToDeg(latRad);
  lonDeg = Math::radToDeg(lonRad);
  hGeo = std::hypot(rho, zi) - radiusN;
}

template<typename FloatType> 
inline FloatType Spheroid<FloatType>::radiusDenom(const FloatType& phiRad) const{
  return FloatType(1) - eccentricitySqr * std::sin(phiRad) * std::sin(phiRad);
}
template<typename FloatType> 
inline FloatType Spheroid<FloatType>::radiusPrime(const FloatType& phiRad) const{
  return majorSemiAxis / std::sqrt( radiusDenom(phiRad) );
}
template<typename FloatType> 
inline FloatType Spheroid<FloatType>::radiusMeridional(const FloatType& phiRad) const{
  return (majorSemiAxis * (FloatType(1) - eccentricitySqr)) / 
    std::pow(radiusDenom(phiRad) , FloatType(1.5));
}  

///////////////////////////////////////////////////////////////////////////////////////
// Karney-Kruger UTM equations from: 
// Deakin, R.E.; Hunter, M.N.; Karney, C.F.F. A fresh look at the UTM projection.
// Very accurate, less than 5 nanometer error

template<typename FloatType> 
inline FloatType 
Spheroid<FloatType>::conformalLatitude(const FloatType& latRad) const
{
// Equations (88) and (89)
// Used on the Forward projection (geoToUtm)
  const FloatType x = std::tan(latRad);
  const FloatType y = x * x;
  const FloatType z = (eccentricity * x) / std::sqrt(FloatType(1) + y);
  const FloatType sigma = std::sinh(eccentricity * std::atanh(z));
  const FloatType r = x * std::sqrt(FloatType(1) + sigma * sigma);
  const FloatType s = sigma * std::sqrt(FloatType(1) + y);
  return std::atan(r - s);
}

template<typename FloatType> 
inline FloatType 
Spheroid<FloatType>::conformalLatitude(const FloatType& u_a, const FloatType& v_a) const
{
// Equation (128)
// Used on the Backward projection (utmToGeo)
  const FloatType x = std::sinh(v_a) * std::sinh(v_a) + std::cos(u_a) * std::cos(u_a);
  const FloatType y = std::sin(u_a) / std::sqrt(x);
  return std::atan(y);
}

template<typename FloatType>
inline FloatType Spheroid<FloatType>::utmCentralMeridian(
  const FloatType& lonDeg, int &utmZone) const
{
  utmZone = 1 + static_cast<int>((lonDeg + FloatType(180)) / FloatType(6));
  const FloatType lonMer = static_cast<FloatType>(((utmZone -1) * 6) - 180 + 3);
  return Math::degToRad(lonMer);
}

template<typename FloatType>
void Spheroid<FloatType>::geoToUtm(
  const FloatType& latDeg,
  const FloatType& lonDeg,
  FloatType& utmE,
  FloatType& utmN,
  int& utmZone,
  char& hemisphere) const
{
  const FloatType phiRad = Math::degToRad(latDeg);
  const FloatType lmbRad = Math::degToRad(lonDeg);
  const FloatType cnfLat = conformalLatitude(phiRad);
  const FloatType lmb0   = utmCentralMeridian(lonDeg, utmZone);
  const FloatType omega  = lmbRad - lmb0;

  const FloatType x = 
    std::tan(cnfLat) * std::tan(cnfLat) + 
    std::cos(omega)  * std::cos(omega); 

  // Gauss-Schreiber coordinates - Equation (84)
  const FloatType u = majorSemiAxis * std::atan(std::tan(cnfLat) / std::cos(omega));
  const FloatType v = majorSemiAxis * std::asinh(std::sin(omega) / std::sqrt(x));
  const FloatType u_a = u / majorSemiAxis;
  const FloatType v_a = v / majorSemiAxis;

  // Krueger Equations (108)
  FloatType sumX(0);
  FloatType sumY(0);
  FloatType t(2);

  for (int i = 0; i < 8; 
    i++, t += FloatType(2))
  {
    sumX += alpha[i] * 
      std::cos(t * u_a) * std::sinh(t * v_a);
    sumY += alpha[i] * 
      std::sin(t * u_a) * std::cosh(t * v_a);
  }

  // TM to UTM
  utmE = FloatType(0.9996) * rectifyingRadius * (v_a + sumX);
  utmN = FloatType(0.9996) * rectifyingRadius * (u_a + sumY);

  utmE += FloatType(5E+05);
  utmN += latDeg < FloatType(0) ? FloatType(1E+07) : FloatType(0);
  hemisphere = latDeg < FloatType(0) ? 'S' : 'N';
}

template<typename FloatType> 
void Spheroid<FloatType>::utmToGeo(
  const FloatType& utmE,
  const FloatType& utmN,
  const int& utmZone,
  const char& hemisphere,
  FloatType& latDeg,
  FloatType& lonDeg) const
{
  FloatType tmX(utmE);
  FloatType tmY(utmN);

  // UTM to TM
  tmX -= FloatType(5E+05);
  tmY -= hemisphere == 'S' ? FloatType(1E+07) : FloatType(0);

  const FloatType X_A = tmX / (FloatType(0.9996) * rectifyingRadius);
  const FloatType Y_A = tmY / (FloatType(0.9996) * rectifyingRadius);

  // Equation (127)
  FloatType sumX(0);
  FloatType sumY(0);
  FloatType t(2);

  for (int i = 0; i < 8; 
    i++, t += FloatType(2))
  {
    sumX += beta[i] * 
      std::cos(t * Y_A) * std::sinh(t * X_A);
    sumY += beta[i] * 
      std::sin(t * Y_A) * std::cosh(t * X_A);
  }

  const FloatType v_a = X_A + sumX;
  const FloatType u_a = Y_A + sumY;

  // Equation (128)
  const FloatType cnfLat = conformalLatitude(u_a, v_a);
  const FloatType omega = std::atan(std::sinh(v_a) / std::cos(u_a));

  // Newton-Raphson solution for latitude
  FloatType aproxT = std::tan(cnfLat);
  FloatType initialT = aproxT;

  for (int i = 0; i < 25; i++)
  {
    // helpers
    FloatType oneT = std::sqrt(FloatType(1) + aproxT * aproxT);
    FloatType oneEccen = FloatType(1) - eccentricity * eccentricity;

    // Equation (130)
    FloatType sigma = std::sinh(eccentricity * 
      std::atanh(eccentricity * aproxT / oneT));

    FloatType oneSigma = std::sqrt(FloatType(1) + sigma * sigma);

    // Equation (132)
    FloatType funcAproxT = aproxT * oneSigma - sigma * oneT - initialT;

    // Equation (133)
    FloatType derivAproxT = (oneSigma * oneT - sigma * aproxT) * 
      (oneEccen * oneT) / (FloatType(1) + oneEccen * aproxT * aproxT);
    
    // Equation (131)
    FloatType newAproxT = aproxT - funcAproxT / derivAproxT;

    // Iterating
    FloatType deltaT = std::fabs(aproxT - newAproxT);
    aproxT = newAproxT;

    if (deltaT < std::numeric_limits<FloatType>::epsilon())
      break;
  }

  const FloatType lambda0 = ((utmZone - 1) * 6) -177;

  latDeg = Math::radToDeg(std::atan(aproxT));
  lonDeg = Math::radToDeg(omega) + lambda0;
}

template<typename FloatType> 
void Spheroid<FloatType>::computeUtmParams()
{
// Most part will be evaluated at compile time.
// Also, this is only computed once at each instatiation.

  using FP = FloatType;

  // Third flattening (up to the 8th power)
  // This is only needed to compute alpha, beta and the rectifying radius
  
  // Equation (9)
  FP n1 = flattening / (FP(2) - flattening); 
  FP n2 = n1 * n1;    
  FP n3 = n1 * n2;    
  FP n4 = n2 * n2;
  FP n5 = n2 * n3;    
  FP n6 = n3 * n3;    
  FP n7 = n3 * n4;    
  FP n8 = n4 * n4;
  
  // Equation (41)
  rectifyingRadius = (majorSemiAxis / (FP(1) + n1)) * ( 
    FP(1) +
    FP(1) / FP(    4) * n2 + 
    FP(1) / FP(   64) * n4 + 
    FP(1) / FP(  256) * n6 + 
    FP(1) / FP(16384) * n8 
  );

  // Equation (62)
  // Krueger eq, but extended to order n^8

  alpha[0] = // alpha2
    n1 * (FP(       1) / FP(       2)) - 
    n2 * (FP(       2) / FP(       3)) + 
    n3 * (FP(       5) / FP(      16)) + 
    n4 * (FP(      41) / FP(     180)) -
    n5 * (FP(     127) / FP(     288)) +
    n6 * (FP(    7891) / FP(   37800)) +
    n7 * (FP(   72161) / FP(  387072)) -
    n8 * (FP(18975107) / FP(50803200));

  alpha[1] = // alpha4
    n2 * (FP(       13) / FP(       48)) - 
    n3 * (FP(        3) / FP(        5)) + 
    n4 * (FP(      557) / FP(     1440)) +
    n5 * (FP(      281) / FP(      630)) -
    n6 * (FP(  1983433) / FP(  1935360)) +
    n7 * (FP(    13769) / FP(    28800)) +
    n8 * (FP(148003883) / FP(174182400));

  alpha[2] = // alpha6
    n3 * (FP(      61) / FP(     240)) -
    n4 * (FP(     103) / FP(     140)) +
    n5 * (FP(   15061) / FP(   26880)) +
    n6 * (FP(  167603) / FP(  181440)) -
    n7 * (FP(67102379) / FP(29030400)) +
    n8 * (FP(79682431) / FP(79833600));

  alpha[3] = // alpha8
    n4 * (FP(      49561) / FP(    161280)) -
    n5 * (FP(        179) / FP(       168)) +
    n6 * (FP(    6601661) / FP(   7257600)) +
    n7 * (FP(      97445) / FP(     49896)) -
    n8 * (FP(40176129013) / FP(7664025600));

  alpha[4] = // alpha10
    n5 * (FP(     34729) / FP(    80640)) -
    n6 * (FP(   3418889) / FP(  1995840)) +
    n7 * (FP(  14644087) / FP(  9123840)) +
    n8 * (FP(2605413599) / FP(622702080));

  alpha[5] = // alpha12
    n6 * (FP(   212378941) / FP(  319334400)) -
    n7 * (FP(    30705481) / FP(   10378368)) +
    n8 * (FP(175214326799) / FP(58118860800));

  alpha[6] = // alpha14
    n7 * (FP( 1522256789) / FP(1383782400)) -
    n8 * (FP(16759934899) / FP(3113510400));

  alpha[7] = // alpha16
    n8 * (FP(1424729850961) / FP(743921418240));


  // Equation (64)
  // Krueger eq, but extended to order n^8

  beta[0] = - // beta2
    n1 * (FP(      1) / FP(       2)) +
    n2 * (FP(      2) / FP(       3)) -
    n3 * (FP(     37) / FP(      96)) + 
    n4 * (FP(      1) / FP(     360)) +
    n5 * (FP(     81) / FP(     512)) -
    n6 * (FP(  96199) / FP(  604800)) +
    n7 * (FP(5406467) / FP(38707200)) -
    n8 * (FP(7944359) / FP(67737600));

  beta[1] = - // beta4
    n2 * (FP(       1) / FP(       48)) -
    n3 * (FP(       1) / FP(       15)) +
    n4 * (FP(     437) / FP(     1140)) -
    n5 * (FP(      46) / FP(      105)) +
    n6 * (FP( 1118711) / FP(  3870720)) -
    n7 * (FP(   51841) / FP(  1209600)) -
    n8 * (FP(24749483) / FP(348364800));

  beta[2] = - // beta6
    n3 * (FP(     17) / FP(     480)) +
    n4 * (FP(     37) / FP(     840)) +
    n5 * (FP(    209) / FP(    4480)) -
    n6 * (FP(   5569) / FP(   90720)) -
    n7 * (FP(9261899) / FP(58060800)) +
    n8 * (FP(6457463) / FP(17740800));

  beta[3] = - // beta8
    n4 * (FP(     4397) / FP(    161280)) +
    n5 * (FP(       11) / FP(       504)) +
    n6 * (FP(   830251) / FP(   7257600)) -
    n7 * (FP(   466511) / FP(   2494800)) -
    n8 * (FP(324154477) / FP(7664025600));

  beta[4] = - // beta10
    n5 * (FP(    4583) / FP(   161280)) +
    n6 * (FP(  108847) / FP(  3991680)) +
    n7 * (FP( 8005831) / FP( 63866880)) -
    n8 * (FP(22894433) / FP(124540416));

  beta[5] = - // beta12
    n6 * (FP(  20648693) / FP(  638668800)) +
    n7 * (FP(  16363163) / FP(  518918400)) +
    n8 * (FP(2204645983) / FP(12915302400));
    
  beta[6] = - // beta14
    n7 * (FP(219941297) / FP( 5535129600)) +
    n8 * (FP(497323811) / FP(12454041600));

  beta[7] = - // beta16
    n8 * (FP(191773887257) / FP(3719607091200));
}

///////////////////////////////////////////////////////////////////////////////////////
// Inverse Geodesic problem using Thadeus Vincenty formula (1975-76) 
// This may fail for antipodal points (or take thousands of iterations)

template<typename FloatType>
FloatType Spheroid<FloatType>::distance(
  const FloatType& latDeg1,
  const FloatType& lonDeg1,
  const FloatType& latDeg2,
  const FloatType& lonDeg2)
const {
// TODO: Implement the Karney solution 
  const FloatType phi1 = Math::degToRad(latDeg1);
  const FloatType lmb1 = Math::degToRad(lonDeg1);

  const FloatType phi2 = Math::degToRad(latDeg2);
  const FloatType lmb2 = Math::degToRad(lonDeg2);

  // reduced latitudes (latitude on the auxiliary sphere)
  const FloatType U1 = std::atan((FloatType(1) - flattening)* std::tan(phi1));
  const FloatType U2 = std::atan((FloatType(1) - flattening)* std::tan(phi2));
  const FloatType L = lmb2 - lmb1;

  // Iterative differencein longitude of the points on the auxiliary sphere
  FloatType lambdaOld = L;
  FloatType lambda = L;

  FloatType sigmaSin, sigmaCos, sigma;
  FloatType twoSigmaCos, twoSigmaCosSqr;
  FloatType alphaCosSqr;
  FloatType h0, h1;

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
    FloatType alphaSin = std::cos(U1) * std::cos(U2) * std::sin(lambda) / sigmaSin;
    alphaCosSqr = FloatType(1) - alphaSin * alphaSin;
    twoSigmaCos = sigmaCos - 
      (FloatType(2) * std::sin(U1) * std::sin(U2) / alphaCosSqr);
    twoSigmaCosSqr = twoSigmaCos * twoSigmaCos ;

    h0 = FloatType(4) + flattening * (FloatType(4) - FloatType(3) * alphaCosSqr);
    FloatType C = (flattening / FloatType(16)) * alphaCosSqr * h0;
    
    h1 = FloatType(2) * twoSigmaCosSqr - FloatType(1); 
    h1 = sigma + C * sigmaSin * (twoSigmaCos + C * sigmaCos * h1);

    // new approximation for lambda
    lambdaOld = lambda;
    lambda = L + (FloatType(1) - C) * flattening * alphaSin * h1;

    if (std::abs(lambdaOld - lambda) < std::numeric_limits<FloatType>::epsilon())
      break;
  }

  h0 = majorSemiAxis / minorSemiAxis;
  FloatType u2 = alphaCosSqr * (h0 * h0 - FloatType(1));

  #if 0
  // Original formulation Vincenty 1975
  FloatType A = FloatType(1) + (u2 / FloatType(16384)) * 
    (FloatType(4096) + u2 * (u2 * (FloatType(320) - FloatType(175) * u2) - FloatType(768)));
  FloatType B = u2 / FloatType(1024) * 
    (FloatType( 256) + u2 * (u2 * (FloatType( 74) - FloatType( 47) * u2) - FloatType(128)));
  #else 
  // Simpler formulation from Vincenty 1976
  FloatType oneK1 = std::sqrt(FloatType(1) + u2);
  FloatType k1 = (oneK1 - FloatType(1)) / (oneK1 + FloatType(1));
  FloatType A = (FloatType(1) + (k1 * k1 / FloatType(4))) / (FloatType(1) - k1);
  FloatType B = k1 * (FloatType(1) - (FloatType(3) / FloatType(8)) * k1 * k1);
  #endif

  h0 = sigmaCos * (FloatType(2) * twoSigmaCosSqr - FloatType(1));
  h1 = B / FloatType(6) * twoSigmaCos * 
    (FloatType(4) * sigmaSin * sigmaSin -FloatType(3)) * 
    (FloatType(4) * twoSigmaCosSqr -FloatType(3));

  FloatType sigmaDelta = B * sigmaSin * (twoSigmaCos + (B / FloatType(4)) * (h0 - h1));
  return minorSemiAxis * A * (sigma - sigmaDelta);
}

} // ns geo

#endif 

