
#include <iostream>
#include "Geodesy.h"

bool testDistance() {
  geo::Spheroid ref;
    
  double PPTE_ecef[3] = {
     3687624.3674,
    -4620818.6827,
    -2386880.3805
  };
  double BRAZ_ecef[3] = {
     4115014.0848,
    -4550641.5491,
    -1741444.0190
  };

  double PPTE_geod[3] = {0};
  double BRAZ_geod[3] = {0};

  ref.ecefToGeod(
    PPTE_ecef[0], PPTE_ecef[1], PPTE_ecef[2],
    PPTE_geod[0], PPTE_geod[1], PPTE_geod[2]
  );

  ref.ecefToGeod(
    BRAZ_ecef[0], BRAZ_ecef[1], BRAZ_ecef[2],
    BRAZ_geod[0], BRAZ_geod[1], BRAZ_geod[2]
  );

  double distMeters = ref.distance(
    PPTE_geod[0], PPTE_geod[1], 
    BRAZ_geod[0], BRAZ_geod[1]);

  return std::abs( distMeters - 777678.312 ) < 0.001;
}

bool testUtm() {
  geo::Spheroid ref;
    
  double PPTE_ecef_1[3] = {
     3687624.3674,
    -4620818.6827,
    -2386880.3805
  };
  double PPTE_ecef_2[3] = {0};

  double PPTE_geod_1[3] = {0};
  double PPTE_geod_2[3] = {0};
  double PPTE_utm[2] = {0};
  int utmZone;
  char hemisphere;

  ref.ecefToGeod(
    PPTE_ecef_1[0], PPTE_ecef_1[1], PPTE_ecef_1[2],
    PPTE_geod_1[0], PPTE_geod_1[1], PPTE_geod_1[2]
  );

  ref.geoToUtm(
    PPTE_geod_1[0], PPTE_geod_1[1],
    PPTE_utm[0], PPTE_utm[1], utmZone, hemisphere
  );

  ref.utmToGeo(
    PPTE_utm[0], PPTE_utm[1], utmZone, hemisphere,     
    PPTE_geod_2[0], PPTE_geod_2[1]
  );

  // we need to copy the height (not used in the UTM method)
  PPTE_geod_2[2] = PPTE_geod_1[2];

  ref.geodToEcef(
    PPTE_geod_2[0], PPTE_geod_2[1], PPTE_geod_2[2],
    PPTE_ecef_2[0], PPTE_ecef_2[1], PPTE_ecef_2[2]  
  );

  return 
    std::abs( PPTE_ecef_2[0] - PPTE_ecef_1[0] ) < 0.001 &&
    std::abs( PPTE_ecef_2[1] - PPTE_ecef_1[1] ) < 0.001 &&
    std::abs( PPTE_ecef_2[2] - PPTE_ecef_1[2] ) < 0.001;
}

template<class TestCallBack>
void runTest(TestCallBack callBack, const char *testName)
{
  if (callBack()) 
    std::cout << "[  OK  ] " << testName << std::endl;
  else 
    std::cout << "[ FAIL ] " << testName << std::endl;

}

int main(int argc, char **argv)
{
  runTest( testDistance, "Geodetic distance" );
  runTest( testUtm,      "UTM projection"    );

  return 0;
}

