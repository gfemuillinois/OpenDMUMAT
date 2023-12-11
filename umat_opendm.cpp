// Implementation of OpenDM as a UMAT
// by: Bryce Mazurowski (brycepm2@gmail.com)
//

#include <iostream>
#include <string>
#include <ctime>
#include <fstream>
#include <math.h>

using std::log;
using std::cout;
using std::endl;

// declarations 
/** @brief Runs simple matPt test on strain controlled load
 */
int runMatPtTest(const int modelVal);

/** @brief Runs TaylorCheck (tm Nathan Shauer) to test
 * material tangent matrix
 */
double taylorCheck(const int modelVal);

/** @brief Calculate L2 Norm of vector
 * sqrt(v_i^2)
 */
double l2Norm(const double *v1);

/** @brief Matrix Vector Product
 * A_{ij} b_{j}
 */
void matVecProd(const double* A, const double* b,
                const double alpha, double* c);

void printVec(const double* v);
void printMat(const double* A);

/** @brief run umat to estimate numerical tangent
 */
double numTangCalc(const int modelVal);

/** @brief Estimates material tangent via perturbed stress
 */
void numericalTangent(const double* stress0, const double* stressPert,
                      const double strainInc, double* numTang);


// /********************************************************************/
// /********************************************************************/

extern "C" {
void runInitOpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
                        double* strain, double* dstrain, double* stress,
                        double* ddsdde, double* sse, double* spd, double* scd);
}

/********************************************************************/
/********************************************************************/

extern "C" void umat(double *stress, double *statev, double *ddsdde,
                     double *sse, double *spd, 
                     double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
                     double *stran, double *dstran, double *time, double *dtime, double *temp,
                     double *dtemp, double *predef, double *dpred, char *cmname, int *ndi,
                     int *nshr, int *ntens, int *nstatv, double *props, int *nprops,
                     double *coords, double *drot, double *pnewdt, double *celent, double *dfgrd0,
                     double *dfgrd1, int *noel, int *npt, int *layer, int *kspt,
                     int *kstep, int *kinc) {
  // SetUp and run openDM model from UMAT call
  runInitOpenDMModel(props, nprops, statev, nstatv, stran, dstran, stress, ddsdde, sse, spd, scd);
}
/********************************************************************/
/********************************************************************/

int main(int argc, char* argv[]) {
  // model choice: 0 - default(2), 2 - 2Param, 4 - 4Param
  int modelVal = 0;
  if ( argc != 1 ) {
    modelVal = std::stoi(argv[1]);
  } else {
    modelVal = 4;
  }
  std::cout << "Using model: " << modelVal << std::endl;
  if (true) {
    const int test = runMatPtTest(modelVal);
    return test;
  } else {
    const double convVal = taylorCheck(modelVal);
  }
  return 0;
}
/********************************************************************/
/********************************************************************/

int runMatPtTest(const int modelVal) {
  // elastic material definitions
  double E11, E22, E33, nu12, nu23, nu13, G12, G23, G13;
  // damage model parameters
  double hs1_11, hs1_12, hs1_13,
    hs2_22, hs2_12, hs2_23,
    hs4_11, hs4_12, hs4_13, hs4_16,
    hs5_11, hs5_12, hs5_13, hs5_16,
    b1, b2, b6,
    y01, y02, y04, y05,
    yc1, yc2, yc4, yc5,
    pe1, pe2, pe4, pe5,
    dc1, dc2, dc4, dc5;
  // applied load
  double loadApplied;
  double maxStrain;
  int nInc = 100;
  const bool pythonTest = true;
  if ( !pythonTest ) {
    // set elasticity
    E11 = 12.416e6; E22 = E11; E33 = 7.513e6;
    nu12 = 0.1106; nu23 = 0.1707; nu13 = nu23;
    G12 = 3.332e6; G23 = 0.670e6; G13 = G23;
    // damage model parameters
    hs1_11 = 2.5; hs1_12 = hs1_11; hs1_13 = hs1_11;
    hs2_22 = 2.5; hs2_12 = hs2_22; hs2_23 = hs2_22;
    hs4_11 = 2.5; hs4_12 = hs4_11; hs4_13 = hs4_11; hs4_16 = hs4_11;
    hs5_11 = 2.5; hs5_12 = hs5_11; hs5_13 = hs5_11; hs5_16 = hs4_11;
    b1 = 0.2; b2 = b1; b6 = b1;
    y01 = 0.005; y02 = y01; y04 = y01; y05 = y01;
    yc1 = 500; yc2 = yc1; yc4 = yc1; yc5 = yc1;
    pe1 = 1.5; pe2 = pe1; pe4 = pe1; pe5 = pe1;
    dc1 = 2.0; dc2 = dc1; dc4 = dc1; dc5 = dc1;
    // load
    maxStrain = 0.003;
  } else {
    // Elastic Properties
    E11 = 335000;       //Young Modulus [MPa]
    E22 = 160000;       //Young Modulus [MPa]
    E33 = 335000;       //Young Modulus [MPa]
    nu12 = 0.17;         //Poisson ratio
    nu23 = 0.00;         //Poisson ratio
    nu13 = 0.00;         //Poisson ratio
    G12 = 95360;        //Shear modulus [MPa]
    G23 = 95360;        //Shear modulus [MPa]
    G13 = 95360;        //Shear modulus [MPa]
    
    // Material State Variables
    hs1_11 = 1.2;        //Damage Effect
    hs1_12 = 1.0;        //Damage Effect
    hs1_13 = 1.0;        //Damage Effect
    hs2_22 = 20.0;       //Damage Effect
    hs2_12 = 1.0;        //Damage Effect
    hs2_23 = 1.0;        //Damage Effect
    hs4_11 = 2.0;        //Damage Effect
    hs4_12 = 0.5;        //Damage Effect
    hs4_13 = 0.5;        //Damage Effect
    hs4_16 = 1.0;        //Damage Effect
    hs5_11 = 2.0;        //Damage Effect
    hs5_12 = 0.5;        //Damage Effect
    hs5_13 = 0.5;        //Damage Effect
    hs5_16 = 1.0;        //Damage Effect
    b1 = 0.3;           //Traction / Shear Coupling 1
    b2 = 0.3;           //Traction / Shear Coupling 2
    b6 = 0.3;           //Traction / Shear Coupling 2
    y01 = 0.3;         //Damage Thresholds 1
    y02 = 0.04;         //Damage Thresholds 2
    y04 = 0.04;         //Damage Thresholds 1
    y05 = 0.04;         //Damage Thresholds 2
    yc1 = 16;           //Damage Evolution Celerity 1
    yc2 = 0.22;         //Damage Evolution Celerity 2
    yc4 = 10;           //Damage Evolution Celerity 1
    yc5 = 10;          //Damage Evolution Celerity 2
    pe1 = 1.15;         //Damage Evolution Exponents 1
    pe2 = 2.4;          //Damage Evolution Exponents 2
    pe4 = 1.1;        //Damage Evolution Exponents 1
    pe5 = 1.1;          //Damage Evolution Exponents 2
    dc1 = 5.5;          //Damage Saturations 1
    dc2 = 1.0;         //Damage Saturations 2
    dc4 = 8;          //Damage Saturations 1
    dc5 = 8;         //Damage Saturations 2
    maxStrain = 0.01;
  }

  // need to set a bunch of stuff
  double strain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double stress[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsdde[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  char cmname[80] = "OpenDM";
  int nprops, nstatv;
  // Set size to be largest possible
  double props[42];
  // elastic props
  props[0] = E11; props[1] = E22; props[2] = E33;
  props[3] = nu12; props[4] = nu23; props[5] = nu13;
  props[6] = G12; props[7] = G23; props[8] = G13;
  if (modelVal == 2) {
    cmname[6] = '2';
    nprops = 21; // Number of properties
    nstatv = 13; // number of statVars
    double props2Param[] = {hs1_11, hs2_22, b1,  b2, y01, y02,
                            yc1, yc2, pe1, pe2, dc1, dc2};
    for (int iProp = 0; iProp < 12; iProp++) {
      props[9+iProp] = props2Param[iProp];
    }
  } else if (modelVal == 4) {
    cmname[6] = '4';
    nprops = 42; // Number of properties
    nstatv = 17; // Number stateVars
    double props4Param[] = {hs1_11, hs2_12, hs1_13,
                            hs2_22, hs2_12, hs2_23,
                            hs4_11, hs4_12, hs4_13, hs4_16,
                            hs5_11, hs5_12, hs5_13, hs5_16,
                            b1, b2, b6,
                            y01, y02, y04, y05,
                            yc1, yc2, yc4, yc5,
                            pe1, pe2, pe4, pe5,
                            dc1, dc2, dc4, dc5};
    for (int iProp = 0; iProp < 33; iProp++) {
      props[9+iProp] = props4Param[iProp];
    }
  }
  double sse, spd, scd, rpl;
  // Current deltaT for which solution is seeking
  double dTime = 0.001;
  double time[] = {0.0, 0.001};

  int ndi = 3;            // 3 direct stress/strain components
  int nshr = 3;           // 3 shear stress/strain components
  int ntens = ndi + nshr; // 6 total componenets. ndi+nshr
  int noel = 1;
  int npt = 0; // Integration point number
  double ddsddt = 0.;
  double drplde = 0.;
  double drpldt = 0.;
  double temp = 0.;
  double dTemp = 0.;
  double predef = 0.;
  double dpred = 0.;
  double coords = 0.;
  double drot = 0.;
  double pnewdt = 0.;
  double dfgrd0 = 0.;
  double dfgrd1 = 0.;
  double celent;
  int layer, kspt, kstep, kinc;
  layer = kspt = kstep = kinc = 0;

  // setup a loop to call umat each step
  double stepSize = maxStrain / (nInc);
  std::ofstream myFile;
  std::string modelNum = std::to_string(modelVal);
  for (int runInt = 0; runInt < 4; ++runInt) {
    cout << "Running " << runInt << endl;
    // file for each run
    std::string runStr = std::to_string(runInt);
    std::string fileName = "umatTestOut_"+modelNum+"_Run"+runStr+".csv";
    std::string sep = ", ";
    myFile.open(fileName);
    myFile << "step, "
           << "strain11, strain22, strain12, "
           << "stress11, stress22, stress12, "
           << "d1, d2, d4, d5" << std::endl;
    // zero out SVs for each new load
    double statev[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int iPull = 0; iPull <= nInc; iPull++) {
      // Increment dstrain
      switch(runInt) {
      case 0:
        // Aligned Tension
        dstrain[0] = iPull * stepSize;
        dstrain[1] = 0.0;
        dstrain[3] = 0.0;
        break;
      case 1:
        // 90 degree tension
        dstrain[0] = 0.0;
        dstrain[1] = iPull*stepSize;
        dstrain[3] = 0.0;
        break;
      case 2:
        // +/- 45 tension
        dstrain[0] = iPull*0.5*stepSize;
        dstrain[1] = iPull*0.5*stepSize;
        dstrain[3] = iPull*stepSize;
        break;
      case 3:
        // +/- 45 tension
        dstrain[0] = iPull*0.5*stepSize;
        dstrain[1] = iPull*0.5*stepSize;
        dstrain[3] = -1.0*iPull*stepSize;
        break;
      }
      // run model
      clock_t start = clock();
      (umat)(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, &ddsddt, &drplde,
             &drpldt, strain, dstrain, time, &dTime, &temp, &dTemp, &predef,
             &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
             &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
             &layer, &kspt, &kstep, &kinc);

      clock_t end = clock();
      std::cout << "step = " << iPull << " Time = " << (end - start) << std::endl;
      std::cout << "Strain: " << std::endl;
      printVec(dstrain);
      std::cout << "Stress: " << std::endl;
      printVec(stress);
      myFile << iPull << sep
             << dstrain[0] << sep << dstrain[1] << sep << dstrain[3] << sep
             << stress[0] << sep << stress[1] << sep << stress[3] << sep 
             << statev[13] << sep << statev[14] << sep << statev[15] << sep << statev[16] << sep
             << std::endl;
    }
    myFile.close();
  }
  return 0;
 }
/********************************************************************/
/********************************************************************/

double taylorCheck(const int modelVal) {
  // elastic material definitions
  double E11, E22, E33, nu12, nu23, nu13, G12, G23, G13;
  // damage model parameters
  double hs1_11, hs1_12, hs1_13,
    hs2_22, hs2_12, hs2_23,
    hs4_11, hs4_12, hs4_13, hs4_16,
    hs5_11, hs5_12, hs5_13, hs5_16,
    b1, b2, b6,
    y01, y02, y04, y05,
    yc1, yc2, yc4, yc5,
    pe1, pe2, pe4, pe5,
    dc1, dc2, dc4, dc5;
  const bool pythonTest = true;
  if ( !pythonTest ) {
    // set elasticity
    E11 = 12.416e6; E22 = E11; E33 = 7.513e6;
    nu12 = 0.1106; nu23 = 0.1707; nu13 = nu23;
    G12 = 3.332e6; G23 = 0.670e6; G13 = G23;
    // damage model parameters
    hs1_11 = 2.5; hs1_12 = hs1_11; hs1_13 = hs1_11;
    hs2_22 = 2.5; hs2_12 = hs2_22; hs2_23 = hs2_22;
    hs4_11 = 2.5; hs4_12 = hs4_11; hs4_13 = hs4_11; hs4_16 = hs4_11;
    hs5_11 = 2.5; hs5_12 = hs5_11; hs5_13 = hs5_11; hs5_16 = hs4_11;
    b1 = 0.2; b2 = b1; b6 = b1;
    y01 = 0.005; y02 = y01; y04 = y01; y05 = y01;
    yc1 = 500; yc2 = yc1; yc4 = yc1; yc5 = yc1;
    pe1 = 1.5; pe2 = pe1; pe4 = pe1; pe5 = pe1;
    dc1 = 2.0; dc2 = dc1; dc4 = dc1; dc5 = dc1;
  } else {
    // Elastic Properties
    E11 = 335000;       //Young Modulus [MPa]
    E22 = 180000;       //Young Modulus [MPa]
    E33 = 335000;       //Young Modulus [MPa]
    nu12 = 0.17;         //Poisson ratio
    nu23 = 0.11;         //Poisson ratio
    nu13 = 0.11;         //Poisson ratio
    G12 = 95360;        //Shear modulus [MPa]
    G23 = 95360;        //Shear modulus [MPa]
    G13 = 95360;        //Shear modulus [MPa]
    
    // Material State Variables
    hs1_11 = 1.8;        //Damage Effect
    hs1_12 = 1.0;        //Damage Effect
    hs1_13 = 1.0;        //Damage Effect
    hs2_22 = 20.0;        //Damage Effect
    hs2_12 = 1.0;        //Damage Effect
    hs2_23 = 1.0;        //Damage Effect
    hs4_11 = 2.0;        //Damage Effect
    hs4_12 = 0.5;        //Damage Effect
    hs4_13 = 0.5;        //Damage Effect
    hs4_16 = 1.0;        //Damage Effect
    hs5_11 = 2.0;        //Damage Effect
    hs5_12 = 0.5;        //Damage Effect
    hs5_13 = 0.5;        //Damage Effect
    hs5_16 = 1.0;        //Damage Effect
    b1 = 0.3;           //Traction / Shear Coupling 1
    b2 = 0.3;           //Traction / Shear Coupling 2
    b6 = 0.3;           //Traction / Shear Coupling 2
    y01 = 0.3;         //Damage Thresholds 1
    y02 = 0.02;         //Damage Thresholds 2
    y04 = 0.04;         //Damage Thresholds 1
    y05 = 0.04;         //Damage Thresholds 2
    yc1 = 10;           //Damage Evolution Celerity 1
    yc2 = 1;         //Damage Evolution Celerity 2
    yc4 = 10;           //Damage Evolution Celerity 1
    yc5 = 10;          //Damage Evolution Celerity 2
    pe1 = 1.2;         //Damage Evolution Exponents 1
    pe2 = 3.2;          //Damage Evolution Exponents 2
    pe4 = 1.1;        //Damage Evolution Exponents 1
    pe5 = 1.1;          //Damage Evolution Exponents 2
    dc1 = 3;          //Damage Saturations 1
    dc2 = 8;         //Damage Saturations 2
    dc4 = 8;          //Damage Saturations 1
    dc5 = 8;         //Damage Saturations 2
  }

  // need to set a bunch of stuff
  char cmname[80] = "OpenDM";
  int nprops, nstatv;
  // Set size to be largest possible
  double props[42];
  // elastic props
  props[0] = E11; props[1] = E22; props[2] = E33;
  props[3] = nu12; props[4] = nu23; props[5] = nu13;
  props[6] = G12; props[7] = G23; props[8] = G13;
  if (modelVal == 2) {
    cmname[6] = '2';
    nprops = 21; // Number of properties
    nstatv = 13; // number of statVars
    double props2Param[] = {hs1_11, hs2_22, b1,  b2, y01, y02,
                            yc1, yc2, pe1, pe2, dc1, dc2};
    for (int iProp = 0; iProp < 12; iProp++) {
      props[9+iProp] = props2Param[iProp];
    }
  } else if (modelVal == 4) {
    cmname[6] = '4';
    nprops = 42; // Number of properties
    nstatv = 17; // Number stateVars
    double props4Param[] = {hs1_11, hs2_12, hs1_13,
                            hs2_22, hs2_12, hs2_23,
                            hs4_11, hs4_12, hs4_13, hs4_16,
                            hs5_11, hs5_12, hs5_13, hs5_16,
                            b1, b2, b6,
                            y01, y02, y04, y05,
                            yc1, yc2, yc4, yc5,
                            pe1, pe2, pe4, pe5,
                            dc1, dc2, dc4, dc5};
    for (int iProp = 0; iProp < 33; iProp++) {
      props[9+iProp] = props4Param[iProp];
    }
  }
  double sse, spd, scd, rpl;
  // Current deltaT for which solution is seeking
  double dTime = 0.001;
  double time[] = {0.0, 0.001};

  int ndi = 3;            // 3 direct stress/strain components
  int nshr = 3;           // 3 shear stress/strain components
  int ntens = ndi + nshr; // 6 total componenets. ndi+nshr
  int noel = 1;
  int npt = 0; // Integration point number
  double ddsddt = 0.;
  double drplde = 0.;
  double drpldt = 0.;
  double temp = 0.;
  double dTemp = 0.;
  double predef = 0.;
  double dpred = 0.;
  double coords = 0.;
  double drot = 0.;
  double pnewdt = 0.;
  double dfgrd0 = 0.;
  double dfgrd1 = 0.;
  double celent;
  int layer, kspt, kstep, kinc;
  layer = kspt = kstep = kinc = 0;

  // setup a loop to call umat each step
  // initial stuff
  double strain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrain0[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrain1[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrain2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // stepSize
  double alpha1 = 0.1, alpha2 = 0.2;
  // seed random number generator
  srand(clock());
  double r1, r2;
  int min = 1, max = 10;
  const double dStrainFact = 1.0e-4;
  for (int iDelEps = 0; iDelEps < 1; ++iDelEps) {
    r1 = min + (rand() % (max - min + 1));
    strain[iDelEps] = r1*0.0005;
  }
  r2 = min + (rand() % (max - min + 1));
  dstrain0[0] = r2*dStrainFact;
  dstrain1[0] = alpha1*dstrain0[0];
  dstrain2[0] = alpha2*dstrain0[0];
  double stress0[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsdde0[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double statev0[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double stress1[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsdde1[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double statev1[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double stress2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsdde2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double statev2[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // get sig(eps) and dSigdEps(eps)
  clock_t start = clock();
  (umat)(stress0, statev0, ddsdde0, &sse, &spd, &scd, &rpl, &ddsddt, &drplde,
           &drpldt, strain, dstrain0, time, &dTime, &temp, &dTemp, &predef,
           &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
           &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
           &layer, &kspt, &kstep, &kinc);
  clock_t end = clock();
  std::cout << "epsilonRun" << std::endl;
  std::cout << "Strain0: "  << std::endl;
  printVec(strain);
  std::cout << "Stress: "  << std::endl;
  printVec(stress0);
  std::cout << "Damage: " << statev0[13] << " " << statev0[14] << " "
            << statev0[15] << " " << statev0[16]  << std::endl;
  std::cout << "Elapsed Time: " << (end - start) << std::endl;

  // run model with alpha1 = 1.0
  start = clock();
  (umat)(stress1, statev1, ddsdde1, &sse, &spd, &scd, &rpl, &ddsddt, &drplde,
         &drpldt, strain, dstrain1, time, &dTime, &temp, &dTemp, &predef,
         &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
         &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
         &layer, &kspt, &kstep, &kinc);
  end = clock();
  std::cout << "alpha1Run" << std::endl;
  std::cout << "dStrain: "  << std::endl;
  printVec(dstrain1);
  std::cout << "Stress: "  << std::endl;
  printVec(stress1);
  std::cout << "Damage: " << statev0[13] << " " << statev0[14] << " "
            << statev0[15] << " " << statev0[16]  << std::endl;
  std::cout << "Elapsed Time: " << (end - start) << std::endl;

  // run model with alpha2 = 2.0
  start = clock();
  (umat)(stress2, statev2, ddsdde2, &sse, &spd, &scd, &rpl, &ddsddt, &drplde,
         &drpldt, strain, dstrain2, time, &dTime, &temp, &dTemp, &predef,
         &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
         &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
         &layer, &kspt, &kstep, &kinc);
  
  end = clock();
  std::cout << "alpha2Run" << std::endl;
  std::cout << "dStrain: "  << std::endl;
  printVec(dstrain2);
  std::cout << "Stress: "  << std::endl;
  printVec(stress2);
  std::cout << "Damage: " << statev0[13] << " " << statev0[14] << " "
            << statev0[15] << " " << statev0[16]  << std::endl;
  std::cout << "Elapsed Time: " << (end - start) << std::endl;

  double jDotEps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  matVecProd(ddsdde0, dstrain0, (alpha1 - alpha2), jDotEps);
  std::cout << "jDotEps = " << std::endl;
  printVec(jDotEps);

  double lhsVec[6];
  for (int iComp = 0; iComp < 6; ++iComp) {
    lhsVec[iComp] = stress2[iComp] - stress1[iComp] + jDotEps[iComp];
  }
  const double norm = l2Norm(lhsVec);
  std::cout << norm << std::endl;
  const double convExp = (log(norm))/
                          (log(alpha2) - log(alpha1));
  std::cout << "Convergence Exp = " << convExp << std::endl;
  return 0.0;
}
/********************************************************************/
/********************************************************************/

double numTangCalc(const int modelVal) {
  // elastic material definitions
  double E11, E22, E33, nu12, nu23, nu13, G12, G23, G13;
  // damage model parameters
  double hs1_11, hs1_12, hs1_13,
    hs2_22, hs2_12, hs2_23,
    hs4_11, hs4_12, hs4_13, hs4_16,
    hs5_11, hs5_12, hs5_13, hs5_16,
    b1, b2, b6,
    y01, y02, y04, y05,
    yc1, yc2, yc4, yc5,
    pe1, pe2, pe4, pe5,
    dc1, dc2, dc4, dc5;
  const bool pythonTest = true;
  if ( !pythonTest ) {
    // set elasticity
    E11 = 12.416e6; E22 = E11; E33 = 7.513e6;
    nu12 = 0.1106; nu23 = 0.1707; nu13 = nu23;
    G12 = 3.332e6; G23 = 0.670e6; G13 = G23;
    // damage model parameters
    hs1_11 = 2.5; hs1_12 = hs1_11; hs1_13 = hs1_11;
    hs2_22 = 2.5; hs2_12 = hs2_22; hs2_23 = hs2_22;
    hs4_11 = 2.5; hs4_12 = hs4_11; hs4_13 = hs4_11; hs4_16 = hs4_11;
    hs5_11 = 2.5; hs5_12 = hs5_11; hs5_13 = hs5_11; hs5_16 = hs4_11;
    b1 = 0.2; b2 = b1; b6 = b1;
    y01 = 0.005; y02 = y01; y04 = y01; y05 = y01;
    yc1 = 500; yc2 = yc1; yc4 = yc1; yc5 = yc1;
    pe1 = 1.5; pe2 = pe1; pe4 = pe1; pe5 = pe1;
    dc1 = 2.0; dc2 = dc1; dc4 = dc1; dc5 = dc1;
  } else {
    // Elastic Properties
    E11 = 335000;       //Young Modulus [MPa]
    E22 = 180000;       //Young Modulus [MPa]
    E33 = 335000;       //Young Modulus [MPa]
    nu12 = 0.17;         //Poisson ratio
    nu23 = 0.11;         //Poisson ratio
    nu13 = 0.11;         //Poisson ratio
    G12 = 95360;        //Shear modulus [MPa]
    G23 = 95360;        //Shear modulus [MPa]
    G13 = 95360;        //Shear modulus [MPa]
    
    // Material State Variables
    hs1_11 = 1.8;        //Damage Effect
    hs1_12 = 1.0;        //Damage Effect
    hs1_13 = 1.0;        //Damage Effect
    hs2_22 = 20.0;        //Damage Effect
    hs2_12 = 1.0;        //Damage Effect
    hs2_23 = 1.0;        //Damage Effect
    hs4_11 = 2.0;        //Damage Effect
    hs4_12 = 0.5;        //Damage Effect
    hs4_13 = 0.5;        //Damage Effect
    hs4_16 = 1.0;        //Damage Effect
    hs5_11 = 2.0;        //Damage Effect
    hs5_12 = 0.5;        //Damage Effect
    hs5_13 = 0.5;        //Damage Effect
    hs5_16 = 1.0;        //Damage Effect
    b1 = 0.3;           //Traction / Shear Coupling 1
    b2 = 0.3;           //Traction / Shear Coupling 2
    b6 = 0.3;           //Traction / Shear Coupling 2
    y01 = 0.3;         //Damage Thresholds 1
    y02 = 0.02;         //Damage Thresholds 2
    y04 = 0.04;         //Damage Thresholds 1
    y05 = 0.04;         //Damage Thresholds 2
    yc1 = 10;           //Damage Evolution Celerity 1
    yc2 = 1;         //Damage Evolution Celerity 2
    yc4 = 10;           //Damage Evolution Celerity 1
    yc5 = 10;          //Damage Evolution Celerity 2
    pe1 = 1.2;         //Damage Evolution Exponents 1
    pe2 = 3.2;          //Damage Evolution Exponents 2
    pe4 = 1.1;        //Damage Evolution Exponents 1
    pe5 = 1.1;          //Damage Evolution Exponents 2
    dc1 = 3;          //Damage Saturations 1
    dc2 = 8;         //Damage Saturations 2
    dc4 = 8;          //Damage Saturations 1
    dc5 = 8;         //Damage Saturations 2
  }

  // need to set a bunch of stuff
  char cmname[80] = "OpenDM";
  int nprops, nstatv;
  // Set size to be largest possible
  double props[42];
  // elastic props
  props[0] = E11; props[1] = E22; props[2] = E33;
  props[3] = nu12; props[4] = nu23; props[5] = nu13;
  props[6] = G12; props[7] = G23; props[8] = G13;
  if (modelVal == 2) {
    cmname[6] = '2';
    nprops = 21; // Number of properties
    nstatv = 13; // number of statVars
    double props2Param[] = {hs1_11, hs2_22, b1,  b2, y01, y02,
                            yc1, yc2, pe1, pe2, dc1, dc2};
    for (int iProp = 0; iProp < 12; iProp++) {
      props[9+iProp] = props2Param[iProp];
    }
  } else if (modelVal == 4) {
    cmname[6] = '4';
    nprops = 42; // Number of properties
    nstatv = 17; // Number stateVars
    double props4Param[] = {hs1_11, hs2_12, hs1_13,
                            hs2_22, hs2_12, hs2_23,
                            hs4_11, hs4_12, hs4_13, hs4_16,
                            hs5_11, hs5_12, hs5_13, hs5_16,
                            b1, b2, b6,
                            y01, y02, y04, y05,
                            yc1, yc2, yc4, yc5,
                            pe1, pe2, pe4, pe5,
                            dc1, dc2, dc4, dc5};
    for (int iProp = 0; iProp < 33; iProp++) {
      props[9+iProp] = props4Param[iProp];
    }
  }
  double sse, spd, scd, rpl;
  // Current deltaT for which solution is seeking
  double dTime = 0.001;
  double time[] = {0.0, 0.001};

  int ndi = 3;            // 3 direct stress/strain components
  int nshr = 3;           // 3 shear stress/strain components
  int ntens = ndi + nshr; // 6 total componenets. ndi+nshr
  int noel = 1;
  int npt = 0; // Integration point number
  double ddsddt = 0.;
  double drplde = 0.;
  double drpldt = 0.;
  double temp = 0.;
  double dTemp = 0.;
  double predef = 0.;
  double dpred = 0.;
  double coords = 0.;
  double drot = 0.;
  double pnewdt = 0.;
  double dfgrd0 = 0.;
  double dfgrd1 = 0.;
  double celent;
  int layer, kspt, kstep, kinc;
  layer = kspt = kstep = kinc = 0;

  // setup a loop to call umat each step
  // initial stuff
  double strain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dstrainPert[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // seed random number generator
  srand(clock());
  double r1, r2;
  int min = 1, max = 10;
  const double dStrainFact = 1.0e-8;
  r2 = min + (rand() % (max - min + 1));
  for (int iDelEps = 0; iDelEps < 6; ++iDelEps) {
    r1 = min + (rand() % (max - min + 1));
    strain[iDelEps] = r1*0.0001;
    dstrain[iDelEps] = 0.0;
  }
  
  double stress[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsdde[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double statev[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double ddsddePert[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double statevPert[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // get sig(eps) and dSigdEps(eps)
  clock_t start = clock();
  (umat)(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, &ddsddt, &drplde,
         &drpldt, strain, dstrain, time, &dTime, &temp, &dTemp, &predef,
         &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
         &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
         &layer, &kspt, &kstep, &kinc);
  clock_t end = clock();
  std::cout << "epsilonRun" << std::endl;
  std::cout << "Strain: "  << std::endl;
  printVec(strain);
  std::cout << "Stress: "  << std::endl;
  printVec(stress);
  std::cout << "Damage: " << statev[13] << " " << statev[14] << " "
            << statev[15] << " " << statev[16]  << std::endl;
  std::cout << "Elapsed Time: " << (end - start) << std::endl;

  
  double stressPert[6];
  double allStressPert[36];
  // run perturbations
  for (int iPert = 0; iPert < 6; iPert++) {
    // Preload stresses and strains
    for (int iStress = 0; iStress < 6; ++iStress) {
      stressPert[iStress] = stress[iStress];
      dstrainPert[iStress] = 0.0;
    }
    dstrainPert[iPert] = r2*dStrainFact;
    // setup stateVars
    for (int iSV = 0; iSV < 17; ++iSV) {
      statevPert[iSV] = statev[iSV];
    }
    start = clock();
    (umat)(stressPert, statevPert, ddsddePert, &sse, &spd,
           &scd, &rpl, &ddsddt, &drplde,
           &drpldt, strain, dstrainPert, time, &dTime, &temp, &dTemp, &predef,
           &dpred, cmname, &ndi, &nshr, &ntens, &nstatv, props, &nprops,
           &coords, &drot, &pnewdt, &celent, &dfgrd0, &dfgrd1, &noel, &npt,
           &layer, &kspt, &kstep, &kinc);
    end = clock();
    std::cout << "alpha1Run" << std::endl;
    std::cout << "dStrain: "  << std::endl;
    printVec(dstrainPert);
    std::cout << "Stress: "  << std::endl;
    printVec(stressPert);
    std::cout << "Damage: " << statevPert[13] << " " << statevPert[14] << " "
              << statevPert[15] << " " << statevPert[16]  << std::endl;
    std::cout << "Elapsed Time: " << (end - start) << std::endl;
    // unload stresses
    for (int iStress = 0; iStress < 6; ++iStress) {
      allStressPert[6*iPert+iStress] = stressPert[iStress] - stress[iStress];
    }
    std::cout << "StressPert: "  << std::endl;
    printVec(&allStressPert[6*iPert]);
  }
  // estimate tangent
  double numTang[36];
  numericalTangent(stress, allStressPert, r2*dStrainFact, numTang);
  std::cout << "Numerical Tangent" << std::endl;
  printMat(numTang);
  std::cout << "Math Tangent" << std::endl;
  printMat(ddsdde);
  
  return 0.0;
}
/********************************************************************/
/********************************************************************/

double l2Norm(const double* v1) {
  // Calc ||v||_2
  double norm = 0.0;
  for (int iComp = 0; iComp < 6; ++iComp) {
    norm += v1[iComp]*v1[iComp];
  }
  return std::sqrt(norm);
}
// /********************************************************************/
// /********************************************************************/

void matVecProd(const double* A, const double* b,
                const double alpha, double* c) {
  // c_{i} = alpha*A_{ij} b_{j}
  for (int iRow = 0; iRow < 6; ++iRow) {
    for (int jCol = 0; jCol < 6; ++jCol) {
      c[iRow] += alpha*A[6*iRow+jCol]*b[jCol];
    }
  }
}
/********************************************************************/
/********************************************************************/

void printVec(const double* v) {
  for (int iRow = 0; iRow < 6; ++iRow) {
    std::cout << v[iRow] << " ";
  }
  std::cout << std::endl;
}
/********************************************************************/
/********************************************************************/

void printMat(const double* A) {
  for (int iRow = 0; iRow < 6; ++iRow) {
    for (int jCol = 0; jCol < 6; ++jCol) {
      std::cout << A[6*iRow+jCol] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
// /********************************************************************/
// /********************************************************************/

void numericalTangent(const double* stress0, const double* stressPert,
                      const double strainInc, double* numTang) {
  for (int iStress = 0; iStress < 6; ++iStress) {
    for (int jPert = 0; jPert < 6; ++jPert) {
      numTang[6*iStress+jPert] = stressPert[iStress+6*jPert]/strainInc;
    }
  }
}
