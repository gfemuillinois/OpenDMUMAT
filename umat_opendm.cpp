// Implementation of OpenDM as a UMAT
// by: Bryce Mazurowski (brycepm2@gmail.com)
//

#include <iostream>
#include <string>
#include <ctime>
#include <fstream>

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
    loadApplied = 10000;
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
    loadApplied = 3350;
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
  double statev[17] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
    for (int iProp = 0; iProp < 31; iProp++) {
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
  double stepSize = loadApplied / (nInc);
  double S11 = 1.0 / E11, S12 = -nu12 / E11, S22 = 1.0 / E22, S44 = 1.0 / G12;
  std::ofstream myFile;
  std::string modelNum = std::to_string(modelVal);
  std::string fileName = "umatTestOut_"+modelNum+".csv";
  myFile.open(fileName);
  myFile << "step, strain, stress, d1" << std::endl;
  myFile << 0 << ", " << 0.0 << ", " << 0.0 << std::endl;
  const bool tensRun = true;
  for (int iPull = 1; iPull <= nInc; iPull++) {
    // Increment dstrain
    if (tensRun) {
      dstrain[0] = iPull * stepSize * (S11);
      // dstrain[1] = iPull * stepSize * (S12);
    } else {
      dstrain[0] = iPull * stepSize * (S11 + S12);
      dstrain[1] = iPull * stepSize * (S12 + S22);
      dstrain[3] = iPull * stepSize * (S44);
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
    std::cout << "Strain: " << (*dstrain) << std::endl;
    std::cout << "Stress: " << (*stress) << std::endl;
    myFile << iPull << ", " << dstrain[0] << ", " << stress[0] << ", "
           << statev[13] << std::endl;
   }
   myFile.close();
 }
 /********************************************************************/
 /********************************************************************/
