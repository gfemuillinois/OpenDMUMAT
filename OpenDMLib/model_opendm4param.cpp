// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// Implementation file for OpenDM Model with 4 damage parameters
// This is derived from OpenDMModel class
// It captures damage in the 0 & 90 directions in plane of composite
// as well as shear damage
//

#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "model_opendm4param.hpp"

using std::sqrt;
using std::cout;
using std::endl;

/********************************************************************/
/********************************************************************/

OpenDMModel4Param::OpenDMModel4Param(double* props, int* nprops,
                     double* statev, int* nstatv)
    : OpenDMModel(props, nprops, statev, nstatv, 4 /*nDamageVals*/)
{
  // Set OpenDM model params
  unpackParams(props);

  // Make transformation matrices for shear params
  createTEpsMats();
  
  // Make H matrices without stress activation
  createHMats();

  // initialize spectral form work areas 
  // set in posPartStrainD*, reused in dEpsD*PlusDEps
  eValsD1 = Matrix3d::Zero(); eValsD2 = Matrix3d::Zero();
  eVectsD1 = Matrix3d::Zero();
  eVectsD2 = Matrix3d::Zero();
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::unpackParams(double* props) {
  // props[9:41] =
  // hs1_11, hs1_12, hs1_13,
  // hs2_22, hs2_12, hs2_23,
  // hs4_11, hs4_12, hs4_13, hs4_16,
  // hs5_11, hs5_12, hs5_13, hs5_16
  // b_1, b_2, b_6,
  // y01, y02, y04, y05, yc1, yc2, yc4, yc5,
  // pe1, pe2, pe4, pe5, dc1, dc2, dc4, dc5
  // unpack props
  hs1 = Eigen::Map<Eigen::VectorXd>(props+9,3);
  hs2 = Eigen::Map<Eigen::VectorXd>(props+12,3);
  hs4 = Eigen::Map<Eigen::VectorXd>(props+15,4);
  hs5 = Eigen::Map<Eigen::VectorXd>(props+19,4);
  // std::cout << "hs1 = " << hs1 << "\n"
            // << "hs2 = " << hs2 << "\n"
            // << "hs4 = " << hs4 << "\n"
            // << "hs5 = " << hs5 << std::endl;
  b = Eigen::Map<Eigen::VectorXd>(props+23,3);
  y0 = Eigen::Map<Eigen::VectorXd>(props+26,4);
  yc = Eigen::Map<Eigen::VectorXd>(props+30,4);
  // std::cout << "yc = " << yc << std::endl;
  pe = Eigen::Map<Eigen::VectorXd>(props+34,4);
  // std::cout << "pe = " << pe << std::endl;
  dc = Eigen::Map<Eigen::VectorXd>(props+38,4);
  // std::cout << "dc = " << dc << std::endl;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::createTEpsMats() {
  // Transformation matrix
  // +45 transform
  double invSqrt2 = 1.0/std::sqrt(2.0);
  Teps_p45(0,0) = 0.5; Teps_p45(0,1) = 0.5; Teps_p45(0,3) = 0.5; 
  Teps_p45(1,0) = 0.5; Teps_p45(1,1) = 0.5; Teps_p45(1,3) = -0.5;
  Teps_p45(2,2) = 1.0;
  Teps_p45(3,0) = -1.0; Teps_p45(3,1) = 1.0;
  Teps_p45(4,4) = invSqrt2; Teps_p45(4,5) = invSqrt2;
  Teps_p45(5,4) = -invSqrt2; Teps_p45(5,5) = invSqrt2;

  // -45 transform
  Teps_n45(0,0) = 0.5; Teps_n45(0,1) = 0.5; Teps_n45(0,3) = -0.5;
  Teps_n45(1,0) = 0.5; Teps_n45(1,1) = 0.5; Teps_n45(1,3) = 0.5;
  Teps_n45(2,2) = 1.0;
  Teps_n45(3,0) = 1.0; Teps_n45(3,1) = -1.0;
  Teps_n45(4,4) = invSqrt2; Teps_n45(4,5) = -invSqrt2;
  Teps_n45(5,4) = invSqrt2; Teps_n45(5,5) = invSqrt2;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::createHMats() {
  // Make H's
  // H1
  // Mode I
  H1 = Matrix6d::Zero();
  H1(0,0) = hs1(0)*S0(0,0);
  // Mode III
  H1(3,3) = hs1(1)*S0(3,3);
  // Mode II
  H1(4,4) = hs1(2)*S0(4,4);

  // H2
  // Mode I
  H2 = Matrix6d::Zero();
  H2(1,1) = hs2(0)*S0(1,1);
  // Mode II
  H2(3,3) = hs2(1)*S0(3,3);
  // Mode III
  H2(5,5) = hs2(2)*S0(5,5);

  // transform S0 for H4 & H5
  // NOTE: It so happens that:
  Matrix6d S_p45, S_n45;
  S_p45 = Teps_p45*S0*Teps_p45.transpose();
  S_n45 = Teps_n45*S0*Teps_n45.transpose();

  // H4
  // Mode I
  H4 = Matrix6d::Zero();
  H4(0,0) = hs4(0)*S_p45(0,0);
  // Mode II
  H4(3,3) = hs4(1)*S_p45(3,3);
  // Mode III
  H4(4,4) = hs4(2)*S_p45(4,4);
  // Transformation induced pieces
  H4(0,3) = hs4(3)*S_p45(0,3);
  H4(0,4) = hs4(3)*S_p45(0,4);
  H4(3,0) = hs4(3)*S_p45(3,0);
  H4(4,0) = hs4(3)*S_p45(4,0);
  // Bring back to matCoords
  H4 = Teps_n45*H4*Teps_n45.transpose();

  // H5
  // Mode I
  H5 = Matrix6d::Zero();
  H5(0,0) = hs5(0)*S_n45(0,0);
  // Mode II
  H5(3,3) = hs5(1)*S_n45(3,3);
  // Mode III
  H5(4,4) = hs5(2)*S_n45(4,4);
  // Transformation induced pieces
  H5(0,3) = hs5(3)*S_n45(0,3);
  H5(0,4) = hs5(3)*S_n45(0,4);
  H5(3,0) = hs5(3)*S_n45(3,0);
  H5(4,0) = hs5(3)*S_n45(4,0);
  // Bring back to matCoords
  H5 = Teps_p45*H5*Teps_p45.transpose();
}
/********************************************************************/
/********************************************************************/

VectorXd OpenDMModel4Param::calcDrivingForces(const Vector6d& epsStar) {
  // EQ 46-55 in OpenDM-4 Parameter Damage CLT doc
  Vector6d epsD1 = Vector6d::Zero(), epsD1Plus = Vector6d::Zero(),
    epsD2 = Vector6d::Zero(), epsD2Plus = Vector6d::Zero();
  epsD1(0) = epsStar(0); epsD1(3) = epsStar(3); epsD1(4) = epsStar(4);
  epsD2(1) = epsStar(1); epsD2(3) = epsStar(3); epsD2(5) = epsStar(5);
  posPartStrainD1(epsD1, epsD1Plus);
  posPartStrainD2(epsD2, epsD2Plus);
  
  // Base driving force
  double z1 = 0.5*(epsD1Plus(0)*C0(0,0)*epsD1Plus(0)
           + b(0)*epsD1Plus(3)*C0(3,3)*epsD1Plus(3)
           + b(1)*epsD1Plus(4)*C0(4,4)*epsD1Plus(4));
  double z2 = 0.5*(epsD2Plus(1)*C0(1,1)*epsD2Plus(1)
           + b(0)*epsD2Plus(3)*C0(3,3)*epsD2Plus(3)
           + b(1)*epsD2Plus(5)*C0(5,5)*epsD2Plus(5));
  double z6 = 0.25*(epsD1Plus(0)*C0(0,0)*epsD1Plus(3)
            + epsD2Plus(1)*C0(1,1)*epsD2Plus(3)
            + b(2)*C0(3,3)*(epsD1Plus(3)*epsD1Plus(0) + 
                            epsD2Plus(3)*epsD2Plus(1)));

  // Driving forces for each damage variable
  double y1 = z1 - std::abs(z6);
  double y2 = z2 - std::abs(z6);
  double y4 = macaulayBracketPlus(z6);
  double y5 = -1.0*macaulayBracketMinus(z6); // must be positive quantity

  VectorXd yMax(4);
  yMax(0) = y1 > yMaxSave(0) ? y1 : yMaxSave(0);
  yMax(1) = y2 > yMaxSave(1) ? y2 : yMaxSave(1);
  yMax(2) = y4 > yMaxSave(2) ? y4 : yMaxSave(2);
  yMax(3) = y5 > yMaxSave(3) ? y5 : yMaxSave(3);
  // std::cout << "z = " << z1 << " " << z2 << " "
            // << z6 << std::endl;
  // std::cout << "yMax = " << yMax << std::endl;
  return yMax;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::posPartStrainD1(const Vector6d& epsD1,
                                        Vector6d& epsD1Plus) {

  // get strains
  const double e11 = epsD1(0), gam12 = epsD1(3), gam13 = epsD1(4);

  // pick given strain state
  const bool hasE11 = (e11 != 0.0), hasE12 = (gam12 != 0.0),
    hasE13 = (gam13 != 0.0);
  // useful terms
  const double rootE11Gam12Gam13 =
    std::sqrt(e11*e11 + gam12*gam12 + gam13*gam13);
  if (hasE11 && hasE12 && hasE13) {
    // case 1 - all D1 strains present
    eValsD1(1,1) = 0.5*(e11 - rootE11Gam12Gam13);
    eValsD1(2,2) = 0.5*(e11 + rootE11Gam12Gam13);

    // v1
    const double v1Denom = std::sqrt(gam12*gam12 + gam13*gam13);
    eVectsD1(1,0) = -gam13*std::abs(gam12)/(gam12*v1Denom);
    eVectsD1(2,0) = std::abs(gam12)/v1Denom;

    // v2
    const double v2Denom = std::sqrt(gam12*gam12 + gam13*gam13 +
                                     (e11 - rootE11Gam12Gam13)*(e11 - rootE11Gam12Gam13));
    eVectsD1(0,1) = (e11 - rootE11Gam12Gam13)*std::abs(gam13)/(gam13*v2Denom);
    eVectsD1(1,1) = gam12*std::abs(gam13)/(gam13*v2Denom);
    eVectsD1(2,1) = std::abs(gam13)/v2Denom;

    // v3
    const double v3Denom = std::sqrt(gam12*gam12 + gam13*gam13 +
                                     (e11 + rootE11Gam12Gam13)*(e11 + rootE11Gam12Gam13));
    eVectsD1(0,2) = (e11 + rootE11Gam12Gam13)*std::abs(gam13)/(gam13*v3Denom);
    eVectsD1(1,2) = gam12*std::abs(gam13)/(gam13*v3Denom);
    eVectsD1(2,2) = std::abs(gam13)/v3Denom;
  } else if (hasE11 && hasE12 && !hasE13) {
    // case 2 - E13 is zero
    eValsD1(1,1) = 0.5*(e11 - rootE11Gam12Gam13);
    eValsD1(2,2) = 0.5*(e11 + rootE11Gam12Gam13);

    // v1
    eVectsD1(2,0) = 1.0;

    // v2
    const double v2Denom = std::sqrt(gam12*gam12 +
                                     (e11 - rootE11Gam12Gam13)*(e11 - rootE11Gam12Gam13));
    eVectsD1(0,1) = (e11 - rootE11Gam12Gam13)*std::abs(gam12)/(gam12*v2Denom);
    eVectsD1(1,1) = std::abs(gam12)/(v2Denom);

    // v3
    const double v3Denom = std::sqrt(gam12*gam12 +
                                     (e11 + rootE11Gam12Gam13)*(e11 + rootE11Gam12Gam13));
    eVectsD1(0,2) = (e11 + rootE11Gam12Gam13)*std::abs(gam12)/(gam12*v3Denom);
    eVectsD1(1,2) = std::abs(gam12)/(v3Denom);
  } else if (hasE11 && !hasE12 && hasE13) {
    // case 3 - E12 is zero
    eValsD1(1,1) = 0.5*(e11 - rootE11Gam12Gam13);
    eValsD1(2,2) = 0.5*(e11 + rootE11Gam12Gam13);

    // v1
    eVectsD1(1,0) = 1.0;

    // v2
    const double v2Denom = std::sqrt(gam13*gam13 +
                                     (e11 - rootE11Gam12Gam13)*(e11 - rootE11Gam12Gam13));
    eVectsD1(0,1) = (e11 - rootE11Gam12Gam13)*std::abs(gam13)/(gam13*v2Denom);
    eVectsD1(2,1) = std::abs(gam13)/(v2Denom);

    // v3
    const double v3Denom = std::sqrt(gam13*gam13 +
                                     (e11 + rootE11Gam12Gam13)*(e11 + rootE11Gam12Gam13));
    eVectsD1(0,2) = (e11 + rootE11Gam12Gam13)*std::abs(gam13)/(gam13*v3Denom);
    eVectsD1(2,2) = std::abs(gam13)/(v3Denom);
  } else if (!hasE11 && hasE12 && hasE13) {
    // case 4 - E11 is zero
    eValsD1(1,1) = -0.5*rootE11Gam12Gam13;
    eValsD1(2,2) = 0.5*rootE11Gam12Gam13;

    // v1
    eVectsD1(1,0) = -1.0*gam13*std::abs(gam12)/(gam12*rootE11Gam12Gam13);
    eVectsD1(2,0) = std::abs(gam12)/(rootE11Gam12Gam13);

    // v2
    const double denom = std::sqrt(2.0*(gam12*gam12 + gam13*gam13));
    eVectsD1(0,1) = -1.0*rootE11Gam12Gam13*std::abs(gam13)/(gam13*denom);
    eVectsD1(1,1) = gam12*std::abs(gam13)/(gam13*denom);
    eVectsD1(2,1) = std::abs(gam13)/denom;

    //v3
    eVectsD1(0,2) = rootE11Gam12Gam13*std::abs(gam13)/(gam13*denom);
    eVectsD1(1,2) = gam12*std::abs(gam13)/(gam13*denom);
    eVectsD1(2,2) = std::abs(gam13)/denom;
    
  } else if (hasE11 && !hasE12 && !hasE13) {
    // case 5 - E12 and E13 are zero
    eValsD1(1,1) = e11;
    // v1
    eVectsD1(1,0) = 1.0;
    // v2
    eVectsD1(0,1) = 1.0;
    // v3
    eVectsD1(2,2) = 1.0;
  } else if (!hasE11 && hasE12 && !hasE13) {
    // case 6 - E11 and E13 are zero
    eValsD1(1,1) = -0.5*gam12;
    eValsD1(2,2) = 0.5*gam12;
    // v1
    eVectsD1(2,0) = 1.0;
    // v2
    eVectsD1(0,1) = -0.707106781186547;
    eVectsD1(1,1) = 0.707106781186547;
    // v3
    eVectsD1(0,2) = 0.707106781186547;
    eVectsD1(1,2) = 0.707106781186547;
  } else if (!hasE11 && !hasE12 && hasE13) {
    // case 7 - E11 and E12 are zero
    eValsD1(1,1) = -0.5*gam13;
    eValsD1(2,2) = 0.5*gam13;
    // v1
    eVectsD1(1,0) = 1.0;
    // v2
    eVectsD1(0,1) = -0.707106781186547;
    eVectsD1(2,1) = 0.707106781186547;
    // v3
    eVectsD1(0,2) = 0.707106781186547;
    eVectsD1(2,2) = 0.707106781186547;
  } else {
      // no strains, set eVects to global coord
      eVectsD1(0, 0) = 1.0;
      eVectsD1(1, 1) = 1.0;
      eVectsD1(2, 2) = 1.0;
  }

  // positive part of eValsD1
  eValsD1(0,0) = macaulayBracketPlus(eValsD1(0,0));
  eValsD1(1,1) = macaulayBracketPlus(eValsD1(1,1));
  eValsD1(2,2) = macaulayBracketPlus(eValsD1(2,2));

  // transform from spectral -> global
  Matrix3d epsilonPlus = eVectsD1*eValsD1*eVectsD1.transpose();

  // Put back into Abaqus notation
  epsD1Plus(0) = epsilonPlus(0,0);
  epsD1Plus(1) = epsilonPlus(1,1);
  epsD1Plus(2) = epsilonPlus(2,2);
  epsD1Plus(3) = 2.0*epsilonPlus(0,1);
  epsD1Plus(4) = 2.0*epsilonPlus(0,2);
  epsD1Plus(5) = 2.0*epsilonPlus(1,2);
}

/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::posPartStrainD2(const Vector6d& epsD2,
                                        Vector6d& epsD2Plus) {

  // get strains
  const double e22 = epsD2(1), gam12 = epsD2(3), gam23 = epsD2(5);

  // pick given strain state
  const bool hasE22 = (e22 != 0.0), hasE12 = (gam12 != 0.0),
    hasE23 = (gam23 != 0.0);
  // useful terms
  const double rootE22Gam12Gam23 =
    std::sqrt(e22*e22 + gam12*gam12 + gam23*gam23);
  if (hasE22 && hasE12 && hasE23) {
    // case 1 - all D2 strains present
    eValsD2(1,1) = 0.5*(e22 - rootE22Gam12Gam23);
    eValsD2(2,2) = 0.5*(e22 + rootE22Gam12Gam23);

    // v1
    const double v1Denom = std::sqrt(gam12*gam12 + gam23*gam23);
    eVectsD2(1,0) = -gam23*std::abs(gam12)/(gam12*v1Denom);
    eVectsD2(2,0) = std::abs(gam12)/v1Denom;

    // v2
    const double v2Denom = std::sqrt(gam12*gam12 + gam23*gam23 +
                                     (e22 - rootE22Gam12Gam23)*(e22 - rootE22Gam12Gam23));
    eVectsD2(0,1) = gam12*std::abs(gam23)/(gam23*v2Denom);
    eVectsD2(1,1) = (e22 - rootE22Gam12Gam23)*std::abs(gam23)/(gam23*v2Denom);
    eVectsD2(2,1) = std::abs(gam23)/v2Denom;

    // v3
    const double v3Denom = std::sqrt(gam12*gam12 + gam23*gam23 +
                                     (e22 + rootE22Gam12Gam23)*(e22 + rootE22Gam12Gam23));
    eVectsD2(0,2) = gam12*std::abs(gam23)/(gam23*v3Denom);
    eVectsD2(1,2) = (e22 + rootE22Gam12Gam23)*std::abs(gam23)/(gam23*v3Denom);
    eVectsD2(2,2) = std::abs(gam23)/v3Denom;
  } else if (hasE22 && hasE12 && !hasE23) {
    // case 2 - E23 is zero
    eValsD2(1,1) = 0.5*(e22 - rootE22Gam12Gam23);
    eValsD2(2,2) = 0.5*(e22 + rootE22Gam12Gam23);

    // v1
    eVectsD2(2,0) = 1.0;

    // v2
    const double v2Denom = std::sqrt(gam12*gam12 +
                                     (e22 + rootE22Gam12Gam23)*(e22 + rootE22Gam12Gam23));
    eVectsD2(0,1) = -(e22 + rootE22Gam12Gam23)*std::abs(gam12)/(gam12*v2Denom);
    eVectsD2(1,1) = std::abs(gam12)/(v2Denom);

    // v3
    const double v3Denom = std::sqrt(gam12*gam12 +
                                     (e22 - rootE22Gam12Gam23)*(e22 - rootE22Gam12Gam23));
    eVectsD2(0,2) = -(e22 - rootE22Gam12Gam23)*std::abs(gam12)/(gam12*v3Denom);
    eVectsD2(1,2) = std::abs(gam12)/(v3Denom);
  } else if (hasE22 && !hasE12 && hasE23) {
    // case 3 - E12 is zero
    eValsD2(1,1) = 0.5*(e22 - rootE22Gam12Gam23);
    eValsD2(2,2) = 0.5*(e22 + rootE22Gam12Gam23);

    // v1
    eVectsD2(1,0) = 1.0;

    // v2
    const double v2Denom = std::sqrt(gam23*gam23 +
                                     (e22 - rootE22Gam12Gam23)*(e22 - rootE22Gam12Gam23));
    eVectsD2(1,1) = (e22 - rootE22Gam12Gam23)*std::abs(gam23)/(gam23*v2Denom);
    eVectsD2(2,1) = std::abs(gam23)/(v2Denom);

    // v3
    const double v3Denom = std::sqrt(gam23*gam23 +
                                     (e22 + rootE22Gam12Gam23)*(e22 + rootE22Gam12Gam23));
    eVectsD2(1,2) = (e22 + rootE22Gam12Gam23)*std::abs(gam23)/(gam23*v3Denom);
    eVectsD2(2,2) = std::abs(gam23)/(v3Denom);
  } else if (!hasE22 && hasE12 && hasE23) {
    // case 4 - E22 is zero
    eValsD2(1,1) = -0.5*rootE22Gam12Gam23;
    eValsD2(2,2) = 0.5*rootE22Gam12Gam23;

    // v1
    eVectsD2(1,0) = -1.0*gam23*std::abs(gam12)/(gam12*rootE22Gam12Gam23);
    eVectsD2(2,0) = std::abs(gam12)/(rootE22Gam12Gam23);

    // v2
    const double denom = std::sqrt(2.0*(gam12*gam12 + gam23*gam23));
    eVectsD2(0,1) = gam12*std::abs(gam23)/(gam23*denom);
    eVectsD2(1,1) = -1.0*rootE22Gam12Gam23*std::abs(gam23)/(gam23*denom);
    eVectsD2(2,1) = std::abs(gam23)/denom;

    //v3
    eVectsD2(0,2) = gam12*std::abs(gam23)/(gam23*denom);
    eVectsD2(1,2) = rootE22Gam12Gam23*std::abs(gam23)/(gam23*denom);
    eVectsD2(2,2) = std::abs(gam23)/denom;
    
  } else if (hasE22 && !hasE12 && !hasE23) {
    // case 5 - E12 and E23 are zero
    eValsD2(1,1) = e22;
    // v1
    eVectsD2(0,0) = 1.0;
    // v2
    eVectsD2(1,1) = 1.0;
    // v3
    eVectsD2(2,2) = 1.0;
  } else if (!hasE22 && hasE12 && !hasE23) {
    // case 6 - E22 and E13 are zero
    eValsD2(1,1) = -0.5*gam12;
    eValsD2(2,2) = 0.5*gam12;
    // v1
    eVectsD2(2,0) = 1.0;
    // v2
    eVectsD2(0,1) = -0.707106781186547;
    eVectsD2(1,1) = 0.707106781186547;
    // v3
    eVectsD2(0,2) = 0.707106781186547;
    eVectsD2(1,2) = 0.707106781186547;
  } else if (!hasE22 && !hasE12 && hasE23) {
    // case 7 - E22 and E12 are zero
    eValsD2(1,1) = -0.5*gam23;
    eValsD2(2,2) = 0.5*gam23;
    // v1
    eVectsD2(1,0) = 1.0;
    // v2
    eVectsD2(1,1) = -0.707106781186547;
    eVectsD2(2,1) = 0.707106781186547;
    // v3
    eVectsD2(1,2) = 0.707106781186547;
    eVectsD2(2,2) = 0.707106781186547;
  } else {
    // no strains, set eVects to global coord
    eVectsD1(0, 0) = 1.0;
    eVectsD1(1, 1) = 1.0;
    eVectsD1(2, 2) = 1.0;
  }


  // positive part of eValsD2
  eValsD2(0,0) = macaulayBracketPlus(eValsD2(0,0));
  eValsD2(1,1) = macaulayBracketPlus(eValsD2(1,1));
  eValsD2(2,2) = macaulayBracketPlus(eValsD2(2,2));

  // transform from spectral -> global
  Matrix3d epsilonPlus = eVectsD2*eValsD2*eVectsD2.transpose();

  // Put back into Abaqus notation
  epsD2Plus(0) = epsilonPlus(0,0);
  epsD2Plus(1) = epsilonPlus(1,1);
  epsD2Plus(2) = epsilonPlus(2,2);
  epsD2Plus(3) = 2.0*epsilonPlus(0,1);
  epsD2Plus(4) = 2.0*epsilonPlus(0,2);
  epsD2Plus(5) = 2.0*epsilonPlus(1,2);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::calcDEpsD1PlusDEps(const Vector6d& epsD1,
                                           Matrix6d& dEpsD1PlusDEps) const {
  // This funciton calcs derivative of epsilonD1+ wrt epsilon
  // It uses a number of simplifications, but should be sufficiently capable
  // ASSUMPTIONS:
  // eValsD1(0,0) = 0.0
  // eValsD1(i,j) = 0.0 for all i != j
  //
  // following this basic setup:
  // epsilonD1+ = L_{ij} <\hat{\epsilon}_j^{D1}>_+

  // get strains
  const double epsilon11 = epsD1(0), gamma12 = epsD1(3), gamma13 = epsD1(4);

  // pick given strain state
  const bool hasE11 = (epsilon11 != 0.0), hasE12 = (gamma12 != 0.0),
    hasE13 = (gamma13 != 0.0);

  // Transformation matrix spectral -> Global
  // since I know I have limited \hat{\epsilon_i^{D1+}}
  // I only need part of this
  double L12 = eVectsD1(1,0)*eVectsD1(1,0), L13 = eVectsD1(2,0)*eVectsD1(2,0),
    L22 = eVectsD1(1,1)*eVectsD1(1,1), L23 = eVectsD1(2,1)*eVectsD1(2,1),
    L32 = eVectsD1(1,2)*eVectsD1(1,2), L33 = eVectsD1(2,2)*eVectsD1(2,2),
    L42 = 2.0*eVectsD1(1,0)*eVectsD1(1,1), L43 = 2.0*eVectsD1(2,0)*eVectsD1(2,1),
    L52 = 2.0*eVectsD1(1,0)*eVectsD1(1,2), L53 = 2.0*eVectsD1(2,0)*eVectsD1(2,2),
    L62 = 2.0*eVectsD1(1,1)*eVectsD1(1,2), L63 = 2.0*eVectsD1(2,1)*eVectsD1(2,2);
  // initialize d L_{ij}/d \epsilon_{k}
  Vector6d dL12dEps = Vector6d::Zero(), dL13dEps = Vector6d::Zero(), 
    dL22dEps = Vector6d::Zero(), dL23dEps = Vector6d::Zero(), 
    dL32dEps = Vector6d::Zero(), dL33dEps = Vector6d::Zero(), 
    dL42dEps = Vector6d::Zero(), dL43dEps = Vector6d::Zero(), 
    dL52dEps = Vector6d::Zero(), dL53dEps = Vector6d::Zero(), 
    dL62dEps = Vector6d::Zero(), dL63dEps = Vector6d::Zero();

  // derivative of eVals
  Vector6d dEVal2dEps = Vector6d::Zero(), dEVal3dEps = Vector6d::Zero();

  // useful terms
  const double eps11Sq = epsilon11*epsilon11, gam12Sq = gamma12*gamma12,
    gam13Sq = gamma13*gamma13;
  const double rootE11Gam12Gam13 =
    std::sqrt(eps11Sq + gam12Sq + gam13Sq);
  if (hasE11 && hasE12 && hasE13) {
    // case 1 - all D1 strains present
    // d L_{ij}/d \epsilon_k for this case
    // \epsilon_{11}
    dL22dEps(0) =  0.50000000000000011*gam12Sq*1.0/rootE11Gam12Gam13*(eps11Sq + rootE11Gam12Gam13*(-2.0*epsilon11 + rootE11Gam12Gam13))/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL32dEps(0) =  0.50000000000000011*gam12Sq*1.0/rootE11Gam12Gam13*(-eps11Sq - rootE11Gam12Gam13*(2.0*epsilon11 + rootE11Gam12Gam13))/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL42dEps(0) =  -0.70710678118654757*1.0/rootE11Gam12Gam13*(eps11Sq - rootE11Gam12Gam13*(2.0*epsilon11 - rootE11Gam12Gam13))*fabs(gamma12*gamma13)/(sqrt(1.0*gam12Sq + gam13Sq)*pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    dL52dEps(0) =  0.70710678118654757*1.0/rootE11Gam12Gam13*(eps11Sq + rootE11Gam12Gam13*(2.0*epsilon11 + rootE11Gam12Gam13))*fabs(gamma12*gamma13)/(sqrt(1.0*gam12Sq + gam13Sq)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    dL62dEps(0) =  -1.0000000000000002*epsilon11*pow(gam12Sq, 2)/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)) - 1.0000000000000002*epsilon11*gam12Sq*gam13Sq/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    dL23dEps(0) =  0.50000000000000011*gam13Sq*1.0/rootE11Gam12Gam13*(eps11Sq + rootE11Gam12Gam13*(-2.0*epsilon11 + rootE11Gam12Gam13))/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL33dEps(0) =  0.50000000000000011*gam13Sq*1.0/rootE11Gam12Gam13*(-eps11Sq - rootE11Gam12Gam13*(2.0*epsilon11 + rootE11Gam12Gam13))/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL43dEps(0) =  0.70710678118654757*1.0/rootE11Gam12Gam13*(eps11Sq + rootE11Gam12Gam13*(-2.0*epsilon11 + rootE11Gam12Gam13))*fabs(gamma12*gamma13)/(sqrt(1.0*gam12Sq + gam13Sq)*pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    dL53dEps(0) =  -0.70710678118654757*1.0/rootE11Gam12Gam13*(eps11Sq + rootE11Gam12Gam13*(2.0*epsilon11 + rootE11Gam12Gam13))*fabs(gamma12*gamma13)/(sqrt(1.0*gam12Sq + gam13Sq)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    dL63dEps(0) =  -1.0000000000000002*epsilon11*gam12Sq*gam13Sq/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)) - 1.0000000000000002*epsilon11*pow(gam13Sq, 2)/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    // \gamma_{12}
    dL12dEps(3) =  -2.0*gam13Sq*gamma12/pow(1.0*gam12Sq + gam13Sq, 2) ;
    dL22dEps(3) = gamma12/rootE11Gam12Gam13*(gam12Sq*(0.5*epsilon11 - rootE11Gam12Gam13) + rootE11Gam12Gam13*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL32dEps(3) =  gamma12/rootE11Gam12Gam13*(-0.5*gam12Sq*(epsilon11 + 2.0*rootE11Gam12Gam13) + rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2) ;
    dL42dEps(3) =  -1.4142135623730951*gamma12*1.0/rootE11Gam12Gam13*(gam12Sq + gam13Sq)*(rootE11Gam12Gam13*(gam12Sq + gam13Sq)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma13/gamma12) - rootE11Gam12Gam13*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 - rootE11Gam12Gam13)*(gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq - epsilon11*gam12Sq*rootE11Gam12Gam13 - epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL52dEps(3) =  1.4142135623730951*gamma12*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(-rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma13/gamma12) + rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq + 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 + epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL62dEps(3) =  gamma12*pow(rootE11Gam12Gam13, -2.0)*(1.0000000000000002*gam12Sq*rootE11Gam12Gam13*(0.5*epsilon11 - 1.0*rootE11Gam12Gam13)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) - 1.0000000000000002*gam12Sq*rootE11Gam12Gam13*(0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) + 2.0000000000000004*pow(rootE11Gam12Gam13, 2.0)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));

    dL13dEps(3) =  2.0*gam13Sq*gamma12/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam13Sq + 1.0*pow(gam13Sq, 2)) ;
    dL23dEps(3) =  gam13Sq*gamma12*1.0/rootE11Gam12Gam13*(0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL33dEps(3) =  gam13Sq*gamma12*1.0/rootE11Gam12Gam13*(-0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL43dEps(3) =  1.4142135623730951*gamma12*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma13/gamma12) - rootE11Gam12Gam13*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 - 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq - 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 - epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL53dEps(3) =  -1.4142135623730951*gamma12*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(-rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma13/gamma12) + rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq + 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 + epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL63dEps(3) =  gam13Sq*gamma12*1.0/rootE11Gam12Gam13*((0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) - (0.50000000000000011*epsilon11 + 1.0000000000000002*rootE11Gam12Gam13)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));
    // \gamma_{13}
    dL12dEps(4) =  2.0*gam12Sq*gamma13/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam13Sq + 1.0*pow(gam13Sq, 2));
    dL22dEps(4) =  gam12Sq*gamma13*1.0/rootE11Gam12Gam13*(0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL32dEps(4) =  gam12Sq*gamma13*1.0/rootE11Gam12Gam13*(-0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL42dEps(4) =  -1.4142135623730951*gamma13*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12/gamma13) - rootE11Gam12Gam13*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 - 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq - 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 - epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL52dEps(4) =  1.4142135623730951*gamma13*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(-rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12/gamma13) + rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq + 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 + epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL62dEps(4) =  gam12Sq*gamma13*1.0/rootE11Gam12Gam13*((0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) - (0.50000000000000011*epsilon11 + 1.0000000000000002*rootE11Gam12Gam13)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));

    dL13dEps(4) =  -2.0*gam12Sq*gamma13/pow(1.0*gam12Sq + gam13Sq, 2);
    dL23dEps(4) =  gamma13*1.0/rootE11Gam12Gam13*(gam13Sq*(0.50000000000000011*epsilon11 - 1.0000000000000002*rootE11Gam12Gam13) + rootE11Gam12Gam13*(0.50000000000000011*eps11Sq - 1.0000000000000002*epsilon11*rootE11Gam12Gam13 + 0.50000000000000011*gam12Sq + 0.50000000000000011*gam13Sq + 0.50000000000000011*pow(rootE11Gam12Gam13, 2.0)))/pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL33dEps(4) =  gamma13*1.0/rootE11Gam12Gam13*(-0.50000000000000011*gam13Sq*(1.0*epsilon11 + 2.0*rootE11Gam12Gam13) + 1.0000000000000002*rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 2);
    dL43dEps(4) = 1.4142135623730951*gamma13*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12/gamma13) - rootE11Gam12Gam13*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 - 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq - 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 - epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL53dEps(4) = -1.4142135623730951*gamma13*1.0/rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(-rootE11Gam12Gam13*(1.0*gam12Sq + gam13Sq)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12/gamma13) + rootE11Gam12Gam13*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*fabs(gamma12*gamma13) + (0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))/pow(0.5*eps11Sq*gam12Sq + 0.5*eps11Sq*gam13Sq + 1.0*epsilon11*gam12Sq*rootE11Gam12Gam13 + epsilon11*gam13Sq*rootE11Gam12Gam13 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam13Sq + 0.5*gam12Sq*pow(rootE11Gam12Gam13, 2.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(rootE11Gam12Gam13, 2.0), 5.0/2.0);
    dL63dEps(4) =  gamma13*pow(rootE11Gam12Gam13, -2.0)*(1.0000000000000002*gam13Sq*rootE11Gam12Gam13*(0.5*epsilon11 - 1.0*rootE11Gam12Gam13)*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) - 1.0000000000000002*gam13Sq*rootE11Gam12Gam13*(0.5*epsilon11 + 1.0*rootE11Gam12Gam13)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)) + 2.0000000000000004*pow(rootE11Gam12Gam13, 2.0)*(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0))*(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0)))/(pow(0.5*eps11Sq - epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0)*pow(0.5*eps11Sq + epsilon11*rootE11Gam12Gam13 + 0.5*gam12Sq + 0.5*gam13Sq + 0.5*pow(rootE11Gam12Gam13, 2.0), 3.0/2.0));

    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(0) = -0.5*epsilon11*1.0/rootE11Gam12Gam13 + 0.5;
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE11Gam12Gam13;
      dEVal2dEps(4) = -0.5*gamma13*1.0/rootE11Gam12Gam13;
    }
    
    if (eValsD1(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(0) = 0.5*epsilon11*1.0/rootE11Gam12Gam13 + 0.5;
      dEVal3dEps(3) = 0.5*gamma12*1.0/rootE11Gam12Gam13;
      dEVal3dEps(4) = 0.5*gamma13*1.0/rootE11Gam12Gam13;
    }
  } else if (hasE11 && hasE12 && !hasE13) {
    // case 2 - E13 is zero
    // \epsilon_{11}
    dL22dEps(0) = gam12Sq*pow(eps11Sq + gam12Sq, -0.5)*(2.0*eps11Sq + sqrt(eps11Sq + gam12Sq)*(-4.0*epsilon11 + 2.0*sqrt(eps11Sq + gam12Sq)))/pow(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0), 2);
    dL32dEps(0) = gam12Sq*(-2.0*eps11Sq + sqrt(eps11Sq + gam12Sq)*(-4.0*epsilon11 - 2.0*sqrt(eps11Sq + gam12Sq)))*pow(eps11Sq + gam12Sq, -0.5)/pow(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0), 2);
    dL62dEps(0) = -8.0*epsilon11*pow(gam12Sq, 2)*sqrt(eps11Sq + gam12Sq)/(sqrt(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0))*sqrt(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0))*(1.0*pow(eps11Sq, 2)*sqrt(eps11Sq + gam12Sq) + 2.0*eps11Sq*gam12Sq*sqrt(eps11Sq + gam12Sq) - 2.0*eps11Sq*pow(eps11Sq + gam12Sq, 1.5) + 1.0*pow(gam12Sq, 2)*sqrt(eps11Sq + gam12Sq) + 2.0*gam12Sq*pow(eps11Sq + gam12Sq, 1.5) + 1.0*pow(eps11Sq + gam12Sq, 2.5)));
    // \gamma_{12}
    dL22dEps(3) = gamma12*pow(eps11Sq + gam12Sq, -0.5)*(1.0*gam12Sq*(2.0*epsilon11 - 4.0*sqrt(eps11Sq + gam12Sq)) + 2.0*sqrt(eps11Sq + gam12Sq)*(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0)))/pow(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0), 2);
    dL32dEps(3) = gamma12*pow(eps11Sq + gam12Sq, -0.5)*(-1.0*gam12Sq*(2.0*epsilon11 + 4.0*sqrt(eps11Sq + gam12Sq)) + 2.0*sqrt(eps11Sq + gam12Sq)*(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0)))/pow(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0), 2);
    dL62dEps(3) = gamma12*(4.0*pow(eps11Sq, 2)*pow(eps11Sq + gam12Sq, 1.0) + 8.0*eps11Sq*gam12Sq*pow(eps11Sq + gam12Sq, 1.0) - 8.0*eps11Sq*pow(eps11Sq + gam12Sq, 2.0) - 4.0*pow(gam12Sq, 2)*pow(eps11Sq + gam12Sq, 1.0) + 4.0*pow(eps11Sq + gam12Sq, 3.0))/(sqrt(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0))*sqrt(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam12Sq) + gam12Sq + pow(eps11Sq + gam12Sq, 1.0))*(1.0*pow(eps11Sq, 2)*pow(eps11Sq + gam12Sq, 1.0) + 2.0*eps11Sq*gam12Sq*pow(eps11Sq + gam12Sq, 1.0) - 2.0*eps11Sq*pow(eps11Sq + gam12Sq, 2.0) + 1.0*pow(gam12Sq, 2)*pow(eps11Sq + gam12Sq, 1.0) + 2.0*gam12Sq*pow(eps11Sq + gam12Sq, 2.0) + 1.0*pow(eps11Sq + gam12Sq, 3.0)));

    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(0) = -0.5*epsilon11*1.0/rootE11Gam12Gam13 + 0.5;
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE11Gam12Gam13;
    }
  } else if (hasE11 && !hasE12 && hasE13) {
    // case 3 - E12 is zero
    // \epsilon_{11}
    dL23dEps(0) = gam13Sq*pow(eps11Sq + gam13Sq, -0.5)*(2.0*eps11Sq + sqrt(eps11Sq + gam13Sq)*(-4.0*epsilon11 + 2.0*sqrt(eps11Sq + gam13Sq)))/pow(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0), 2);
    dL33dEps(0) = gam13Sq*(-2.0*eps11Sq + sqrt(eps11Sq + gam13Sq)*(-4.0*epsilon11 - 2.0*sqrt(eps11Sq + gam13Sq)))*pow(eps11Sq + gam13Sq, -0.5)/pow(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0), 2);
    dL63dEps(0) = -8.0*epsilon11*pow(gam13Sq, 2)*sqrt(eps11Sq + gam13Sq)/(sqrt(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0))*sqrt(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0))*(1.0*pow(eps11Sq, 2)*sqrt(eps11Sq + gam13Sq) + 2.0*eps11Sq*gam13Sq*sqrt(eps11Sq + gam13Sq) - 2.0*eps11Sq*pow(eps11Sq + gam13Sq, 1.5) + 1.0*pow(gam13Sq, 2)*sqrt(eps11Sq + gam13Sq) + 2.0*gam13Sq*pow(eps11Sq + gam13Sq, 1.5) + 1.0*pow(eps11Sq + gam13Sq, 2.5)));
    // \gamma_{13}
    dL23dEps(4) = gamma13*pow(eps11Sq + gam13Sq, -0.5)*(1.0*gam13Sq*(2.0*epsilon11 - 4.0*sqrt(eps11Sq + gam13Sq)) + 2.0*sqrt(eps11Sq + gam13Sq)*(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0)))/pow(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0), 2);
    dL33dEps(4) = gamma13*pow(eps11Sq + gam13Sq, -0.5)*(-1.0*gam13Sq*(2.0*epsilon11 + 4.0*sqrt(eps11Sq + gam13Sq)) + 2.0*sqrt(eps11Sq + gam13Sq)*(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0)))/pow(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0), 2);
    dL63dEps(4) = gamma13*(4.0*pow(eps11Sq, 2)*pow(eps11Sq + gam13Sq, 1.0) + 8.0*eps11Sq*gam13Sq*pow(eps11Sq + gam13Sq, 1.0) - 8.0*eps11Sq*pow(eps11Sq + gam13Sq, 2.0) - 4.0*pow(gam13Sq, 2)*pow(eps11Sq + gam13Sq, 1.0) + 4.0*pow(eps11Sq + gam13Sq, 3.0))/(sqrt(eps11Sq - 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0))*sqrt(eps11Sq + 2*epsilon11*sqrt(eps11Sq + gam13Sq) + gam13Sq + pow(eps11Sq + gam13Sq, 1.0))*(1.0*pow(eps11Sq, 2)*pow(eps11Sq + gam13Sq, 1.0) + 2.0*eps11Sq*gam13Sq*pow(eps11Sq + gam13Sq, 1.0) - 2.0*eps11Sq*pow(eps11Sq + gam13Sq, 2.0) + 1.0*pow(gam13Sq, 2)*pow(eps11Sq + gam13Sq, 1.0) + 2.0*gam13Sq*pow(eps11Sq + gam13Sq, 2.0) + 1.0*pow(eps11Sq + gam13Sq, 3.0)));

    // eVal3
    if (eValsD1(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(0) = 0.5*epsilon11*1.0/rootE11Gam12Gam13 + 0.5;
      dEVal3dEps(4) = 0.5*gamma13*1.0/rootE11Gam12Gam13;
    }


  } else if (!hasE11 && hasE12 && hasE13) {
    // case 4 - E11 is zero
    // \gamma_{12}
    dL12dEps(3) = -2.0*gam13Sq*gamma12/pow(1.0*gam12Sq + gam13Sq, 2) ;
    dL22dEps(3) = gamma12*(-4.0*gam12Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 2*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);
    dL32dEps(3) = gamma12*(-4.0*gam12Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 2*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);
    dL42dEps(3) = 0.17677669529663689*gamma12*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma13/gamma12) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL52dEps(3) = 0.17677669529663689*gamma12*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma13/gamma12) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL62dEps(3) = gamma12*(-8.0*gam12Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 4.0*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);
    dL13dEps(3) = 2.0*gam13Sq*gamma12/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam13Sq + 1.0*pow(gam13Sq, 2)) ;
    dL23dEps(3) = -4.0*gam13Sq*gamma12/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2) ;
    dL33dEps(3) = -4.0*gam13Sq*gamma12/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2) ;
    dL43dEps(3) = -0.17677669529663689*gamma12*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma13/gamma12) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL53dEps(3) = -0.17677669529663689*gamma12*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma13/gamma12) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL63dEps(3) = -8.0*gam13Sq*gamma12/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2);
    // \gamma_{13}
    dL12dEps(4) = 2.0*gam12Sq*gamma13/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam13Sq + 1.0*pow(gam13Sq, 2));
    dL22dEps(4) = -4.0*gam12Sq*gamma13/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2);
    dL32dEps(4) = -4.0*gam12Sq*gamma13/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2);
    dL42dEps(4) = 0.17677669529663689*gamma13*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12/gamma13) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL52dEps(4) = 0.17677669529663689*gamma13*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12/gamma13) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL62dEps(4) = -8.0*gam12Sq*gamma13/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2);
    dL13dEps(4) = -2.0*gam12Sq*gamma13/pow(1.0*gam12Sq + gam13Sq, 2);
    dL23dEps(4) = gamma13*(-4.0*gam13Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 2.0*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);
    dL33dEps(4) = gamma13*(-4.0*gam13Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 2.0*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);
    dL43dEps(4) = -0.17677669529663689*gamma13*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12/gamma13) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL53dEps(4) = -0.17677669529663689*gamma13*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam13Sq)*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12/gamma13) + 4.0*(1.0*gam12Sq + gam13Sq)*fabs(gamma12*gamma13) + 2.0*(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0))*fabs(gamma12*gamma13))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam13Sq + 0.5*gam12Sq*pow(gam12Sq + gam13Sq, 1.0) + 0.5*pow(gam13Sq, 2) + 0.5*gam13Sq*pow(gam12Sq + gam13Sq, 1.0), 5.0/2.0);
    dL63dEps(4) = gamma13*(-8.0*gam13Sq*(gam12Sq + gam13Sq + pow(gam12Sq + gam13Sq, 1.0)) + 4.0*pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam13Sq + pow(gam12Sq + gam13Sq, 1.0), 3);

    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE11Gam12Gam13;
      dEVal2dEps(4) = -0.5*gamma13*1.0/rootE11Gam12Gam13;
    }
    
    if (eValsD1(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(3) = 0.5*gamma12*1.0/rootE11Gam12Gam13;
      dEVal3dEps(4) = 0.5*gamma13*1.0/rootE11Gam12Gam13;
    }

  } else if (hasE11 && !hasE12 && !hasE13) {
    // case 5 - E12 and E13 are zero
    // in this case we are in spectral form so L is constant
    // 
    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(0) = 1.0;
    }
 
  } else if (!hasE11 && hasE12 && !hasE13) {
    // case 6 - E11 and E13 are zero
    // Here again L is constant
    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(3) = -0.5;
    }
    
    if (eValsD1(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(3) = 0.5;
    }
  } else if (!hasE11 && !hasE12 && hasE13) {
    // case 7 - E11 and E12 are zero
    // Here again L is constant
    // eVal2
    if (eValsD1(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(4) = -0.5;
    }
    
    if (eValsD1(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D1}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(4) = 0.5;
    }
  }

  // Full derivs
  // Row 1
  // \epsilon_{11}
  dEpsD1PlusDEps(0,0) = dL12dEps(0) * eValsD1(1,1) + dL13dEps(0) * eValsD1(2,2) +
                         L12 * dEVal2dEps(0) + L13 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(0,3) = dL12dEps(3) * eValsD1(1,1) + dL13dEps(3) * eValsD1(2,2) +
                         L12 * dEVal2dEps(3) + L13 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(0,4) = dL12dEps(4) * eValsD1(1,1) + dL13dEps(4) * eValsD1(2,2) +
                         L12 * dEVal2dEps(4) + L13 * dEVal3dEps(4);
  // Row 2
  // \epsilon_{11}
  dEpsD1PlusDEps(1,0) = dL22dEps(0) * eValsD1(1,1) + dL23dEps(0) * eValsD1(2,2) +
                         L22 * dEVal2dEps(0) + L23 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(1,3) = dL22dEps(3) * eValsD1(1, 1) + dL23dEps(3) * eValsD1(2, 2) +
                         L22 * dEVal2dEps(3) + L23 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(1,4) = dL22dEps(4) * eValsD1(1, 1) + dL23dEps(4) * eValsD1(2, 2) +
                         L22 * dEVal2dEps(4) + L23 * dEVal3dEps(4);
  // Row 3
  // \epsilon_{11}
  dEpsD1PlusDEps(2,0) = dL32dEps(0) * eValsD1(1, 1) + dL33dEps(0) * eValsD1(2, 2) +
                         L32 * dEVal2dEps(0) + L33 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(2,3) = dL32dEps(3) * eValsD1(1, 1) + dL33dEps(3) * eValsD1(2, 2) +
                         L32 * dEVal2dEps(3) + L33 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(2,4) = dL32dEps(4) * eValsD1(1, 1) + dL33dEps(4) * eValsD1(2, 2) +
                         L32 * dEVal2dEps(4) + L33 * dEVal3dEps(4);
  // Row 4
  // \epsilon_{11}
  dEpsD1PlusDEps(3,0) = dL42dEps(0) * eValsD1(1, 1) + dL43dEps(0) * eValsD1(2, 2) +
                         L42 * dEVal2dEps(0) + L43 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(3,3) = dL42dEps(3) * eValsD1(1, 1) + dL43dEps(3) * eValsD1(2, 2) +
                         L42 * dEVal2dEps(3) + L43 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(3,4) = dL42dEps(4) * eValsD1(1, 1) + dL43dEps(4) * eValsD1(2, 2) +
                         L42 * dEVal2dEps(4) + L43 * dEVal3dEps(4);
  // Row 5
  // \epsilon_{11}
  dEpsD1PlusDEps(4,0) = dL52dEps(0) * eValsD1(1, 1) + dL53dEps(0) * eValsD1(2, 2) +
                         L52 * dEVal2dEps(0) + L53 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(4,3) = dL52dEps(3) * eValsD1(1, 1) + dL53dEps(3) * eValsD1(2, 2) +
                         L52 * dEVal2dEps(3) + L53 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(4,4) = dL52dEps(4) * eValsD1(1, 1) + dL53dEps(4) * eValsD1(2, 2) +
                         L52 * dEVal2dEps(4) + L53 * dEVal3dEps(4);
  // Row 6
  // \epsilon_{11}
  dEpsD1PlusDEps(5,0) = dL62dEps(0) * eValsD1(1, 1) + dL63dEps(0) * eValsD1(2, 2) +
                         L62 * dEVal2dEps(0) + L63 * dEVal3dEps(0);
  // \gamma_{12}
  dEpsD1PlusDEps(5,3) = dL62dEps(3) * eValsD1(1, 1) + dL63dEps(3) * eValsD1(2, 2) +
                         L62 * dEVal2dEps(3) + L63 * dEVal3dEps(3);
  // \gamma_{13}
  dEpsD1PlusDEps(5,4) = dL62dEps(4) * eValsD1(1,1) + dL63dEps(4) * eValsD1(2, 2) +
                         L62 * dEVal2dEps(4) + L63 * dEVal3dEps(4);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::calcDEpsD2PlusDEps(const Vector6d& epsD2,
                                           Matrix6d& dEpsD2PlusDEps) const {
  // This funciton calcs derivative of epsilonD1+ wrt epsilon
  // It uses a number of simplifications, but should be sufficiently capable
  // ASSUMPTIONS:
  // eValsD2(0,0) = 0.0
  // eValsD2(i,j) = 0.0 for all i != j
  //
  // following this basic setup:
  // epsilonD2+ = L_{ij} <\hat{\epsilon}_{j}^{D2}>_+

  // get strains
  const double epsilon22 = epsD2(1), gamma12 = epsD2(3), gamma23 = epsD2(5);

  // pick given strain state
  const bool hasE22 = (epsilon22 != 0.0), hasE12 = (gamma12 != 0.0),
    hasE23 = (gamma23 != 0.0);

  // Transformation matrix spectral -> Global
  // since I know I have limited \hat{\epsilon_i^{D2+}}
  // I only need part of this
  double L12 = eVectsD2(1,0)*eVectsD2(1,0), L13 = eVectsD2(2,0)*eVectsD2(2,0),
    L22 = eVectsD2(1,1)*eVectsD2(1,1), L23 = eVectsD2(2,1)*eVectsD2(2,1),
    L32 = eVectsD2(1,2)*eVectsD2(1,2), L33 = eVectsD2(2,2)*eVectsD2(2,2),
    L42 = 2.0*eVectsD2(1,0)*eVectsD2(1,1), L43 = 2.0*eVectsD2(2,0)*eVectsD2(2,1),
    L52 = 2.0*eVectsD2(1,0)*eVectsD2(1,2), L53 = 2.0*eVectsD2(2,0)*eVectsD2(2,2),
    L62 = 2.0*eVectsD2(1,1)*eVectsD2(1,2), L63 = 2.0*eVectsD2(2,1)*eVectsD2(2,2);
  // initialize d L_{ij}/d \epsilon_{k}
  Vector6d dL12dEps = Vector6d::Zero(), dL13dEps = Vector6d::Zero(), 
    dL22dEps = Vector6d::Zero(), dL23dEps = Vector6d::Zero(), 
    dL32dEps = Vector6d::Zero(), dL33dEps = Vector6d::Zero(), 
    dL42dEps = Vector6d::Zero(), dL43dEps = Vector6d::Zero(), 
    dL52dEps = Vector6d::Zero(), dL53dEps = Vector6d::Zero(), 
    dL62dEps = Vector6d::Zero(), dL63dEps = Vector6d::Zero();

  // derivative of eVals
  Vector6d dEVal2dEps = Vector6d::Zero(), dEVal3dEps = Vector6d::Zero();

  // useful terms
  const double eps22Sq = epsilon22*epsilon22, gam12Sq = gamma12*gamma12,
    gam23Sq = gamma23*gamma23;
  const double rootE22Gam12Gam23 =
    std::sqrt(eps22Sq + gam12Sq + gam23Sq);
  if (hasE22 && hasE12 && hasE23) {
    // case 1 - all D2 strains present
    // d L_{ij}/d \epsilon_k for this case
    // \epsilon_{11}
    dL22dEps(1) = pow(rootE22Gam12Gam23, -2.0)*(eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*eps22Sq + rootE22Gam12Gam23*(-1.0000000000000002*epsilon22 + 0.50000000000000011*rootE22Gam12Gam23)) - 1.0000000000000002*eps22Sq*rootE22Gam12Gam23*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + pow(rootE22Gam12Gam23, 2.0)*(2.0000000000000004*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + pow(rootE22Gam12Gam23, 2.0)*(-1.0000000000000002*epsilon22*(1.0*eps22Sq + rootE22Gam12Gam23*(-2.0*epsilon22 + rootE22Gam12Gam23)) + rootE22Gam12Gam23*(0.50000000000000011*eps22Sq + rootE22Gam12Gam23*(-1.0000000000000002*epsilon22 + 0.50000000000000011*rootE22Gam12Gam23))))/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2) ;
    dL32dEps(1) = pow(rootE22Gam12Gam23, -2.0)*(-eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*eps22Sq + rootE22Gam12Gam23*(1.0000000000000002*epsilon22 + 0.50000000000000011*rootE22Gam12Gam23)) + 1.0000000000000002*eps22Sq*rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + pow(rootE22Gam12Gam23, 2.0)*(2.0000000000000004*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) - pow(rootE22Gam12Gam23, 2.0)*(1.0000000000000002*epsilon22*(1.0*eps22Sq + rootE22Gam12Gam23*(2.0*epsilon22 + rootE22Gam12Gam23)) + rootE22Gam12Gam23*(0.50000000000000011*eps22Sq + rootE22Gam12Gam23*(1.0000000000000002*epsilon22 + 0.50000000000000011*rootE22Gam12Gam23))))/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2) ;
    dL62dEps(1) = 1.0*epsilon22*(-eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) - eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) + gam12Sq*pow(rootE22Gam12Gam23, 4.0) + gam23Sq*pow(rootE22Gam12Gam23, 4.0))/(sqrt(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*sqrt(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*(0.25*pow(eps22Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*eps22Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam12Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam23Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(rootE22Gam12Gam23, 6.0))) ;
    dL23dEps(1) = 0.50000000000000011*gam23Sq*1.0/rootE22Gam12Gam23*(eps22Sq + rootE22Gam12Gam23*(-2.0*epsilon22 + rootE22Gam12Gam23))/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2) ;
    dL33dEps(1) = 0.50000000000000011*gam23Sq*1.0/rootE22Gam12Gam23*(-eps22Sq - rootE22Gam12Gam23*(2.0*epsilon22 + rootE22Gam12Gam23))/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2) ;
    dL43dEps(1) = 0.70710678118654757*1.0/rootE22Gam12Gam23*(eps22Sq + rootE22Gam12Gam23*(-2.0*epsilon22 + rootE22Gam12Gam23))*fabs(gamma12*gamma23)/(sqrt(1.0*gam12Sq + gam23Sq)*pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)) ;
    dL53dEps(1) = -0.70710678118654757*1.0/rootE22Gam12Gam23*(eps22Sq + rootE22Gam12Gam23*(2.0*epsilon22 + rootE22Gam12Gam23))*fabs(gamma12*gamma23)/(sqrt(1.0*gam12Sq + gam23Sq)*pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)) ;
    dL63dEps(1) = -1.0000000000000002*epsilon22*gam12Sq*gam23Sq/(pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)*pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)) - 1.0000000000000002*epsilon22*pow(gam23Sq, 2)/(pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)*pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)) ;
    // \gamma_{12}
    dL22dEps(3) = gamma12*pow(rootE22Gam12Gam23, -2.0)*(eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23) - 1.0000000000000002*epsilon22*rootE22Gam12Gam23*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + pow(rootE22Gam12Gam23, 2.0)*(-1.0000000000000002*epsilon22*(1.0*epsilon22 - 2.0*rootE22Gam12Gam23) + rootE22Gam12Gam23*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)) + 1.0000000000000002*pow(rootE22Gam12Gam23, 2.0)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL32dEps(3) = gamma12*pow(rootE22Gam12Gam23, -2.0)*(-eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23) + 1.0000000000000002*epsilon22*rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) - pow(rootE22Gam12Gam23, 2.0)*(1.0000000000000002*epsilon22*(1.0*epsilon22 + 2.0*rootE22Gam12Gam23) + rootE22Gam12Gam23*(0.50000000000000011*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23)) + 1.0000000000000002*pow(rootE22Gam12Gam23, 2.0)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL62dEps(3) = gamma12*(-0.5*pow(eps22Sq, 2)*pow(rootE22Gam12Gam23, 2.0) - 2.0*eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) - 2.0*eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*pow(gam12Sq, 2)*pow(rootE22Gam12Gam23, 2.0) - 1.0*gam12Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*pow(gam23Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(rootE22Gam12Gam23, 6.0))/(sqrt(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*sqrt(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*(0.25*pow(eps22Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*eps22Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam12Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam23Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(rootE22Gam12Gam23, 6.0)));
    dL13dEps(3) = 2.0*gam23Sq*gamma12/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam23Sq + 1.0*pow(gam23Sq, 2));
    dL23dEps(3) = gam23Sq*gamma12*1.0/rootE22Gam12Gam23*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL33dEps(3) = gam23Sq*gamma12*1.0/rootE22Gam12Gam23*(-0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2) ;
    dL43dEps(3) = 1.4142135623730951*gamma12*1.0/rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma23/gamma12) - rootE22Gam12Gam23*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12*gamma23) + (0.5*epsilon22 - 1.0*rootE22Gam12Gam23)*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23))*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))/pow(0.5*eps22Sq*gam12Sq + 0.5*eps22Sq*gam23Sq - 1.0*epsilon22*gam12Sq*rootE22Gam12Gam23 - epsilon22*gam23Sq*rootE22Gam12Gam23 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam23Sq + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 2.0), 5.0/2.0);
    dL53dEps(3) = -1.4142135623730951*gamma12*1.0/rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(-rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma23/gamma12) + rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12*gamma23) + (0.5*epsilon22 + 1.0*rootE22Gam12Gam23)*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23))*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))/pow(0.5*eps22Sq*gam12Sq + 0.5*eps22Sq*gam23Sq + 1.0*epsilon22*gam12Sq*rootE22Gam12Gam23 + epsilon22*gam23Sq*rootE22Gam12Gam23 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam23Sq + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 2.0), 5.0/2.0);
    dL63dEps(3) = gam23Sq*gamma12*1.0/rootE22Gam12Gam23*((0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) - (0.50000000000000011*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/(pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)*pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0));
    // \gamma_{23}
    dL22dEps(5) = gamma23*pow(rootE22Gam12Gam23, -2.0)*(eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23) - 1.0000000000000002*epsilon22*rootE22Gam12Gam23*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + pow(rootE22Gam12Gam23, 2.0)*(-1.0000000000000002*epsilon22*(1.0*epsilon22 - 2.0*rootE22Gam12Gam23) + rootE22Gam12Gam23*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23)) + 1.0000000000000002*pow(rootE22Gam12Gam23, 2.0)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL32dEps(5) = gamma23*pow(rootE22Gam12Gam23, -2.0)*(-eps22Sq*rootE22Gam12Gam23*(0.50000000000000011*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23) + 1.0000000000000002*epsilon22*rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) - pow(rootE22Gam12Gam23, 2.0)*(1.0000000000000002*epsilon22*(1.0*epsilon22 + 2.0*rootE22Gam12Gam23) + rootE22Gam12Gam23*(0.50000000000000011*epsilon22 + 1.0000000000000002*rootE22Gam12Gam23)) + 1.0000000000000002*pow(rootE22Gam12Gam23, 2.0)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL62dEps(5) = gamma23*(-0.5*pow(eps22Sq, 2)*pow(rootE22Gam12Gam23, 2.0) - 2.0*eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) - 2.0*eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*pow(gam12Sq, 2)*pow(rootE22Gam12Gam23, 2.0) - 1.0*gam12Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*pow(gam23Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(rootE22Gam12Gam23, 6.0))/(sqrt(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*sqrt(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*(0.25*pow(eps22Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*eps22Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) - 0.5*eps22Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam12Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*gam23Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(gam23Sq, 2)*pow(rootE22Gam12Gam23, 2.0) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 4.0) + 0.25*pow(rootE22Gam12Gam23, 6.0)));
    dL13dEps(5) = -2.0*gam12Sq*gamma23/pow(1.0*gam12Sq + gam23Sq, 2);
    dL23dEps(5) = gamma23*1.0/rootE22Gam12Gam23*(gam23Sq*(0.50000000000000011*epsilon22 - 1.0000000000000002*rootE22Gam12Gam23) + rootE22Gam12Gam23*(0.50000000000000011*eps22Sq - 1.0000000000000002*epsilon22*rootE22Gam12Gam23 + 0.50000000000000011*gam12Sq + 0.50000000000000011*gam23Sq + 0.50000000000000011*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL33dEps(5) = gamma23*1.0/rootE22Gam12Gam23*(-0.50000000000000011*gam23Sq*(1.0*epsilon22 + 2.0*rootE22Gam12Gam23) + 1.0000000000000002*rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 2);
    dL43dEps(5) = 1.4142135623730951*gamma23*1.0/rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12/gamma23) - rootE22Gam12Gam23*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12*gamma23) + (0.5*epsilon22 - 1.0*rootE22Gam12Gam23)*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23))*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))/pow(0.5*eps22Sq*gam12Sq + 0.5*eps22Sq*gam23Sq - 1.0*epsilon22*gam12Sq*rootE22Gam12Gam23 - epsilon22*gam23Sq*rootE22Gam12Gam23 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam23Sq + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 2.0), 5.0/2.0);
    dL53dEps(5) = -1.4142135623730951*gamma23*1.0/rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(-rootE22Gam12Gam23*(1.0*gam12Sq + gam23Sq)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12/gamma23) + rootE22Gam12Gam23*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*fabs(gamma12*gamma23) + (0.5*epsilon22 + 1.0*rootE22Gam12Gam23)*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23))*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))/pow(0.5*eps22Sq*gam12Sq + 0.5*eps22Sq*gam23Sq + 1.0*epsilon22*gam12Sq*rootE22Gam12Gam23 + epsilon22*gam23Sq*rootE22Gam12Gam23 + 0.5*pow(gam12Sq, 2) + 1.0*gam12Sq*gam23Sq + 0.5*gam12Sq*pow(rootE22Gam12Gam23, 2.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(rootE22Gam12Gam23, 2.0), 5.0/2.0);
    dL63dEps(5) = gamma23*pow(rootE22Gam12Gam23, -2.0)*(1.0000000000000002*gam23Sq*rootE22Gam12Gam23*(0.5*epsilon22 - 1.0*rootE22Gam12Gam23)*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) - 1.0000000000000002*gam23Sq*rootE22Gam12Gam23*(0.5*epsilon22 + 1.0*rootE22Gam12Gam23)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)) + 2.0000000000000004*pow(rootE22Gam12Gam23, 2.0)*(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0))*(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0)))/(pow(0.5*eps22Sq - epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0)*pow(0.5*eps22Sq + epsilon22*rootE22Gam12Gam23 + 0.5*gam12Sq + 0.5*gam23Sq + 0.5*pow(rootE22Gam12Gam23, 2.0), 3.0/2.0));

    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(1) = -0.5*epsilon22*1.0/rootE22Gam12Gam23 + 0.5;
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE22Gam12Gam23;
      dEVal2dEps(5) = -0.5*gamma23*1.0/rootE22Gam12Gam23;
    }
    
    if (eValsD2(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(1) = 0.5*epsilon22*1.0/rootE22Gam12Gam23 + 0.5;
      dEVal3dEps(3) = 0.5*gamma12*1.0/rootE22Gam12Gam23;
      dEVal3dEps(5) = 0.5*gamma23*1.0/rootE22Gam12Gam23;
    }
  } else if (hasE22 && hasE12 && !hasE23) {
    // case 2 - E23 is zero
    // \epsilon_{22}
    dL22dEps(1) = gam12Sq*(-2.0*eps22Sq + sqrt(eps22Sq + gam12Sq)*(-4.0*epsilon22 - 2.0*sqrt(eps22Sq + gam12Sq)))*pow(eps22Sq + gam12Sq, -0.5)/pow(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0), 2);
    dL32dEps(1) = gam12Sq*pow(eps22Sq + gam12Sq, -0.5)*(2.0*eps22Sq + sqrt(eps22Sq + gam12Sq)*(-4.0*epsilon22 + 2.0*sqrt(eps22Sq + gam12Sq)))/pow(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0), 2);
    dL62dEps(1) = -8.0*epsilon22*pow(gam12Sq, 2)*sqrt(eps22Sq + gam12Sq)/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0))*(1.0*pow(eps22Sq, 2)*sqrt(eps22Sq + gam12Sq) + 2.0*eps22Sq*gam12Sq*sqrt(eps22Sq + gam12Sq) - 2.0*eps22Sq*pow(eps22Sq + gam12Sq, 1.5) + 1.0*pow(gam12Sq, 2)*sqrt(eps22Sq + gam12Sq) + 2.0*gam12Sq*pow(eps22Sq + gam12Sq, 1.5) + 1.0*pow(eps22Sq + gam12Sq, 2.5)));

    // \gamma_{12}
    dL22dEps(3) = gamma12*pow(eps22Sq + gam12Sq, -0.5)*(-1.0*gam12Sq*(2.0*epsilon22 + 4.0*sqrt(eps22Sq + gam12Sq)) + 2.0*sqrt(eps22Sq + gam12Sq)*(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0)))/pow(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0), 2);
    dL32dEps(3) = gamma12*pow(eps22Sq + gam12Sq, -0.5)*(1.0*gam12Sq*(2.0*epsilon22 - 4.0*sqrt(eps22Sq + gam12Sq)) + 2.0*sqrt(eps22Sq + gam12Sq)*(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0)))/pow(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0), 2);
    dL62dEps(3) = gamma12*(4.0*pow(eps22Sq, 2)*pow(eps22Sq + gam12Sq, 1.0) + 8.0*eps22Sq*gam12Sq*pow(eps22Sq + gam12Sq, 1.0) - 8.0*eps22Sq*pow(eps22Sq + gam12Sq, 2.0) - 4.0*pow(gam12Sq, 2)*pow(eps22Sq + gam12Sq, 1.0) + 4.0*pow(eps22Sq + gam12Sq, 3.0))/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam12Sq) + gam12Sq + pow(eps22Sq + gam12Sq, 1.0))*(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam12Sq, 1.0) + 2.0*eps22Sq*gam12Sq*pow(eps22Sq + gam12Sq, 1.0) - 2.0*eps22Sq*pow(eps22Sq + gam12Sq, 2.0) + 1.0*pow(gam12Sq, 2)*pow(eps22Sq + gam12Sq, 1.0) + 2.0*gam12Sq*pow(eps22Sq + gam12Sq, 2.0) + 1.0*pow(eps22Sq + gam12Sq, 3.0)));

    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(1) = -0.5*epsilon22*1.0/rootE22Gam12Gam23 + 0.5;
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE22Gam12Gam23;
    }
  } else if (hasE22 && !hasE12 && hasE23) {
    // case 3 - E12 is zero
    // \epsilon_{22}
    dL22dEps(1) = gam23Sq*(-2.0*eps22Sq*sqrt(eps22Sq + gam23Sq) + 4.0*epsilon22*pow(eps22Sq + gam23Sq, 1.0) - 2.0*pow(eps22Sq + gam23Sq, 1.5))/(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) + 6.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) - 4.0*pow(epsilon22, 3)*pow(eps22Sq + gam23Sq, 1.5) - 4.0*epsilon22*gam23Sq*pow(eps22Sq + gam23Sq, 1.5) - 4.0*epsilon22*pow(eps22Sq + gam23Sq, 2.5) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0));
    dL32dEps(1) = gam23Sq*(2.0*eps22Sq*sqrt(eps22Sq + gam23Sq) + 4.0*epsilon22*pow(eps22Sq + gam23Sq, 1.0) + 2.0*pow(eps22Sq + gam23Sq, 1.5))/(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) + 6.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) + 4.0*pow(epsilon22, 3)*pow(eps22Sq + gam23Sq, 1.5) + 4.0*epsilon22*gam23Sq*pow(eps22Sq + gam23Sq, 1.5) + 4.0*epsilon22*pow(eps22Sq + gam23Sq, 2.5) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0));
    dL62dEps(1) = 8.0*epsilon22*gam23Sq*(-eps22Sq*pow(eps22Sq + gam23Sq, 1.0) + pow(eps22Sq + gam23Sq, 2.0))/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 2.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0)));
    dL23dEps(1) = gam23Sq*pow(eps22Sq + gam23Sq, -0.5)*(2.0*eps22Sq + sqrt(eps22Sq + gam23Sq)*(-4.0*epsilon22 + 2.0*sqrt(eps22Sq + gam23Sq)))/pow(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0), 2) ;
    dL33dEps(1) = gam23Sq*(-2.0*eps22Sq + sqrt(eps22Sq + gam23Sq)*(-4.0*epsilon22 - 2.0*sqrt(eps22Sq + gam23Sq)))*pow(eps22Sq + gam23Sq, -0.5)/pow(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0), 2);
    dL63dEps(1) = -8.0*epsilon22*pow(gam23Sq, 2)*sqrt(eps22Sq + gam23Sq)/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*(1.0*pow(eps22Sq, 2)*sqrt(eps22Sq + gam23Sq) + 2.0*eps22Sq*gam23Sq*sqrt(eps22Sq + gam23Sq) - 2.0*eps22Sq*pow(eps22Sq + gam23Sq, 1.5) + 1.0*pow(gam23Sq, 2)*sqrt(eps22Sq + gam23Sq) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 1.5) + 1.0*pow(eps22Sq + gam23Sq, 2.5)));

    // \gamma_{23}
    dL22dEps(5) = gamma23*1.0/(eps22Sq + gam23Sq)*(eps22Sq*sqrt(eps22Sq + gam23Sq)*(2.0*epsilon22 - 4.0*sqrt(eps22Sq + gam23Sq)) - 2.0*epsilon22*sqrt(eps22Sq + gam23Sq)*(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0)) + pow(eps22Sq + gam23Sq, 1.0)*(-2.0*epsilon22 + sqrt(eps22Sq + gam23Sq))*(2.0*epsilon22 - 4.0*sqrt(eps22Sq + gam23Sq)) + 2.0*pow(eps22Sq + gam23Sq, 1.0)*(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0)))/pow(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0), 2);
    dL32dEps(5) = gamma23*(-2.0*eps22Sq*pow(eps22Sq + gam23Sq, 1.0) + 2.0*epsilon22*gam23Sq*sqrt(eps22Sq + gam23Sq) - 4.0*epsilon22*pow(eps22Sq + gam23Sq, 1.5) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 2.0*pow(eps22Sq + gam23Sq, 2.0))/(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) + 6.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) + 4.0*pow(epsilon22, 3)*pow(eps22Sq + gam23Sq, 1.5) + 4.0*epsilon22*gam23Sq*pow(eps22Sq + gam23Sq, 1.5) + 4.0*epsilon22*pow(eps22Sq + gam23Sq, 2.5) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0));
    dL62dEps(5) = gamma23*(-4.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) - 16.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 4.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 4.0*pow(eps22Sq + gam23Sq, 3.0))/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 2.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0)));
    dL23dEps(5) = gamma23*pow(eps22Sq + gam23Sq, -0.5)*(1.0*gam23Sq*(2.0*epsilon22 - 4.0*sqrt(eps22Sq + gam23Sq)) + 2.0*sqrt(eps22Sq + gam23Sq)*(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0)))/pow(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0), 2);
    dL33dEps(5) = gamma23*pow(eps22Sq + gam23Sq, -0.5)*(-1.0*gam23Sq*(2.0*epsilon22 + 4.0*sqrt(eps22Sq + gam23Sq)) + 2.0*sqrt(eps22Sq + gam23Sq)*(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0)))/pow(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0), 2);
    dL63dEps(5) = gamma23*(4.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 8.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 8.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) - 4.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 4.0*pow(eps22Sq + gam23Sq, 3.0))/(sqrt(eps22Sq - 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*sqrt(eps22Sq + 2*epsilon22*sqrt(eps22Sq + gam23Sq) + gam23Sq + pow(eps22Sq + gam23Sq, 1.0))*(1.0*pow(eps22Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*eps22Sq*gam23Sq*pow(eps22Sq + gam23Sq, 1.0) - 2.0*eps22Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(gam23Sq, 2)*pow(eps22Sq + gam23Sq, 1.0) + 2.0*gam23Sq*pow(eps22Sq + gam23Sq, 2.0) + 1.0*pow(eps22Sq + gam23Sq, 3.0)));

    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(1) = -0.5*epsilon22*1.0/rootE22Gam12Gam23 + 0.5;
      dEVal3dEps(5) = -0.5*gamma23*1.0/rootE22Gam12Gam23;
    }
    // eVal3
    if (eValsD2(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(1) = 0.5*epsilon22*1.0/rootE22Gam12Gam23 + 0.5;
      dEVal3dEps(5) = 0.5*gamma23*1.0/rootE22Gam12Gam23;
    }


  } else if (!hasE22 && hasE12 && hasE23) {
    // case 4 - E22 is zero
    // \gamma_{12}
    dL22dEps(3) = 2.0*gamma12*(gam12Sq + gam23Sq - pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL32dEps(3) = 2.0*gamma12*(gam12Sq + gam23Sq - pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL62dEps(3) = 4.0*gamma12*(-gam12Sq - gam23Sq + pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL13dEps(3) = 2.0*gam23Sq*gamma12/(1.0*pow(gam12Sq, 2) + 2.0*gam12Sq*gam23Sq + 1.0*pow(gam23Sq, 2));
    dL23dEps(3) = -4.0*gam23Sq*gamma12/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL33dEps(3) = -4.0*gam23Sq*gamma12/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL43dEps(3) = -0.17677669529663689*gamma12*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma23/gamma12) + 4.0*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23) + 2.0*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12*gamma23))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam23Sq + 0.5*gam12Sq*pow(gam12Sq + gam23Sq, 1.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(gam12Sq + gam23Sq, 1.0), 5.0/2.0);
    dL53dEps(3) = -0.17677669529663689*gamma12*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma23/gamma12) + 4.0*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23) + 2.0*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12*gamma23))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam23Sq + 0.5*gam12Sq*pow(gam12Sq + gam23Sq, 1.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(gam12Sq + gam23Sq, 1.0), 5.0/2.0);
    dL63dEps(3) = -8.0*gam23Sq*gamma12/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    // \gamma_{23}
    dL22dEps(5) = 2.0*gamma23*(gam12Sq + gam23Sq - pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL32dEps(5) = 2.0*gamma23*(gam12Sq + gam23Sq - pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL62dEps(5) = 4.0*gamma23*(-gam12Sq - gam23Sq + pow(gam12Sq + gam23Sq, 1.0))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2);
    dL13dEps(5) = -2.0*gam12Sq*gamma23/pow(1.0*gam12Sq + gam23Sq, 2);
    dL23dEps(5) = gamma23*(-4.0*gam23Sq*(gam12Sq + gam23Sq + pow(gam12Sq + gam23Sq, 1.0)) + 2.0*pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 3);
    dL33dEps(5) = gamma23*(-4.0*gam23Sq*(gam12Sq + gam23Sq + pow(gam12Sq + gam23Sq, 1.0)) + 2.0*pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 3);
    dL43dEps(5) = -0.17677669529663689*gamma23*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12/gamma23) + 4.0*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23) + 2.0*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12*gamma23))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam23Sq + 0.5*gam12Sq*pow(gam12Sq + gam23Sq, 1.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(gam12Sq + gam23Sq, 1.0), 5.0/2.0);
    dL53dEps(5) = -0.17677669529663689*gamma23*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*(-2.0*(1.0*gam12Sq + gam23Sq)*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12/gamma23) + 4.0*(1.0*gam12Sq + gam23Sq)*fabs(gamma12*gamma23) + 2.0*(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0))*fabs(gamma12*gamma23))/pow(0.5*pow(gam12Sq, 2) + gam12Sq*gam23Sq + 0.5*gam12Sq*pow(gam12Sq + gam23Sq, 1.0) + 0.5*pow(gam23Sq, 2) + 0.5*gam23Sq*pow(gam12Sq + gam23Sq, 1.0), 5.0/2.0);
    dL63dEps(5) = gamma23*(-8.0*gam23Sq*(gam12Sq + gam23Sq + pow(gam12Sq + gam23Sq, 1.0)) + 4.0*pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 2))/pow(gam12Sq + 1.0*gam23Sq + pow(gam12Sq + gam23Sq, 1.0), 3);

    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(3) = -0.5*gamma12*1.0/rootE22Gam12Gam23;
      dEVal2dEps(5) = -0.5*gamma23*1.0/rootE22Gam12Gam23;
    }
    
    if (eValsD2(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(3) = 0.5*gamma12*1.0/rootE22Gam12Gam23;
      dEVal3dEps(5) = 0.5*gamma23*1.0/rootE22Gam12Gam23;
    }

  } else if (hasE22 && !hasE12 && !hasE23) {
    // case 5 - E12 and E23 are zero
    // in this case we are in spectral form so L is constant
    // 
    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(1) = 1.0;
    }
 
  } else if (!hasE22 && hasE12 && !hasE23) {
    // case 6 - E22 and E23 are zero
    // Here again L is constant
    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(3) = -0.5;
    }
    
    if (eValsD2(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(3) = 0.5;
    }
  } else if (!hasE22 && !hasE12 && hasE23) {
    // case 7 - E22 and E12 are zero
    // Here again L is constant
    // eVal2
    if (eValsD2(1,1) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal2dEps(5) = -0.5;
    }
    
    if (eValsD2(2,2) != 0.0) {
      // if this is false
      // d<\hat{\epsilon_{j}^D2}>+ / d \hat{\epsilon_{j}} == 0.0
      dEVal3dEps(5) = 0.5;
    }
  }

  // Full derivs
  // Row 1
  // \epsilon_{11}
  dEpsD2PlusDEps(0,1) = dL12dEps(1) * eValsD2(1, 1) + dL13dEps(1) * eValsD2(2, 2) +
                         L12 * dEVal2dEps(1) + L13 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(0,3) = dL12dEps(3) * eValsD2(1, 1) + dL13dEps(3) * eValsD2(2, 2) +
                         L12 * dEVal2dEps(3) + L13 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(0,5) = dL12dEps(5) * eValsD2(1, 1) + dL13dEps(5) * eValsD2(2, 2) +
                         L12 * dEVal2dEps(5) + L13 * dEVal3dEps(5);
  // Row 2
  // \epsilon_{11}
  dEpsD2PlusDEps(1,1) = dL22dEps(1) * eValsD2(1, 1) + dL23dEps(1) * eValsD2(2, 2) +
                         L22 * dEVal2dEps(1) + L23 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(1,3) = dL22dEps(3) * eValsD2(1, 1) + dL23dEps(3) * eValsD2(2, 2) +
                         L22 * dEVal2dEps(3) + L23 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(1,5) = dL22dEps(5) * eValsD2(1, 1) + dL23dEps(5) * eValsD2(2, 2) +
                         L22 * dEVal2dEps(5) + L23 * dEVal3dEps(5);
  // Row 3
  // \epsilon_{11}
  dEpsD2PlusDEps(2,1) = dL32dEps(1) * eValsD2(1, 1) + dL33dEps(1) * eValsD2(2, 2) +
                         L32 * dEVal2dEps(1) + L33 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(2,3) = dL32dEps(3) * eValsD2(1, 1) + dL33dEps(3) * eValsD2(2, 2) +
                         L32 * dEVal2dEps(3) + L33 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(2,5) = dL32dEps(5) * eValsD2(1, 1) + dL33dEps(5) * eValsD2(2, 2) +
                         L32 * dEVal2dEps(5) + L33 * dEVal3dEps(5);
  // Row 4
  // \epsilon_{11}
  dEpsD2PlusDEps(3,1) = dL42dEps(1) * eValsD2(1, 1) + dL43dEps(1) * eValsD2(2, 2) +
                         L42 * dEVal2dEps(1) + L43 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(3,3) = dL42dEps(3) * eValsD2(1, 1) + dL43dEps(3) * eValsD2(2, 2) +
                         L42 * dEVal2dEps(3) + L43 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(3,5) = dL42dEps(5) * eValsD2(1, 1) + dL43dEps(5) * eValsD2(2, 2) +
                         L42 * dEVal2dEps(5) + L43 * dEVal3dEps(5);
  // Row 5
  // \epsilon_{11}
  dEpsD2PlusDEps(4,1) = dL52dEps(1) * eValsD2(1, 1) + dL53dEps(1) * eValsD2(2, 2) +
                         L52 * dEVal2dEps(1) + L53 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(4,3) = dL52dEps(3) * eValsD2(1, 1) + dL53dEps(3) * eValsD2(2, 2) +
                         L52 * dEVal2dEps(3) + L53 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(4,5) = dL52dEps(5) * eValsD2(1, 1) + dL53dEps(5) * eValsD2(2, 2) +
                         L52 * dEVal2dEps(5) + L53 * dEVal3dEps(5);
  // Row 6
  // \epsilon_{11}
  dEpsD2PlusDEps(5,1) = dL62dEps(1) * eValsD2(1, 1) + dL63dEps(1) * eValsD2(2, 2) +
                         L62 * dEVal2dEps(1) + L63 * dEVal3dEps(1);
  // \gamma_{12}
  dEpsD2PlusDEps(5,3) = dL62dEps(3) * eValsD2(1, 1) + dL63dEps(3) * eValsD2(2, 2) +
                         L62 * dEVal2dEps(3) + L63 * dEVal3dEps(3);
  // \gamma_{23}
  dEpsD2PlusDEps(5,5) = dL62dEps(5) * eValsD2(1, 1) + dL63dEps(5) * eValsD2(2,2) +
                         L62 * dEVal2dEps(5) + L63 * dEVal3dEps(5);


}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::computeSEff(const Vector6d& stressEst,
                                    const VectorXd& dVals,
                                    Matrix6d& Seff) {
  // TODO: stress activation!!!
  Seff = S0 +
    dVals(0)*H1 + dVals(1)*H2 + dVals(2)*H4 + dVals(3)*H5;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel4Param::computeMatTang(const Matrix6d& Ceff,
                                       const Vector6d& epsStar,
                                       const Vector6d& epsStarMac,
                                       const VectorXd& yMaxVals,
                                       const VectorXd& gVals,
                                       Matrix6d& matTang) const {
  // compute ddsdde
  // material tangent computed analytically
  // lots of math here, hard parts done in sympy
  //
  // ddsdde = dCeff/de*e + Ceff
  // dCeff/de = d(inv(Seff))/d(d)*d(d)/dg*dg/dy*dy/dz*dz/de
  // (sum for each d_i)
  //

  // d(inv(Seff))/d(d) = -inv(Seff)*dSeff/d(d)*inv(Seff)
  Matrix6d dInvSeffdd1 = -1.0*Ceff*H1*Ceff;
  Matrix6d dInvSeffdd2 = -1.0*Ceff*H2*Ceff;
  Matrix6d dInvSeffdd4 = -1.0*Ceff*H4*Ceff;
  Matrix6d dInvSeffdd5 = -1.0*Ceff*H5*Ceff;

  // d(d)/dg
  // TODO: Optimize this?
  double dd1dg1 = pe(0)*dc(0)*pow(gVals(0), pe(0) - 1.0)
    *exp(-1.0*pow(gVals(0),pe(0)));
  double dd2dg2 = pe(1)*dc(1)*pow(gVals(1), pe(1) - 1.0)
    *exp(-1.0*pow(gVals(1),pe(1)));
  double dd4dg4 = pe(2)*dc(2)*pow(gVals(2), pe(2) - 1.0)
    *exp(-1.0*pow(gVals(2),pe(2)));
  double dd5dg5 = pe(3)*dc(3)*pow(gVals(3), pe(3) - 1.0)
    *exp(-1.0*pow(gVals(3),pe(3)));

  // dg/dy
  // if yMax < yMaxOld, yMax = yMaxOld, dg/dy = 0.0
  double dg1dy1 = 0.0, dg2dy2 = 0.0, dg4dy4 = 0.0, dg5dy5 = 0.0;
  if (yMaxVals(0) > yMaxSave(0)) {
    dg1dy1 = 0.5*1.0/sqrt(yMaxVals(0)*yc(0));
  }
  if (yMaxVals(1) > yMaxSave(1)) {
    dg2dy2 = 0.5*1.0/sqrt(yMaxVals(1)*yc(1));
  }
  if (yMaxVals(2) > yMaxSave(2)) {
    dg4dy4 = 0.5*1.0/sqrt(yMaxVals(2)*yc(2));
  }
  if (yMaxVals(3) > yMaxSave(3)) {
    dg5dy5 = 0.5*1.0/sqrt(yMaxVals(3)*yc(3));
  }

  // dy/dz
  double dy1dz1 = 1.0, dy1dz6 = -1.0, dy2dz2 = 1.0, dy2dz6 = -1.0;
  double dy4dz6 = yMaxVals(2) > 0.0 ? 1.0 : 0.0,
    dy5dz6 = yMaxVals(3) > 0.0 ? 1.0 : 0.0;
  
  // dz/dEpsD1/D2
  Vector6d epsD1 = Vector6d::Zero(), epsD2 = Vector6d::Zero();
  // [eps11, 0, 0, gam12, gam13, 0]
  epsD1(0) = epsStar(0); epsD1(3) = epsStar(3); epsD1(4) = epsStar(4);
  // [0, eps22, 0, gam12, 0, gam23]
  epsD2(1) = epsStar(1); epsD2(3) = epsStar(3); epsD2(5) = epsStar(5);
  Vector6d dz1dEpsD1 = Vector6d::Zero(), dz2dEpsD2 = Vector6d::Zero(),
    dz6dEpsD1 = Vector6d::Zero(), dz6dEpsD2 = Vector6d::Zero();
  // z1
  dz1dEpsD1(0) = C0(0,0)*epsD1(0);
  dz1dEpsD1(3) = b(0)*C0(3,3)*epsD1(3);
  dz1dEpsD1(4) = b(1)*C0(4,4)*epsD1(4);
  // z2
  dz2dEpsD2(0) = C0(1,1)*epsD2(1);
  dz2dEpsD2(3) = b(0)*C0(3,3)*epsD2(3);
  dz2dEpsD2(5) = b(1)*C0(5,5)*epsD2(5);
  //z6 - epsD1
  dz6dEpsD1(0) = 0.25*(C0(0,0) + b(2)*C0(3,3))*epsD1(3);
  dz6dEpsD1(3) = 0.25*(C0(0,0) + b(2)*C0(3,3))*epsD1(0);
  //z6 - epsD2
  dz6dEpsD2(1) = 0.25*(C0(1,1) + b(2)*C0(3,3))*epsD2(3);
  dz6dEpsD2(3) = 0.25*(C0(1,1) + b(2)*C0(3,3))*epsD2(1);

  //dEpsD*_i / dEps_j
  Matrix6d dEpsD1dEps = Matrix6d::Zero(), dEpsD2dEps = Matrix6d::Zero(); 
  calcDEpsD1PlusDEps(epsD1, dEpsD1dEps);
  calcDEpsD2PlusDEps(epsD2, dEpsD2dEps);
  
  // MatrixVectorProd
  Eigen::Tensor<double,3> dInvSeff1de(6,6,6), dInvSeff2de(6,6,6), dInvSeff4de(6,6,6), dInvSeff5de(6,6,6);
  dInvSeff1de.setZero(); dInvSeff4de.setZero(); dInvSeff4de.setZero(); dInvSeff5de.setZero();
  Vector6d dz1dEps = dEpsD1dEps.transpose()*dz1dEpsD1;
  Vector6d dz2dEps = dEpsD2dEps.transpose()*dz2dEpsD2;
  Vector6d dz6dEps = dEpsD1dEps.transpose()*dz6dEpsD1 + dEpsD2dEps.transpose()*dz6dEpsD2;
  // NOTE: dz1/dEpsD1 -> 6x1, dEpsD1/dEps -> 6x6
  // dEpsD1_i/dEps_j * dz1/dEpsD1_i is what I want
  for (int iInd = 0; iInd < 6; iInd++) {
    for (int jInd = 0; jInd < 6; jInd++) {
      for (int kInd = 0; kInd < 6; kInd++) {
        dInvSeff1de(iInd,jInd,kInd) = dInvSeffdd1(iInd,jInd)*dd1dg1*
          dg1dy1*(dy1dz1*dz1dEps(kInd) + dy1dz6*dz6dEps(kInd));
        dInvSeff2de(iInd,jInd,kInd) = dInvSeffdd2(iInd, jInd)*dd2dg2*
          dg2dy2*(dy2dz2*dz2dEps(kInd) + dy2dz6*dz6dEps(kInd));
        dInvSeff4de(iInd,jInd,kInd) = dInvSeffdd4(iInd,jInd)*dd4dg4*
          dg4dy4*dy4dz6*dz6dEps(kInd);
        dInvSeff5de(iInd,jInd,kInd) = dInvSeffdd5(iInd,jInd)*dd5dg5*
          dg5dy5*dy5dz6*dz6dEps(kInd);
      }
    }
  }
  
  // dCeff/de
  Eigen::Tensor<double,3> dCeffde = dInvSeff1de + dInvSeff2de + dInvSeff4de + dInvSeff5de;
  // dCeff/de . e
  Matrix6d dCeffdEpsDotEps = Matrix6d::Zero();
  for (int iInd = 0; iInd < 6; iInd++) {
    for (int jInd = 0; jInd < 6; jInd++) {
      for (int kInd = 0; kInd < 6; kInd++) {
        dCeffdEpsDotEps(iInd, jInd) += dCeffde(iInd, jInd, kInd)*epsStar(kInd);
      }
    }
  }
  // Material Tangent = dCeff/de . e + Ceff
  matTang = dCeffdEpsDotEps + Ceff;
}
/********************************************************************/
/********************************************************************/


