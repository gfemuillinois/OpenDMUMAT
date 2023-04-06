// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// Implementation file for OpenDM Model class
//

#include <iostream>
#include <math.h>

#include "model_opendm.hpp"
#include "math_opendm.hpp"


// Some things of note:
// Abaqus stress conv: S11, S22, S33, S12, S13, S23
//

/********************************************************************/
/********************************************************************/

OpenDMModel::OpenDMModel(double* props, int* nprops, double* statev,
                         int* nstatv) {

  if ((*nprops) != 21) {
    std::cout << "NOT ENOUGH PROPS!!!" << std::endl;
  }
  if ((*nstatv) != 13) {
    std::cout << "NOT ENOUGH STATEVARS!!!" << std::endl;
  }
  // props[0:9] = E11, E22, E33, nu_12, nu_23, nu_13, G_12, G_23, G_13
  // props[9:20] = h_s1, h_s2, b_1, b_2, y01, y02, yc1, yc2, pe1, pe2, dc1, dc2
  // statev = y_maxi, CeffOld, strainOld
  //
  // Calc elastic compliance tensor
  // e11, e22, e33, g12, g13, g23
  computeS0(props);
  Matrix6d Ctemp = Matrix6d::Zero();
  matrixInverse(S0, Ctemp);
  C0 = Ctemp;

  // Set OpenDM model params
  unpackParams(props);

  // set stateVars
  unpackStateVars(statev);
 
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::computeS0(double* props) {
  // Build orthotropic compliance tensor
  S0 = Matrix6d::Zero();
  double invE11 = 1.0/props[0];
  double invE22 = 1.0/props[1];
  double invE33 = 1.0/props[2];
  double nu12, nu23, nu13;
  nu12 = props[3]; nu23 = props[4]; nu13 = props[5];
  double invG12 = 1.0/props[6];
  double invG23 = 1.0/props[7];
  double invG13 = 1.0/props[8];
  // Set matrix
  S0(0,0) = invE11; S0(1,1) = invE22; S0(2,2) = invE33;
  S0(1,0) = nu12*invE11; S0(0,1) = S0(1,0);
  S0(2,0) = nu13*invE11; S0(0,2) = S0(2,0);
  S0(2,1) = nu23*invE22; S0(1,2) = S0(2,1);
  S0(3,3) = invG12; S0(4,4) = invG13; S0(5,5) = invG23;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::unpackParams(double* props) {
  // props[9:20] = h_s1, h_s2, b_1, b_2, y01, y02,
  // yc1, yc2, pe1, pe2, dc1, dc2
  // unpack props
  hs = {props[9], props[10]};
  b = {props[11], props[12]};
  y0 = {props[13], props[14]};
  yc = {props[15], props[16]};
  pe = {props[17], props[18]};
  dc = {props[19], props[20]};

}
/********************************************************************/
/********************************************************************/

void OpenDMModel::unpackStateVars(double* statev) {
  // statev[0:1] = y1MaxOld, y2MaxOld
  yMaxSave = {statev[0], statev[1]};
  // statev[2:10] = CeffOld11, C22, C33, C23, C13, C12, C44, C55, C66
  CeffOld = Matrix6d::Zero();
  CeffOld(0,0) = statev[2];
  CeffOld(1,1) = statev[3];
  CeffOld(2,2) = statev[4];
  CeffOld(1,2) = statev[5]; CeffOld(2,1) = CeffOld(1,2);
  CeffOld(0,2) = statev[6]; CeffOld(2,0) = CeffOld(0,2);
  CeffOld(0,1) = statev[7]; CeffOld(1,0) = CeffOld(0,1);
  CeffOld(3,3) = statev[8];
  CeffOld(4,4) = statev[9];
  CeffOld(5,5) = statev[10];
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::updateStateVars(const Vector2d& yMax, const Matrix6d& Ceff,
				                  const Vector2d& dVals, double* statev) {
  // statev[0:1] = y1MaxOld, y2MaxOld
  // check max
  if (yMax(0) > yMaxSave(0)) {
    yMaxSave(0) = yMax(0);
  }
  if (yMax(1) > yMaxSave(1)) {
    yMaxSave(1) = yMax(1);
  }
  statev[0] = yMaxSave[0]; statev[1] = yMaxSave[1];

  // statev[2:10] = Ceff11, C22, C33, C23, C13, C12, C44, C55, C66
  statev[2] = Ceff(0,0);
  statev[3] = Ceff(1,1);
  statev[4] = Ceff(2,2);
  statev[5] = Ceff(1,2);
  statev[6] = Ceff(0,2);
  statev[7] = Ceff(0,1);
  statev[8] = Ceff(3,3);
  statev[9] = Ceff(4,4);
  statev[10] = Ceff(5,5);

  // statev[11:12] = d1, d2
  statev[11] = dVals(0);
  statev[12] = dVals(1);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::setUMATOuts(const Vector6d& sig, const Vector6d& epsStar,
			      const Matrix6d& matTang, double* stress,
			      double* ddsdde, double* sse) const {
  // assign vector positions to pointer outs
  //
  for (int iRow = 0; iRow < 6; iRow++) {
    stress[iRow] = sig(iRow);
    ddsdde[iRow] = matTang(iRow, 0);
    ddsdde[6+iRow] = matTang(iRow, 1);
    ddsdde[12+iRow] = matTang(iRow, 2);
    ddsdde[18+iRow] = matTang(iRow, 3);
    ddsdde[24+iRow] = matTang(iRow, 4);
    ddsdde[30+iRow] = matTang(iRow, 5);
  }

  // strain energy density
  (*sse) = 0.5*sig.dot(epsStar);
  
  }
/********************************************************************/
/********************************************************************/

void OpenDMModel::runModel(double* strain, double* dstrain, double* stress,
			   double* statev, double* ddsdde, double* sse,
			   double* spd, double* scd) {
  // Calculate y, g, d, Ceff, upStress, matTang
  //

  // epsStar = eps - epsTh (formality for now)
  Vector6d epsStar, epsStarMac;
  // Update strain
  for (int i = 0; i < 6; i++) {
    epsStar(i) = strain[i] + dstrain[i];
  }
  // Estimate stress for activation functions
  Vector6d stressEst = CeffOld*epsStar;
  
  // Macauley strains
  epsStarMac = epsStar;
  // Only check normal strains
  for (int iRow = 0; iRow < 3; iRow++) {
    epsStarMac(iRow) = macaulayBracket(epsStar(iRow));
  }
  // update yMax attribute values here
  Vector2d yMaxVals = calcYVals(epsStarMac);
  // get g and d values
  Vector2d gVals = calcGVals(yMaxVals);
  Vector2d dVals = calcDVals(gVals);

  // Make two H matrices
  Matrix6d H1, H2;
  calcH1H2(stressEst, H1, H2);

  // Get Seff, Ceff
  Matrix6d Seff = S0 + dVals(0)*H1 + dVals(1)*H2;
  Matrix6d Ceff = Matrix6d::Zero();
  matrixInverse(Seff, Ceff);

  // calcStress
  Vector6d sig = Ceff*epsStar;
  // calcTangent
  Matrix6d matTang = Matrix6d::Zero();
  computeMatTang(Ceff, H1, H2, epsStar, epsStarMac, gVals, yMaxVals, matTang);

  // Update stateVars (yMax, Ceff)
  updateStateVars(yMaxVals, Ceff, dVals, statev);

  // update stress, tang, energies
  // Neglecting scd, no time effects
  // TODO: Also skipping spd for now
  setUMATOuts(sig, epsStar, matTang, stress, ddsdde, sse);
}
/********************************************************************/
/********************************************************************/

Vector2d OpenDMModel::calcYVals(const Vector6d& epsStarMac) const {
  // EQ 32 in OpenDM-Simple CLT doc
  double E11 = 1.0/S0(0,0);
  double E22 = 1.0/S0(1,1);
  double G12 = C0(3,3);

  double e11 = epsStarMac(0);
  double e22 = epsStarMac(1);
  double g12 = epsStarMac(3);
  
  double y1 = 0.5*(E11*e11*e11 + b(0)*G12*g12*g12);
  double y2 = 0.5*(E22*e22*e22 + b(1)*G12*g12*g12);

  Vector2d yMax;
  yMax(0) = y1 > yMaxSave(0) ? y1 : yMaxSave(0);
  yMax(1) = y2 > yMaxSave(1) ? y2 : yMaxSave(1);
  return yMax;
}
/********************************************************************/
/********************************************************************/

Vector2d OpenDMModel::calcGVals(const Vector2d& yMax) const {
  // EQ 34 from OpenDM-Simple CLT doc
  Vector2d gVals = Vector2d::Zero();
  for (int iVal = 0; iVal < 2; iVal++) {
    double yDiff = macaulayBracket(sqrt(yMax(iVal)) - sqrt(y0(iVal)));
    gVals(iVal) = yDiff/sqrt(yc(iVal));
  }

  return gVals;
}
/********************************************************************/
/********************************************************************/

Vector2d OpenDMModel::calcDVals(const Vector2d& gVals) const {
  // EQ 35 from OpenDM-Simple CLT doc
  Vector2d dVals = Vector2d::Zero();
  for (int iVal = 0; iVal < 2; iVal++) {
    double gExp = -1.0*pow(gVals(iVal), pe(iVal));
    dVals(iVal) = dc(iVal)*(1.0 - exp(gExp));
  }
  return dVals;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::calcH1H2(const Vector6d& stressEst,
			   Matrix6d& H1,
			   Matrix6d& H2) const {
  // zero out inputs
  H1 = Matrix6d::Zero();
  H2 = Matrix6d::Zero();
  // H1
  // Mode I
  H1(0,0) = stressAct(stressEst(0))*S0(0,0);
  // Mode III
  H1(3,3) = hs(0)*S0(3,3);
  // Mode II
  H1(4,4) = hs(0)*S0(4,4);
  // H2
  // Mode I
  H2(1,1) = stressAct(stressEst(1))*S0(1,1);
  // Mode II
  H2(3,3) = hs(1)*S0(3,3);
  // Mode III
  H2(5,5) = hs(1)*S0(5,5);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::computeMatTang(const Matrix6d& Ceff, const Matrix6d& H1,
				 const Matrix6d& H2, const Vector6d& epsStar,
				 const Vector6d& epsStarMac, const Vector2d yMaxVals,
				 const Vector2d gVals,
				 Matrix6d& matTang) const {
  // compute ddsdde
  // material tangent computed analytically
  // Some fancy math in here but it is doable and basically just chain
  // rule
  //
  // ddsdde = dCeff/de*e + Ceff
  // dCeff/de = d(inv(Seff))/d(d)*d(d)/dg*dg/dy*dy/de (sum for each d_i)
  //

  // d(inv(Seff))/d(d) = -inv(Seff)*dSeff/d(d)*inv(Seff)
  Matrix6d dInvSeffdd1 = -1.0*Ceff*H1*Ceff;
  Matrix6d dInvSeffdd2 = -1.0*Ceff*H2*Ceff;

  // d(d)/dg
  // TODO: Optimize this?
  double dd1dg1 = pe(0)*dc(0)*pow(gVals(0), pe(0) - 1.0)*exp(-1.0*pow(gVals(0),pe(0)));
  double dd2dg2 = pe(1)*dc(1)*pow(gVals(1), pe(1) - 1.0)*exp(-1.0*pow(gVals(1),pe(1)));

  // dg/dy
  // if yMax < yMaxOld, yMax = yMaxOld, dg/dy = 0.0
  double dg1dy1 = 0.0, dg2dy2 = 0.0;
  if (yMaxVals(0) > yMaxSave(0)) {
    dg1dy1 = 0.5*1.0/sqrt(yMaxVals(0)*yc(0));
  }
  if (yMaxVals(1) > yMaxSave(1)) {
    dg2dy2 = 0.5*1.0/sqrt(yMaxVals(1)*yc(1));
  }

  // dy/de
  // epsStarMac norm strains are zero if compressive
  Vector6d dy1de = Vector6d::Zero(), dy2de = Vector6d::Zero();
  dy1de(0) = 1.0/S0(0,0)*epsStarMac(0); dy1de(3) = b(0)*1.0/S0(3,3)*epsStarMac(3);
  dy2de(1) = 1.0/S0(1,1)*epsStarMac(1); dy2de(3) = b(1)*1.0/S0(3,3)*epsStarMac(3);

  // Combine scalar derivs to vect
  Vector6d dd1de = Vector6d::Zero(), dd2de = Vector6d::Zero();
  dd1de = dd1dg1*dg1dy1*dy1de;
  dd2de = dd2dg2*dg2dy2*dy2de;

  // MatrixVectorProd
  Vector6d dInvSeff1de = Vector6d::Zero(), dInvSeff2de = Vector6d::Zero();
  dInvSeff1de = dInvSeffdd1*dd1de;
  dInvSeff2de = dInvSeffdd2*dd2de;

  // dCeff/de
  Vector6d dCeffde = dInvSeff1de + dInvSeff2de;

  // Material Tangent = dCeff/de . e + Ceff
  matTang = dCeffde*epsStar.transpose() + Ceff;

}
/********************************************************************/
/********************************************************************/

/** @brief Init OpenDM Object and run model from cpp source
 *  To facilitate simple umat file that links with sharedLib with all of this stuff
 */
extern "C" {
  void runInitOpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
			  double* strain, double* dstrain, double* stress,
			  double* ddsdde, double* sse, double* spd, double* scd) {

    // Create OpenDM model object
    OpenDMModel* p_openDMModel = new OpenDMModel(props, nprops, statev, nstatv);
    // run openDMModel
    p_openDMModel->runModel(strain, dstrain, stress, statev, ddsdde, sse, spd, scd);
    // delete OpenDMModel object
    delete p_openDMModel;
  }
}
