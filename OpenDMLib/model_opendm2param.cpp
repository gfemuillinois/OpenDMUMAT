// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// Implementation file for OpenDM Model with 2 damage parameters
// This is derived from OpenDMModel class
// It captures damage in the 0 & 90 directions in plane of composite
//

#include "model_opendm2param.hpp"
#include <iostream>

/********************************************************************/
/********************************************************************/

OpenDMModel2Param::OpenDMModel2Param(double* props, int* nprops,
				     double* statev, int* nstatv)
    : OpenDMModel(props, nprops, statev, nstatv, 2 /*nDamageVals*/)
{
  // Set OpenDM model params
  unpackParams(props);

  // Build base H1 & H2
  createH1H2();

}
/********************************************************************/
/********************************************************************/

void OpenDMModel2Param::unpackParams(double* props) {
  // props[9:20] = h_s1, h_s2, b_1, b_2, y01, y02,
  // yc1, yc2, pe1, pe2, dc1, dc2
  // unpack props
  hs1 = Eigen::Map<Eigen::VectorXd>(props+9,1);
  hs2 = Eigen::Map<Eigen::VectorXd>(props+10,1);
  b = Eigen::Map<Eigen::VectorXd>(props+11,2);
  y0 = Eigen::Map<Eigen::VectorXd>(props+13,2);
  yc = Eigen::Map<Eigen::VectorXd>(props+15,2);
  pe = Eigen::Map<Eigen::VectorXd>(props+17,2);
  dc = Eigen::Map<Eigen::VectorXd>(props+19,2);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel2Param::createH1H2() {
  // H1 & H2 without stress activation
  H1 = Matrix6d::Zero();
  H2 = Matrix6d::Zero();
  // H1
  // Mode I
  H1(0,0) = S0(0,0);
  // Mode III
  H1(3,3) = hs1(0)*S0(3,3);
  // Mode II
  H1(4,4) = hs1(0)*S0(4,4);
  // H2
  // Mode I
  H2(1,1) = S0(1,1);
  // Mode II
  H2(3,3) = hs2(0)*S0(3,3);
  // Mode III
  H2(5,5) = hs2(0)*S0(5,5);

}
/********************************************************************/
/********************************************************************/

VectorXd OpenDMModel2Param::calcDrivingForces(const Vector6d& epsStarMac) {
  // EQ 32 in OpenDM-Simple CLT doc
  double E11 = 1.0/S0(0,0);
  double E22 = 1.0/S0(1,1);
  double G12 = C0(3,3);

  double e11 = epsStarMac(0);
  double e22 = epsStarMac(1);
  double g12 = epsStarMac(3);
  
  double y1 = 0.5*(E11*e11*e11 + b(0)*G12*g12*g12);
  double y2 = 0.5*(E22*e22*e22 + b(1)*G12*g12*g12);

  VectorXd yMax(2);
  yMax(0) = y1 > yMaxSave(0) ? y1 : yMaxSave(0);
  yMax(1) = y2 > yMaxSave(1) ? y2 : yMaxSave(1);
  return yMax;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel2Param::computeSEff(const Vector6d& stressEst,
				    const VectorXd& dVals,
				    Matrix6d& Seff) {
  // stress activation
  // H1
  // Mode I
  H1(0,0) = stressAct(stressEst(0))*S0(0,0);
  // H2
  // Mode I
  H2(1,1) = stressAct(stressEst(1))*S0(1,1);

  // Seff = S0 + \sum_i d_i H_i
  Seff = S0 + dVals(0)*H1 + dVals(1)*H2;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel2Param::computeMatTang(const Matrix6d& Ceff, const Vector6d& epsStar,
				       const Vector6d& epsStarMac, const VectorXd& yMaxVals,
				       const VectorXd& gVals,
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


