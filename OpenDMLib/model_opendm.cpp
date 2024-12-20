// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// Implementation file for OpenDM Model class
//

#include <iostream>
#include <math.h>

#include "model_opendm.hpp"

// Some things of note:
// Abaqus stress conv: S11, S22, S33, S12, S13, S23
//
using std::cout;
using std::endl;
using std::complex;

#define CHECK_TANG_EVALS 0

#if CHECK_TANG_EVALS
#include <Eigen/Eigenvalues>
inline void checkTangEvals(Matrix6d& matTang, Vector6d& epsStar) {
  Eigen::EigenSolver<Matrix6d> eigensolver(matTang);
  std::complex<double> eVal;
  for (int iEig = 0; iEig < 6; ++iEig) {
    eVal = eigensolver.eigenvalues().col(0)[iEig];
    if (real(eVal) <= 0.1) {
      cout << "epsStar = \n" << epsStar << endl; 
      cout << "matTang = \n" << matTang << endl; 
      cout << "eVals = \n" << eigensolver.eigenvalues().col(0) << endl; 
      // eigenvectors
      Eigen::Matrix<complex<double>, 6, 6> eVect = eigensolver.eigenvectors();
      // transform back from spectral
      Eigen::Matrix<complex<double>, 6, 1>  epsSpect = eVect*epsStar;
      cout << "epsSpect = \n" << epsSpect << endl; 
      break;
      throw std::runtime_error("Negative eVal!!!");
    }
  }
}
#endif


/********************************************************************/
/********************************************************************/

OpenDMModel::OpenDMModel(double* props, int* nprops, double* statev,
                         int* nstatv, const int nDVarsIn)
  : nDamageVars(nDVarsIn) {

  if (nDVarsIn == 2) {
    if ((*nprops) != 25) {
      std::cout << "OpenDMModel::OpenDMModel: NOT ENOUGH PROPS!!!"
        " using 2 param model " << std::endl;
    }
    if ((*nstatv) != 13) {
      std::cout << "NOT ENOUGH STATEVARS!!!"
        " using 2 param model " << std::endl;
    }
  } else if (nDVarsIn == 4) {
    if ((*nprops) != 42) {
      std::cout << "OpenDMModel::OpenDMModel: NOT ENOUGH PROPS!!!"
        " using 4 param model " << std::endl;
    }
    if ((*nstatv) != 29) {
      std::cout << "NOT ENOUGH STATEVARS!!!"
        " using 4 param model " << std::endl;
    }
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
  S0(1,0) = -nu12*invE11; S0(0,1) = S0(1,0);
  S0(2,0) = -nu13*invE11; S0(0,2) = S0(2,0);
  S0(2,1) = -nu23*invE22; S0(1,2) = S0(2,1);
  S0(3,3) = invG12; S0(4,4) = invG13; S0(5,5) = invG23;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::unpackStateVars(double* statev) {
  // statev[0:8] = CeffOld11, C22, C33, C23, C13, C12, C44, C55, C66
  CeffOld = Matrix6d::Zero();
  CeffOld(0,0) = statev[0];
  CeffOld(1,1) = statev[1];
  CeffOld(2,2) = statev[2];
  CeffOld(1,2) = statev[3]; CeffOld(2,1) = CeffOld(1,2);
  CeffOld(0,2) = statev[4]; CeffOld(2,0) = CeffOld(0,2);
  CeffOld(0,1) = statev[5]; CeffOld(1,0) = CeffOld(0,1);
  CeffOld(3,3) = statev[6];
  CeffOld(4,4) = statev[7];
  CeffOld(5,5) = statev[8];
  // statev[9:9+nDamageVars] = yMaxOld
  yMaxSave.resize(nDamageVars);
  for (int iY = 0; iY < nDamageVars; iY++) {
    yMaxSave(iY) = statev[9+iY];
  }
  // I do not need to unpack old damage vars...

}
/********************************************************************/
/********************************************************************/

void OpenDMModel::updateStateVars(const VectorXd& yMax, const Matrix6d& Ceff,
                  const VectorXd& dVals, double* statev) {
  // statev[0:8] = Ceff11, C22, C33, C23, C13, C12, C44, C55, C66
  statev[0] = Ceff(0,0);
  statev[1] = Ceff(1,1);
  statev[2] = Ceff(2,2);
  statev[3] = Ceff(1,2);
  statev[4] = Ceff(0,2);
  statev[5] = Ceff(0,1);
  statev[6] = Ceff(3,3);
  statev[7] = Ceff(4,4);
  statev[8] = Ceff(5,5);

  // statev[9:9+nDamageVars] = y1MaxOld, y2MaxOld,...
  // Start at statev[9] and assign from here
  for (int iPos = 0; iPos < 2*nDamageVars; iPos++) {
    if (iPos < nDamageVars) {
      // check yMax
      if (yMax(iPos) > yMaxSave(iPos)) {
        yMaxSave(iPos) = yMax(iPos);
      }
      // update yMax Values
      statev[9+iPos] = yMaxSave(iPos);
    } else {
      // update damage vars
      statev[9+iPos] = dVals(iPos - nDamageVars);
    }
  }
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

  // epsStar = eps - epsTh (Abq send these by default)
  Vector6d epsStar; 
  // Strains specific to damage modes
  Vector6d epsD1Plus = Vector6d::Zero(), epsD2Plus = Vector6d::Zero();
  // Update strain
  for (int i = 0; i < 6; i++) {
    // NO THERMAL STRAIN AT THIS POINT!!!
    epsStar(i) = strain[i] + dstrain[i];
  }
  // Estimate stress for activation functions
  Vector6d stressEst = CeffOld*epsStar;
  
  // update yMax attribute values here
  VectorXd yMaxVals(nDamageVars),
    gVals(nDamageVars), dVals(nDamageVars);
  
  yMaxVals = calcDrivingForces(epsStar, epsD1Plus, epsD2Plus);
  // get g and d values
  calcGVals(yMaxVals, gVals);
  calcDVals(gVals, dVals);

  // Get Seff, Ceff
  Matrix6d Seff = Matrix6d::Zero();
  computeSEff(stressEst, dVals, Seff);
  Matrix6d Ceff = Matrix6d::Zero();
  matrixInverse(Seff, Ceff);

  // calcStress
  Vector6d sig = Ceff*epsStar;

  // calcTangent
  Matrix6d matTang = Matrix6d::Zero();
#if defined USE_NUM_TANG
  computeNumMatTang(Ceff, epsStar, sig, matTang);
#else
  computeMatTang(Ceff, epsStar, epsD1Plus, epsD2Plus, yMaxVals,
                 gVals, matTang);
#endif
  // cout << "Ceff = \n" << Ceff << endl;
  // cout << "matTang = \n" << matTang << endl;

  // NOTE: as formulated, tangent is not sym, making sym here
  Matrix6d symTang = 0.5*(matTang + matTang.transpose());
#if CHECK_TANG_EVALS
  checkTangEvals(symTang, epsStar);
#endif

  
  // Update stateVars (yMax, Ceff)
  updateStateVars(yMaxVals, Ceff, dVals, statev);

  // update stress, tang, energies
  // Neglecting scd, no time effects
  // TODO: Also skipping spd for now
  setUMATOuts(sig, epsStar, symTang, stress, ddsdde, sse);
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::calcGVals(const VectorXd& yMax, VectorXd& gVals) const {
  // EQ 34 from OpenDM-Simple CLT doc
  for (int iVal = 0; iVal < nDamageVars; iVal++) {
    double yDiff = macaulayBracketPlus(sqrt(yMax(iVal)) - sqrt(y0(iVal)));
    gVals(iVal) = yDiff/sqrt(yc(iVal));
  }
  // std::cout << "g = " << gVals << std::endl;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::calcDVals(const VectorXd& gVals, VectorXd& dVals) const {
  // EQ 35 from OpenDM-Simple CLT doc
  for (int iVal = 0; iVal < nDamageVars; iVal++) {
    double gExp = -1.0*pow(gVals(iVal), pe(iVal));
    dVals(iVal) = dc(iVal)*(1.0 - std::exp(gExp));
  }
  // std::cout << "d = " << dVals << std::endl;
}
/********************************************************************/
/********************************************************************/

void OpenDMModel::computeNumMatTang(const Matrix6d& Ceff,
                                    const Vector6d& epsStar,
                                    const Vector6d& stress,
                                    Matrix6d& numMatTang) {
  // Compute numerical material tangent via repeated stressCalcs
  // at assigned strain perturbations
  //
  const double strainPert = 5.0e-7;
  // initialize strain
  Vector6d epsPert = Vector6d::Zero();
  // loop over perturbations
  for (unsigned int iPert = 0; iPert < 6; ++iPert) {
    // set to step strain
    epsPert = epsStar;
    // perturb strain entry
    epsPert[iPert] += strainPert;
    // cout << "perturbed strain = " << epsPert << endl;
    // update yMax attribute values here
    VectorXd yMaxVals(nDamageVars),
      gVals(nDamageVars), dVals(nDamageVars);
  
    // Strains specific to damage modes
    Vector6d epsD1Plus = Vector6d::Zero(), epsD2Plus = Vector6d::Zero();
    yMaxVals = calcDrivingForces(epsPert, epsD1Plus, epsD2Plus);
    // get g and d values
    calcGVals(yMaxVals, gVals);
    calcDVals(gVals, dVals);
    // cout << "yMaxSave = " << yMaxSave(0) << yMaxSave(1) << yMaxSave(2) << yMaxSave(3) << endl;
    // cout << "yMax = " << yMaxVals(0) << yMaxVals(1) << yMaxVals(2) << yMaxVals(3) << endl;
    // cout << "d = " << dVals(0) << dVals(1) << dVals(2) << dVals(3) << endl;

    // Get Seff, Ceff
    Matrix6d Seff = Matrix6d::Zero();
    computeSEff(stress, dVals, Seff);
    Matrix6d Ceff = Matrix6d::Zero();
    matrixInverse(Seff, Ceff);
    // cout << "Ceff = \n" << Ceff << endl;

    // calcStress
    Vector6d sigPert = Ceff*epsPert;
    // cout << "perturbed stress = " << sigPert << endl;
    // numMatTang components
    for (unsigned int iStress = 0; iStress < 6; ++iStress) {
      numMatTang(iPert, iStress) = (sigPert(iStress) - stress(iStress))/strainPert;
    }
  }

}
/********************************************************************/
/********************************************************************/

