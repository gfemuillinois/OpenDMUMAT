// - Bryce Mazurowski <brycepm2@gmail.com> -
// 
// Class that implements the OpenDM model
// This assumes interface with a UMAT call from somewhere else
//
#ifndef MODEL_OPENDM_H
#define MODEL_OPENDM_H

#include <Eigen/Dense>
#include "math_opendm.hpp"

// 6x6 matrix of doubles
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
// 3x3 matrix of doubles
typedef Eigen::Matrix<double, 3, 3> Matrix3d;
// 6 row col vector of doubles
typedef Eigen::Matrix<double, 6, 1> Vector6d;
// 6 row col vector of doubles
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

class OpenDMModel {

public:
  // Constructor
  OpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
	      const int nDamageVars);

  // default destrutor
  virtual ~OpenDMModel()=default;

  // default copy constructor (should never be used)
  OpenDMModel(const OpenDMModel&)=delete;
  
  // default operator= (should never be used)
  OpenDMModel& operator=(const OpenDMModel&)=delete;

  /** @brief run the openDM model for a given UMAT call
   */
  void runModel(double* strain, double* dstrain, double* stress,
		double* statev, double* ddsdde, double* sse,
		double* spd, double* scd);

protected:
  
  /** @brief get activation function for stress
   * basically just to affect stiffness only in tension
   */
  double stressAct(const double& stressVal) const {
    // Heaviside from stressVal
    return stressVal > 0.0 ? 1.0 : 0.0;
  }

  // number of damage variables for a given model
  // 2Param or 4Param
  const int nDamageVars;


  // compliance tensor
  Matrix6d S0;
  // stiffness tensor
  Matrix6d C0;
  
  // Vectors of all params for each damage var
  VectorXd hs1, hs2, b, y0, yc, pe, dc;
  // vector of yMax (stateVar)
  VectorXd yMaxSave;
  // matrix of prev Ceff (stateVar)
  Matrix6d CeffOld;

private:
  /** @brief computeS0 from UMAT array props
   */
  void computeS0(double* props);

   /** @brief unpack stateVars from UMAT call
   */
  void unpackStateVars(double* statev);

  /** @brief update statevars ater run
   */
  void updateStateVars(const VectorXd& yMax, const Matrix6d& Ceff,
		       const VectorXd& dVals, double* statev);

  /** @brief set outputs from UMAT
   */
  void setUMATOuts(const Vector6d& sig, const Vector6d& epsStar,
		   const Matrix6d& matTang, double* stress,
		   double* ddsdde, double* sse) const;

  /** @brief calculate g_i for given strain
   */
  void calcGVals(const VectorXd& yMax, VectorXd& gVals) const;

  /** @brief calculate damage values for each damage param
   */
  void calcDVals(const VectorXd& gVals, VectorXd& dVals) const;

  /** @brief unpack model parameters from props
   */
  virtual void unpackParams(double* props)=0;

  /** @brief calculate yMax_i for given strain
   */
  virtual VectorXd calcDrivingForces(const Vector6d& epsStarMac)=0;

  /** @brief calculate damage effect tensor H1, H2
   */
  virtual void computeSEff(const Vector6d& stressEst, const VectorXd& dVals,
		Matrix6d& Seff) =0;

  /** @brief compute analytical material tangent stiffness
   */
  virtual void computeMatTang(const Matrix6d& Ceff, const Vector6d& epsStar,
			      const Vector6d& epsStarMac, const VectorXd& yMaxVals,
			      const VectorXd& gVals,
			      Matrix6d& matTang)=0;

};
#endif
