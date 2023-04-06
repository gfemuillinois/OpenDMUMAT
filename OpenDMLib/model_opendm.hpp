// - Bryce Mazurowski <brycepm2@gmail.com> -
// 
// Class that implements the OpenDM model
// This assumes interface with a UMAT call from somewhere else
//

#include <Eigen/Dense>

// 6x6 matrix of doubles
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
// 6 row col vector of doubles
typedef Eigen::Matrix<double, 6, 1> Vector6d;
// 6 row col vector of doubles
typedef Eigen::Matrix<double, 2, 1> Vector2d;

class OpenDMModel {

public:
  // Constructor
  OpenDMModel(double* props, int* nprops, double* statev, int* nstatv);

  /** @brief run the openDM model for a given UMAT call
   */
  void runModel(double* strain, double* dstrain, double* stress,
		double* statev, double* ddsdde, double* sse,
		double* spd, double* scd);

private:
  /** @brief computeS0 from UMAT array props
   */
  void computeS0(double* props);

  /** @brief unpack model parameters from props
   */
  void unpackParams(double* props);

  /** @brief unpack stateVars from UMAT call
   */
  void unpackStateVars(double* statev);

  /** @brief update statevars ater run
   */
  void updateStateVars(const Vector2d& yMax, const Matrix6d& Ceff,
		               const Vector2d& dVals, double* statev);

  /** @brief set outputs from UMAT
   */
  void setUMATOuts(const Vector6d& sig, const Vector6d& epsStar,
		   const Matrix6d& matTang, double* stress,
		   double* ddsdde, double* sse) const;

  /** @brief calculate yMax_i for given strain
   */
  Vector2d calcYVals(const Vector6d& epsStarMac) const;

  /** @brief calculate g_i for given strain
   */
  Vector2d calcGVals(const Vector2d& yMax) const;

  /** @brief calculate damage values for each damage param
   */
  Vector2d calcDVals(const Vector2d& gVals) const;

  /** @brief calculate damage effect tensor H1, H2
   */
  void calcH1H2(const Vector6d& stressEst, Matrix6d& H1,
		Matrix6d& H2) const;

  /** @brief get activation function for stress
   * basically just to affect stiffness only in tension
   */
  double stressAct(const double& stressVal) const {
    // Heaviside from stressVal
    return stressVal > 0.0 ? 1.0 : 0.0;
  }

  /** @brief compute analytical material tangent stiffness
   */
  void computeMatTang(const Matrix6d& Ceff, const Matrix6d& H1,
		      const Matrix6d& H2, const Vector6d& epsStar,
		      const Vector6d& epsStarMac, const Vector2d yMaxVals,
		      const Vector2d gVals,
		      Matrix6d& matTang) const;

  // compliance tensor
  Matrix6d S0;
  // stiffness tensor
  Matrix6d C0;
  
  // Vectors of all params for each damage var
  Vector2d hs, b, y0, yc, pe, dc;
  // vector of yMax (stateVar)
  Vector2d yMaxSave;
  // matrix of prev Ceff (stateVar)
  Matrix6d CeffOld;

};
