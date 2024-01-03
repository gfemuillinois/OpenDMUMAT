// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// definition file for OpenDM Model with 2 damage parameters
// This is derived from OpenDMModel class
// It captures damage in the 0 & 90 directions in plane of composite
//
#ifndef MODEL_OPENDM2PARAM_H
#define MODEL_OPENDM2PARAM_H

#include "model_opendm.hpp"

class OpenDMModel2Param : public OpenDMModel {

public:
  // constructor
  OpenDMModel2Param(double* props, int* nprops, double* statev, int* nstatv);

  // default destructor
  virtual ~OpenDMModel2Param() override =default;

  // delete copy constructor (should never be used)
  OpenDMModel2Param(const OpenDMModel2Param&)=delete;
  
  // delete operator= (should never be used)
  OpenDMModel2Param& operator=(const OpenDMModel2Param&)=delete;

private:

  /** @brief unpack model parameters from props
   */
  virtual void unpackParams(double* props) override;

  /** @brief Create matrix H1 & H2 for calculation
   * of Seff. NOTE this does not include stress activation
   */
  void createH1H2();

  /** @brief calculate yMax_i for given strain
   */
  virtual VectorXd calcDrivingForces(const Vector6d& epsStarMac) override;

  /** @brief calculate damage effect tensor H1, H2
   */
  virtual void computeSEff(const Vector6d& stressEst, const VectorXd& dVals,
			   Matrix6d& Seff) override;

  /** @brief compute analytical material tangent stiffness
   */
  virtual void computeMatTang(const Matrix6d& Ceff, const Vector6d& epsStar,
		      const Vector6d& epsStarMac, const VectorXd& yMaxVals,
		      const VectorXd& gVals,
		      Matrix6d& matTang) override;

  // Matrices to compute Seff
  Matrix6d H1, H2;
  
protected:
  // no protected members

};
#endif
