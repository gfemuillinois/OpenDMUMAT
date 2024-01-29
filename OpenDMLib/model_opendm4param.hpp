// - Bryce Mazurowski <brycepm2@gmail.com> -
//
// definition file for OpenDM Model with 2 damage parameters
// This is derived from OpenDMModel class
// It captures damage in the 0 & 90 directions in plane of composite
//
#ifndef MODEL_OPENDM4PARAM_H
#define MODEL_OPENDM4PARAM_H

#include "model_opendm.hpp"

class OpenDMModel4Param : public OpenDMModel {

public:
  // constructor
  OpenDMModel4Param(double* props, int* nprops, double* statev, int* nstatv);

  // default destructor
  virtual ~OpenDMModel4Param() override =default;

  // delete copy constructor (should never be used)
  OpenDMModel4Param(const OpenDMModel4Param&)=delete;
  
  // delete operator= (should never be used)
  OpenDMModel4Param& operator=(const OpenDMModel4Param&)=delete;

  /** @brief unpack stateVars from UMAT call
   */
  virtual void unpackStateVars(double* statev) override;

  /** @brief update statevars ater run
   */
  virtual void updateStateVars(const VectorXd& yMax, const Matrix6d& Ceff,
               const VectorXd& dVals, double* statev) override;

  /** @brief calculate yMax_i for given strain
   */
  virtual VectorXd calcDrivingForces(const Vector6d& epsStar,
                                     Vector6d& epsD1Plus,
                                     Vector6d& epsD2Plus) override;


private:

  /** @brief unpack model parameters from props
   */
  virtual void unpackParams(double* props) override;

  /** @brief Create transformation matrices to get S/H into
   * +/-45 degree coordinate frames and back
   */
  void createTEpsMats();

  /** @brief Create H_i matrices for computing Seff
   */
  void createHMats();
  
  /** @brief calculate positive part of D1 strains
   */
  void posPartStrainD1(const Vector6d& epsD1, Vector6d& epsD1Plus);

    /** @brief calculate positive part of D2 strains
   */
  void posPartStrainD2(const Vector6d& epsD2, Vector6d& epsD2Plus);

  /** @brief calculate derivative of epsilonD1Plus
   * wrt strain components
   */
  void calcDEpsD1PlusDEps(const Vector6d& epsD1,
                          Matrix6d& dEpsD1PlusDEps) const;

  /** @brief calculate derivative of epsilonD1Plus
   * wrt strain components
   */
  void calcDEpsD2PlusDEps(const Vector6d& epsD1,
                          Matrix6d& dEpsD1PlusDEps) const;
  
  /** @brief calculate damage effect tensor H1, H2
   */
  virtual void computeSEff(const Vector6d& stressEst, const VectorXd& dVals,
               Matrix6d& Seff) override;

  /** @brief compute analytical material tangent stiffness
   */
  virtual void computeMatTang(const Matrix6d& Ceff, const Vector6d& epsStar,
              const Vector6d& epsD1Plus, const Vector6d& epsD2Plus, 
              const VectorXd& yMaxVals,
              const VectorXd& gVals,
              Matrix6d& matTang) override;
  
  // extra vectors for damage params
  VectorXd hs4, hs5;
  
  // Transformation matrices for shear damage
  Matrix6d Teps_p45, Teps_n45;
  
  // Matrices to compute Seff
  Matrix6d H1, H2, H4, H5;

  // eVal/eVect work areas
  Matrix3d eValsD1, eVectsD1, eValsD2, eVectsD2;
  
protected:
  // no protected members

};


#endif
