// Basic math operations tailored for openDM
#ifndef MATH_OPENDM_H
#define MATH_OPENDM_H

#include <Eigen/Dense>

// 3x3 matrix of doubles
typedef Eigen::Matrix<double, 3, 3> Matrix3d;
// 6x6 matrix of doubles
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
// 6 row col vector of doubles
typedef Eigen::Matrix<double, 6, 1> Vector6d;

/********************************************************************/
/********************************************************************/

inline void matrixInverse(const Matrix6d& A, Matrix6d& B) {
  // Computes B = A^{-1} using eigen
  B = A.inverse();
}
/********************************************************************/
/********************************************************************/

inline double macaulayBracketPlus(const double& a) {
  // Calculate macaulay bracket + 
  // <a>_+ = a > 0.0 ? a : 0.0
  return a > 0.0 ? a : 0.0;
}
/********************************************************************/
/********************************************************************/

inline double macaulayBracketMinus(const double& a) {
  // Calculate macaulay bracket + 
  // <a>_+ = a > 0.0 ? a : 0.0
  return a < 0.0 ? a : 0.0;
}
/********************************************************************/
/********************************************************************/

#endif
