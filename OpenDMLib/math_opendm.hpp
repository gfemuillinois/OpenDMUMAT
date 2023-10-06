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
  // Computes B = A^{-1}
  // This code assumes we have a 6x6 matrix composed of a
  // dense 3x3 part in upper left
  // and diagonal part in lower right
  // zeros else

  double a, b, c, d, e, f, g, h, i;
  a = A(0,0);
  b = A(0,1);
  c = A(0,2);
  d = A(1,0);
  e = A(1,1);
  f = A(1,2);
  g = A(2,0);
  h = A(2,1);
  i = A(2,2);
  
  // A^{-1} = 1/det(A) * adj(A)
  // determinant
  double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  
  // adj(A) = cof(A)^{T}
  // Upper Left 
  B(0,0) = (e*i - f*h)/det;
  B(0,1) = (c*h - b*i)/det;
  B(0,2) = (b*f - c*e)/det;
  
  B(1,0) = (f*g - d*i)/det;
  B(1,1) = (a*i - c*g)/det;
  B(1,2) = (c*d - a*f)/det;
  
  B(2,0) = (d*h - e*g)/det;
  B(2,1) = (b*g - a*h)/det;
  B(2,2) = (a*e - b*d)/det;

  // Lower Right
  B(3,3) = 1.0/A(3,3);
  B(4,4) = 1.0/A(4,4);
  B(5,5) = 1.0/A(5,5);
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

inline void posPartStrain(const Vector6d& eps, Vector6d& epsPlus) {
  // TODO: Could optimize this for D1/D2 but probably need to split
  // calculate positive part of strain tensor
  // aka MacBracketPlus on eVals
  // Then transformed back to normal coord
  Matrix3d epsilon;
  epsilon(0,0) = eps(0); epsilon(1,1) = eps(1); epsilon(2,2) = eps(2);
  epsilon(0,1) = 0.5*eps(3); epsilon(1,0) = 0.5*eps(3);
  epsilon(0,2) = 0.5*eps(4); epsilon(0,2) = 0.5*eps(4);
  epsilon(1,2) = 0.5*eps(5); epsilon(1,2) = 0.5*eps(5);

  Eigen::SelfAdjointEigenSolver<Matrix3d> eigensolver(epsilon);
  if (eigensolver.info() != Eigen::Success) abort(); 

  epsilon = Matrix3d::Zero();
  for (int iEVal = 0; iEVal<3; iEVal++) {
    // According to eigen3Docs, this will be explicit form
    epsilon(iEVal,iEVal) = macaulayBracketPlus(eigensolver.eigenvalues()(iEVal));
  }
  // ASSUMING Eigenvectors come back in normalized form
  // transform epsilon back to normal coords
  Matrix3d eVects = eigensolver.eigenvectors();
  epsilon = eVects*epsilon*eVects.transpose();
  epsPlus(0) = epsilon(0,0); epsPlus(1) = epsilon(1,1); epsPlus(2) = epsilon(2,2);
  epsPlus(3) = 2.0*epsilon(0,1);
  epsPlus(4) = 2.0*epsilon(0,2);
  epsPlus(5) = 2.0*epsilon(1,2);
}
/********************************************************************/
/********************************************************************/

#endif
