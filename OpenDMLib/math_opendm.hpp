// Basic math operations tailored for openDM
#ifndef MATH_OPENDM_H
#define MATH_OPENDM_H

#include <Eigen/Dense>

/********************************************************************/
/********************************************************************/

void matrixInverse(const Matrix6d& A, Matrix6d& B) {
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
  
  double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  
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

double macaulayBracket(const double& a) {
  // Calculate macaulay bracket + 
  // <a>_+ = a > 0.0 ? a : 0.0
  double b = a > 0.0 ? a : 0.0;
  return b;
}
/********************************************************************/
/********************************************************************/
#endif
