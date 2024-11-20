// Implementation of OpenDM as a UMAT
// by: Bryce Mazurowski (brycepm2@gmail.com)
//

/** @brief interface function to C++ UMAT library
 */
extern "C" {
void runInitOpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
                        double* strain, double* dstrain, double* stress,
                        double* ddsdde, double* sse, double* spd, double* scd);
}

/********************************************************************/
/********************************************************************/
/** @brief Abaqus UMAT interface 
 */

extern "C" void umat(double *stress, double *statev, double *ddsdde,
                     double *sse, double *spd, 
                     double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
                     double *stran, double *dstran, double *time, double *dtime, double *temp,
                     double *dtemp, double *predef, double *dpred, char *cmname, int *ndi,
                     int *nshr, int *ntens, int *nstatv, double *props, int *nprops,
                     double *coords, double *drot, double *pnewdt, double *celent, double *dfgrd0,
                     double *dfgrd1, int *noel, int *npt, int *layer, int *kspt,
                     int *kstep, int *kinc) {
  // SetUp and run openDM model from UMAT call
  runInitOpenDMModel(props, nprops, statev, nstatv, stran, dstran,
                     stress, ddsdde, sse, spd, scd);
}
/********************************************************************/
/********************************************************************/
