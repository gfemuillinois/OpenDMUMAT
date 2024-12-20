#include <iostream>

#include "model_opendm.hpp"
#include "model_opendm2param.hpp"
#include "model_opendm4param.hpp"

using std::cout;
using std::endl;

/** @brief Init OpenDM Object and run model from cpp source
 *  To facilitate simple umat file that links with sharedLib with all of this stuff
 */
extern "C" {
  void runInitOpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
			  double* strain, double* dstrain, double* stress,
			  double* ddsdde, double* sse, double* spd, double* scd) {
    // Create OpenDM model object
    // TODO: add some error catches!!
    OpenDMModel* p_openDMModel = nullptr;
    if ((*nprops) == 25 && (*nstatv) == 13 ) {
      // create 2 parameter damage model
      p_openDMModel = new OpenDMModel2Param(props, nprops, statev, nstatv);
    } else if ((*nprops) == 42 && (*nstatv) == 29) {
      // create 4 parameter model
      p_openDMModel = new OpenDMModel4Param(props, nprops, statev, nstatv);
    }
    // run openDMModel
    p_openDMModel->runModel(strain, dstrain, stress, statev, ddsdde, sse, spd, scd);
    // delete OpenDMModel object
    delete p_openDMModel;
  }
}
