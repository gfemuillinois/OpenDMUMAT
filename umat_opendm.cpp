// Implementation of OpenDM as a UMAT
//

#include <iostream>
#include <ctime>
#include <fstream>

// /********************************************************************/
// /********************************************************************/

extern "C" {
void runInitOpenDMModel(double* props, int* nprops, double* statev, int* nstatv,
			double* strain, double* dstrain, double* stress,
			double* ddsdde, double* sse, double* spd, double* scd);
}

/********************************************************************/
/********************************************************************/

extern "C" void umat(double *stress, double *statev, double *ddsdde, double *sse, double *spd,
		     double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
		     double *stran, double *dstran, double *time, double *dtime, double *temp,
		     double *dtemp, double *predef, double *dpred, char *cmname, int *ndi,
		     int *nshr, int *ntens, int *nstatv, double *props, int *nprops, 
		     double *coords, double *drot, double *pnewdt, double *celent, double *dfgrd0, 
		     double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, 
		     int *kstep, int *kinc) {
  // SetUp and run openDM model from UMAT call//
  runInitOpenDMModel(props, nprops, statev, nstatv, stran, dstran, stress, ddsdde, sse, spd, scd);
 
}
/********************************************************************/
/********************************************************************/

 int main() {										   
   // elastic material definitions							   
   double E11, E22, E33, nu12, nu23, nu13, G12, G23, G13;				   
   E11 = 12.416e6;									   
   E22 = E11;										   
   E33 = 7.513e6;									   
   nu12 = 0.1106;									   
   nu23 = 0.1707;									  
   nu13 = nu23;									   
   G12 = 3.332e6;									   
   G23 = 0.670e6;									   
   G13 = G23;										   
   // damage model parameters								   
   double hs1, hs2, b1, b2, y01, y02, yc1, yc2, pe1, pe2, dc1, dc2;			   
   hs1 = 2.5;										   
   hs2 = hs1;										   
   b1 = 0.2;										   
   b2 = b1;										   
   y01 = 0.005;									   
   y02 = y01;										   
   yc1 = 500;										   
   yc2 = yc1;										   
   pe1 = 1.5;										   
   pe2 = pe1;										   
   dc1 = 2.0;										   
   dc2 = dc1;										   
   											   
   // need to set a bunch of stuff							   
   double strain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};					   
   double dstrain[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};					   
   double stress[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};					   
   double statev[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};		   
   double ddsdde[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 					   
 		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 					   
 		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 					   
 		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 					   
 		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 					   
 		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 					   
   double props[] = {E11, E22, E33, nu12, nu23, nu13, G12, G23, G13,			   
 		    hs1, hs2, b1, b2, y01, y02, yc1, yc2, pe1, pe2, dc1, dc2};		   
   double sse, spd, scd, rpl;								   
   // Current deltaT for which solution is seeking					   
   double dTime = 0.001;								   
   double time[] = {0.0, 0.001};							   
  											   
   int ndi = 3; // 3 direct stress/strain components					   
   int nshr = 3; // 3 shear stress/strain components					   
   int ntens = ndi+nshr; // 6 total componenets. ndi+nshr				   
   int nstatv = 11;									   
   int nprops = 20; // Number of properties for the user mat				   
   int noel = 1;									
   int npt = 0; // Integration point number						   
   double ddsddt = 0.;								   
   double drplde = 0.;								   
   double drpldt = 0.;								   
   double temp = 0.;									      
   double dTemp = 0.;									   
   double predef = 0.; 								   
   double dpred = 0.;									   
   char cmname[80] = "OpenDM";							   
   double coords = 0.;								   
   double drot = 0.;									   
   double pnewdt = 0.;								   
   double dfgrd0 = 0.;								   
   double dfgrd1 = 0.;								   
   double celent;									   
   int layer, kspt, kstep, kinc;							   
   layer = kspt = kstep = kinc = 0;							   
 											   
   // setup a loop to call umat each step						   
   int nInc = 10;									   
   double stepSize = 10000.0/(nInc);							   
   double S11 = 1.0/E11, S12 = -nu12/E11, S22 = 1.0/E22, S44 = 1.0/G12;		   
   std::ofstream myFile;								   
   myFile.open("umatTestOut_Shear.csv");						   
   myFile << "step, strain, stress" << std::endl;					   
   myFile << 0 << ", " << 0.0 << ", " << 0.0 << std::endl;				   
   for (int iPull = 1; iPull <= nInc; iPull++) {					   
     // Increment dstrain								   
     dstrain[0] = iPull*stepSize*(S11 + S12);						   
     dstrain[1] = iPull*stepSize*(S12 + S22);						   
     dstrain[3] = iPull*stepSize*(S44);						   
     // run model									   
     clock_t start = clock();								   
     (umat)(stress, statev, ddsdde, &sse, &spd, &scd, &rpl, &ddsddt,			   
 	   &drplde, &drpldt, strain, dstrain, time, &dTime, &temp,			   
 	   &dTemp, &predef, &dpred, cmname, &ndi, &nshr, &ntens, &nstatv,		   
 	   props, &nprops, &coords, &drot, &pnewdt, &celent, &dfgrd0,			   
 	   &dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc);				   
 											   
     clock_t end = clock();								   
     std::cout << "step = " << iPull << " Time = " << (end - start) << std::endl;	   
     // std::cout << "Strain: " << (*dstrain) << std::endl;				   
     // std::cout << "Stress: " << (*stress) << std::endl;				   
     myFile << iPull << ", " << dstrain[3] << ", " << stress[3] << std::endl;		   
   }											   
   myFile.close();									   
 											   
 }											   
 /********************************************************************/		   
 /********************************************************************/		   