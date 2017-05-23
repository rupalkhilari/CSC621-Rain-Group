// in myclass.h
#ifndef __SCCALC_HH_INCLUDED__   // if x.h hasn't been included yet...
#define __SCCALC_HH_INCLUDED__

#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>

// #include "Riostream.h"
// #include "TROOT.h"
// #include "TApplication.h"
// #include "TCanvas.h"
// #include "TH1.h"
// #include "TSystem.h"
// #include "TBrowser.h"
// #include "TFile.h"
// #include "TRandom.h"
// #include "TMultiDimFit.h"
// #include "TVectorD.h"
// #include "TMath.h"

#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <math.h>

struct SpineComponent {
  int index1;
  int index2;
  double angle;
  int apex;
};

//Spine Curvature Calculator
class ScCalc
{
private:
    double transMatrix[4][4];
    void loadSpineX(double spine[][3], unsigned length, unsigned spineNumber);
public:
    std::vector<std::vector<int> > spine1vec;
    double (*spine1)[3];
    double (*fit1)[3];
    double (*fit2)[3];
    unsigned spine1Length;
    double* spacing1;
    double* angles1;
    double* angles2;
    std::vector<SpineComponent> compX;
    std::vector<SpineComponent> compY;
    std::vector<double> kappaX;
    std::vector<double> kappaY;
    std::vector<double>  kappa;
    std::vector<double> caX;
    std::vector<double> caY;
    std::vector<double> maxAnglesX;
    std::vector<double> maxAnglesY;
    std::vector<std::vector<int>> colorMap;
    std::vector<std::string> vertebraePoints;
    std::string curveType;
    std::string sagittalModifier;
    std::string lumbarModifier;

    // TRandom* gRandom;
    void printVector(std::vector<std::vector<int> > vec);
    void saveVector(std::vector<std::vector<int> > vec, char* fileName);
    std::vector<std::vector<int> > loadVector(char* fileName);
    void loadSpine1(double spine[][3], unsigned length);
    void loadSpine1(std::vector<std::vector<int>> spine);
    void loadTransofrm(double matrix[4][4]);
    void transformSpine1();
    void compareSpines();
    void printSpine1();
    void printTransform();
    void printAngles();
	//TODO: add mean square distance function
    //TODO: curve fit function

	//TODO: add circular calculation funcion
	//TODO: add dot product over determinant product funcion
	//TODO: add smoothing function
	static double firstDerivative(double currentData[5]);
	static double curveDerivative(double points[5][3]);
	static double curveSecDerivative(double points[9][3]);

    void crateSpineFit(int spineNum, unsigned order = 7);
    bool polynomialfit(int obs, int degree, double *dx, 
        double *dy, double *store);

    void anglesX(double points[][3], double angles[], int npoints);
    void anglesY(double points[][3], double angles[], int npoints);
    void maXanglesX(double points[][3], double angles[], int npoints);
    void getMax3Dangles(double points[][3], double angles[], int npoints);

    // Int_t multidimfit(bool doFit = true);
    // int CompareResults(TMultiDimFit *fit, bool doFit);
    // void makeData(Double_t* x, Double_t& d, Double_t& e);

    // added this
    std::vector<std::vector<int>> loadAnnotationData(char* fileName);
    void loadColorMap(char* fileName);
    std::vector<SpineComponent> getMax3DanglesMod(double points[][3], double angles[], int npoints);
    void analyze3DAngles(std::vector<SpineComponent> comp, int type);
    void crateSpineFitX(int spineNum, unsigned order = 7);
    void crateSpineFitY(int spineNum, unsigned order = 7);
    double getCurvatureAngle(double x1, double y1, double z1, double x2, double y2, double z2, double kappa);
    void getCurvatureAngleMod(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz, int type);
    void doLenkeClassification(std::vector<SpineComponent> comp, std::vector<SpineComponent> results);
    void doLenkeClassificationReloaded(std::vector<SpineComponent> comp, std::vector<SpineComponent> results, int type);
    void doSimpleClassification(std::vector<SpineComponent> comp, std::vector<SpineComponent> results);
    void computeGeometricCurvature2D(double points[][3], unsigned spLength, unsigned order, double *xStore, double *x);
    void getGeometricCurvature(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz);
    void getGeometricCurvatureMod(double points[][3], unsigned spLength, unsigned order, double *xStore, double *yStore, double *dz, int type);
};

#endif
