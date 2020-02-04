#ifndef mvrhd_h
#define mvrhd_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <list>
#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b) 
#define IDEX(a, b) (a>b?1:0) 
using namespace Rcpp;
using namespace Eigen;

//----------------------------------------------------------------**
//***----------------------parameters for penalization---------------**
extern struct Options
{
	int r1;
	int r2;
	int r3;
	int p;
	int q;
	int K;
	int n;
	double eps;
	double eps1;
	double utol;
	double ftol;
	double Pitol;     
	double tau_min;
	double eta;
	double tiny; 
	double rhols; 
	double gamma; 
	int is_LR;	
	int max_step;
    int max_step1;	
	int onStiefel;
	int intercept;
	int pz;
}opts;

extern struct Options_pen
{
	int pen; 
	int nlam;
	int dfmax;
	int isPenU;
	int isPenColumn;
	int isFISC;
	double lam_max;
	double lam_min;
	double alpha;
	double gamma_pen;
	double lambda;
	double gamma_tanh;
	double thresh;
	double delta;
	double max_step;
	double max_step1;
}opts_pen;

//----------------------------------------------------------------**
//***----------------------parameters for penalization---------------**
//typedef struct Options opts;

//typedef struct Options_pen opts_pen;

//----------------------------------------------------------------**
//***----------------------UpTriangularInv------------------------**
MatrixXd UpTriangularInv(MatrixXd A);
//----------------------------------------------------------------**
//***---------------------- Q*R of qr decomposition --------------**
MatrixXd QbyR(MatrixXd Q, MatrixXd R, int isQR);
//----------------------------------------------------------------**
//***--------------------penalty----------------------------------**
double penalties(double z, double v, double lambda, double alpha, double gamma, int penalty);

#endif