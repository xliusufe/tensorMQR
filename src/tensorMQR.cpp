//[[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <list>
#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b) 
#define IDEX(a, b) (a>b?1:0) 
using namespace Rcpp;
using namespace Eigen;

struct Options
{
	int r1;
	int r2;
	int r3;
	int p;
	int q;
	int n;
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
}opts;

struct Options_pen
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
//***----------------------cbind----------------------------------**
MatrixXd cbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n = A.rows(), p1 = A.cols(), p2 = B.cols();
	MatrixXd C = MatrixXd::Constant(n, p1 + p2, 0);
	C.block(0, 0, n, p1) = A;
	C.block(0, p1, n, p2) = B;
	return C;
}
//----------------------------------------------------------------**
//***----------------------rbind----------------------------------**
MatrixXd rbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n1 = A.rows(), n2 = B.rows(), p = A.cols();
	MatrixXd C = MatrixXd::Constant(n1 + n2, p, 0);
	C.block(0, 0, n1, p) = A;
	C.block(n1, 0, n2, p) = B;
	return C;
}

//----------------------------------------------------------------**
//***----------------------UpTriangularInv------------------------**
MatrixXd UpTriangularInv(MatrixXd A)
{
	int i, j, k,n=A.cols();
	MatrixXd B = MatrixXd::Constant(n, n, 0);
	for (i = 0; i < n; i++) B(i, i) = 1;
	for (i = n - 1; i >= 0; i--)//rows
	{
		if (A(i, i) != 1)
			for (j = i; j<n; j++)
				B(i, j) = B(i, j) / A(i, i);
		if (i>0)
		{
			for (j = i; j<n; j++)// columns
				for (k = 0; k<i; k++)// rows
					B(k, j) = B(k, j) - A(k, i) * B(i, j);
		}
	}
	return B;
}

//----------------------------------------------------------------**
//***---------------------- Q*R of qr decomposition --------------**
MatrixXd QbyR(MatrixXd Q, MatrixXd R, int isQR)
{
	//isQR=1 denotes Q*R; otherwise R*Q
	int i, j, k, p = R.cols(),n;
	double temp1;
	MatrixXd A = Q;
	if(isQR){
		n = Q.rows();
		for (k = 0; k < p; k++)
			for (j = 0; j < n; j++) {
				temp1 = 0;
				for (i = 0; i <= k; i++) temp1 += Q(j, i)*R(i, k);
				A(j,k) = temp1;
			}
	}
	else{
		n = Q.cols();
		for (k = 0; k < n; k++)
			for (j = 0; j < p; j++) {
				temp1 = 0;
				for (i = j; i < p; i++) temp1 += R(j, i) * Q(i, k);
				A(j,k) = temp1;
			}
	}
	return A;
}

//----------------------------------------------------------------**
//***---------------------- R*R of Upper triangle ----------------**
MatrixXd tRbyR(MatrixXd R)
{
  int i, j, k, ii, p;
  double temp1;
  MatrixXd A = R;
  p = R.cols();
  
  for (k = 0; k < p; k++) {
    for (j = 0; j < p; j++) {
      temp1 = 0;
	  ii = k < j ? k : j;
      for (i = 0; i <= ii; i++) temp1 += R(i, j) * R(i, k);
      A(j,k) = temp1;
    }
  }
  return A;
}

//----------------------------------------------------------------**
//***----------------------extract columns------------------------**
MatrixXd extractRows(MatrixXd Z, int K, VectorXi b)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>n) stop("The length of index b must be not greater than the number of columns of Z!");
	MatrixXd C = MatrixXd::Constant(K, p, 0);
	for(i=0;i<K;i++) C.row(i) = Z.row(b[i]);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractCols(MatrixXd Z, int K, VectorXi b)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>p) stop("The length of index b must be not greater than the number of columns of Z!");
	MatrixXd C = MatrixXd::Constant(n, K, 0);
	for(i=0;i<K;i++) C.col(i) = Z.col(b[i]);
	return C;
}

//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractCols_by(MatrixXd Z, int min, int K, int by)
{
	int n = Z.rows(), p = Z.cols(), i;
	if(K>p) stop("The length of index b must be not greater than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(n, K, 0);
	for(i=0;i<K;i++) C.col(i) = Z.col(min+by*i);
	return C;
}
//----------------------------------------------------------------**
//***--------------------transfer modal of unfoldings-------------**
// [[Rcpp::export]]
MatrixXd TransferModalUnfoldings(MatrixXd S, int d1, int d2, int r1, int r2, int r3)
{
  int j;
  MatrixXd S1,S_3;
  if (d1 == 3) {
    if (d2 == 1){
      S1 = S.row(0).transpose();
      S1.resize(r1, r2); 
      for (j = 1; j < r3;j++) {
        S_3 = S.row(j).transpose();
        S_3.resize(r1, r2);
        S1 = cbind_rcpp(S1, S_3);
      }
    }	
    if (d2 == 2) {
      S1 = S.block(0, 0, r3, r1).transpose();
      S1.resize(1,r1*r3);
      for (j = 1; j < r2; j++) {
        S_3 = S.block(0, j*r1, r3, r1).transpose();
        S_3.resize(1,r1*r3);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  
  if (d1 == 2) {
    if (d2 == 1) {
      S1 = S.block(0, 0, r2, r1).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1,S.block(0,j*r1,r2,r1).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r2, r1).transpose();
      S1.resize(1, r2*r1);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r1, r2, r1).transpose();
        S_3.resize(1, r2*r1);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
    
  }
  
  if (d1 == 1) {
    if (d2 == 2) {
      S1 = S.block(0, 0, r1, r2).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1, S.block(0, j*r2, r1, r2).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r1, r2);
      S1.resize(1, r1*r2);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r2, r1, r2);
        S_3.resize(1, r1*r2);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  return S1;
}

//----------------------------------------------------------------**
//***----------------------uppertriangleNorm1 of A----------------**
double uppertriangleNorm1(MatrixXd A)
{
	int j, k, p;
	double norm1 = 0, temp1;
	p = A.cols();

	for (k = 0; k < p; k++) {
		temp1 = 0;
		for (j = 0; j <= k; j++) {
			temp1 += fabs(A(j, k));
		}
		if (temp1 > norm1) norm1 = temp1;
	}
	return norm1;
}

//----------------------------------------------------------------**
//***------QRcondition_number for symetric matrix-----------------**
double condition_numberQRSym(MatrixXd R)
{
  return uppertriangleNorm1(R) * uppertriangleNorm1(UpTriangularInv(R));
}

//----------------------------------------------------------------**
//***----------------solve linear system by QR--------------------**
VectorXd solveEquationQR(MatrixXd A, VectorXd b)
{
  //solve linear system Ax = b, A = A.transpose();
  int p = b.size();
  HouseholderQR<MatrixXd> qr;
  qr.compute(A);
  MatrixXd R, Q;
  Q = qr.householderQ();
  R = qr.matrixQR().triangularView<Upper>();
  VectorXd X;
  
  if (condition_numberQRSym(R) > 1e10){
    MatrixXd temp, IDEN = MatrixXd::Identity(p, p);
    temp = A + (IDEN.array()*1e-4).matrix();
    X = temp.colPivHouseholderQr().solve(b);
  }
  else X = QbyR(Q.transpose(),UpTriangularInv(R),0)*b;
  return X;
}

//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
double dotproduct(VectorXd a, VectorXd b, int n){
	double out=0.0;
	for(int i=0;i<n;i++) out += a[i]*b[i];
	return out;
}

//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
VectorXd SEQ(double min, double max, double by)
{
	int leng = static_cast<int>(floor((max - min) / by + pow(3, -12)) + 1);
	VectorXd C = VectorXd::Constant(leng, 0);
	for (int i = 0; i < leng; i++)
		C[i] = min + by * i;
	return C;
}

//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
VectorXd SEQN(int min, int n, int by)
{
	// by is the increment with by >0	
	VectorXd C = VectorXd::Constant(n,0);
	for (int i = 0; i < n; i++)	C[i] = min + by * i;
	return C;
}

//----------------------------------------------------------------**
//***----------------------produce Z -----------------------------**
// [[Rcpp::export]]
MatrixXd produceZ(MatrixXd X){
	int n=X.rows(), p=X.cols();
	MatrixXd Z  = MatrixXd::Constant(n, p*p, 0);
	for (int i = 0; i < n; i++) Z.row(i) = kroneckerProduct(X.row(i), X.row(i));
	return Z;
}

//----------------------------------------------------------------**
//***----------------------extract subvector----------------------**
VectorXd subvector_extract(VectorXd U, int min, int n, int by){
	VectorXd V = VectorXd::Constant(n,0);
	for(int i=0; i<n; i++) V[i] = U[min+by*i];
	return V;
}

//----------------------------------------------------------------**
//***----------------------update subvector-----------------------**
void subvector_updated(VectorXd &U, VectorXd V, int min, int n, int by){
	for(int i=0; i<n; i++) U[min+by*i] = V[i];
}

//----------------------------------------------------------------**
//***----------------------update subvector-----------------------**
MatrixXd myQR(MatrixXd U, int k){
	int i, n=U.rows(), flag=0;
	MatrixXd R, Q, IDEN = MatrixXd::Identity(k, k);
	HouseholderQR<MatrixXd> qr;
	qr.compute(U);
	Q = (qr.householderQ());
	Q = Q.block(0,0,n,k);
	R = qr.matrixQR().triangularView<Upper>();
	for(i=0;i<k;i++) if(R(i,i)<0) {IDEN(i,i) = -1;	flag=1;}
	if(flag) Q = Q*IDEN;
	return Q;
}

//----------------------------------------------------------------**
//***--------------------updateC----------------------------------**
List updateC(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int q = opts.q, r3=opts.r3, j,kp;
	MatrixXd ztilde = Z * kroneckerProduct(B, A);
	MatrixXd StZ = ztilde * S.transpose();
	MatrixXd Cnew = MatrixXd::Constant(q, r3, 0);
	HouseholderQR<MatrixXd> qr;
	qr.compute(StZ);
	MatrixXd R = qr.matrixQR().triangularView<Upper>();
	MatrixXd Q = qr.householderQ();

	
	kp = StZ.cols();
	MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
	if (pow(condition_numberQRSym(R),2) > 1e10){
	  temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
	  for (j = 0; j < q; j++) Cnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
	}
	else
	  for (j = 0; j < q; j++) Cnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();
	
	qr.compute(Cnew);
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	return List::create(Named("Bnew") = VQ.block(0, 0, q, r3), Named("Snew") = QbyR(S, VR.block(0, 0, r3, r3), 0));
}

//----------------------------------------------------------------**
//***--------------------updateS----------------------------------**
MatrixXd updateS0(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C)
{
	int r1=opts.r1, r2=opts.r2, r3=opts.r3, d = r1*r2*r3, k, k1, j, j1;
	MatrixXd ztilde = Z * kroneckerProduct(A, B), vectorS;
	VectorXd U;
	U.setZero(d);
	MatrixXd  V = MatrixXd::Constant(d, d, 0);

	for (k = 0; k < r1*r2; k++) {
		for (j = 0; j < r3; j++) {
			U[k*r3 + j] = ztilde.col(k).transpose()*Y*C.col(j);
			for (k1 = 0; k1 < r1*r2; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(k*r3 + j, k1*r3 + j1) = kroneckerProduct(
						ztilde.col(k1).array()*ztilde.col(k).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
			}
		}
	}
	vectorS = solveEquationQR(V,U);
	vectorS.resize(r3, r1*r2);
	return vectorS;
}

//----------------------------------------------------------------**
//***--------------------updateA----------------------------------**
List updateA(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = opts.r1,r2 = opts.r2, p = opts.p, d = r1 * p, t1,t2,t3,t4,j;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zbw1, zbw2, vectorA;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);

	for (t2 = 0; t2<p; t2++) {
		for (t1 = 0; t1<r1; t1++) {
			Wt1 = W.col(t1);
			for (j = 1; j<r2; j++) Wt1 = cbind_rcpp(Wt1, W.col(j*r1 + t1));
			Zt2 = Z.col(t2);
			for (j = 1; j<p; j++)	Zt2 = cbind_rcpp(Zt2, Z.col(j*p + t2));
			zbw1 = Zt2 * B*(Wt1.transpose());
			tU[t2*r1 + t1] = (Y.array()*zbw1.array()).sum();

			for (t4 = 0; t4<p; t4++) {
				for (t3 = 0; t3<r1; t3++) {
					Wt3 = W.col(t3);
					for (j = 1; j<r2; j++)	Wt3 = cbind_rcpp(Wt3, W.col(j*r1 + t3));
					Zt4 = Z.col(t4);
					for (j = 1; j<p; j++) Zt4 = cbind_rcpp(Zt4, Z.col(j*p + t4));
					zbw2 = Zt4 * B*(Wt3.transpose());
					tV(t2*r1 + t1, t4*r1 + t3) = (zbw1.array()*zbw2.array()).sum();
				}
			}
		}
	}
	vectorA = solveEquationQR(tV,tU);
	vectorA.resize(r1, p);
	
	HouseholderQR<MatrixXd> qr;
	qr.compute(vectorA.transpose());
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	MatrixXd S1new, S1 = TransferModalUnfoldings(S,3,1,r1,r2,opts.r3);
	S1new = QbyR(S1, VR.block(0, 0, r1, r1), 0);
	return List::create(Named("Anew") = VQ.block(0, 0, p, r1), Named("Snew") = TransferModalUnfoldings(S1new,1,3,r1,r2,opts.r3));
}

//----------------------------------------------------------------**
//***--------------------updateB----------------------------------**
List updateB(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = opts.r1,r2 = opts.r2, p = opts.p, d = r2*p, q=opts.q, n=opts.n, t1,t2,t3,t4;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);

	for (t2 = 0; t2<p; t2++) {
		for (t1 = 0; t1<r2; t1++) {
			Wt1 = W.block(0, t1*r1, q, r1);  
			Zt2 = Z.block(0, t2*p, n, p);  
			zaw1 = Zt2 * A*(Wt1.transpose()); 
			tU[t2*r2 + t1] = (Y.array()*zaw1.array()).sum();
			for (t4 = 0; t4<p; t4++) {
				for (t3 = 0; t3<r2; t3++) {
					Wt3 = W.block(0, t3*r1, q, r1); 
					Zt4 = Z.block(0, t4*p, n, p); 
					zaw2 = Zt4 * A*(Wt3.transpose()); 
					tV(t2*r2 + t1, t4*r2 + t3) = (zaw1.array()*zaw2.array()).sum();
				}
			}
		}
	}
	vectorB = solveEquationQR(tV, tU);
	vectorB.resize(r2, p);
	
	HouseholderQR<MatrixXd> qr;
	qr.compute(vectorB.transpose());
	MatrixXd VR = qr.matrixQR().triangularView<Upper>();
	MatrixXd VQ = qr.householderQ();
	MatrixXd S2new, S2 = TransferModalUnfoldings(S,3,2,r1,r2,opts.r3);
	S2new = QbyR(S2, VR.block(0, 0, r1, r1), 0);
	return List::create(Named("Bnew") = VQ.block(0, 0, p, r1), Named("Snew") = TransferModalUnfoldings(S2new,2,3,r1,r2,opts.r3));
}

//----------------------------------------------------------------**
//***--------------------Estimate initial without penalty---------------**
List Estimation_nonsym(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	double  likhd0 = pow(10, 6), likhd1 = 0;
	MatrixXd Dnew,Cnew,Anew,Bnew,Snew;
	VectorXd convergence1;
	List fit;
	int step = 0;
	while (step<opts.max_step1) {
		convergence1 = VectorXd::Constant(4, 1);
		step = step + 1;
		Snew = updateS0(Y, Z, A, B, C);
		Dnew = C * Snew*kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = pow((Y - Z * Dnew.transpose()).norm(),2);
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;
		fit = updateC(Y, Z, A, B, C, S);
		Cnew = fit[0];
		Snew = fit[1];
		Dnew = Cnew * Snew *kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = pow((Y - Z * Dnew.transpose()).norm(),2);

		if (likhd1<likhd0) {
			C = Cnew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;

		fit = updateA(Y, Z, A, B, C, S);
		Anew = fit[0];
		Snew = fit[1];
		Dnew = C * Snew *kroneckerProduct(B.transpose(), Anew.transpose());
		likhd1 = pow((Y - Z * Dnew.transpose()).norm(),2);
		if (likhd1<likhd0) {
			A = Anew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[2]=0;

		fit = updateB(Y, Z, A, B, C, S);
		Bnew = fit[0];
		Snew = fit[1];
		Dnew = C * Snew *kroneckerProduct(Bnew.transpose(), A.transpose());
		likhd1 = pow((Y - Z * Dnew.transpose()).norm(),2);
		if (likhd1<likhd0) {
			B = Bnew;
			S = Snew;
			if ((likhd0 - likhd1) / likhd0<opts.ftol) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[3]=0;
		if(convergence1.sum()==0) break;
	}
	return List::create(Named("likhd") = likhd0, Named("A") = A, Named("B") = B, Named("C") = C, Named("S") = S);
}

//----------------------------------------------------------------**
//***--------------------updateV----------------------------------**
List updateV(MatrixXd Y, MatrixXd Z, MatrixXd U, MatrixXd S)
{
	int q = opts.q,j,kp;
	MatrixXd ztilde = Z * kroneckerProduct(U, U);
	MatrixXd StZ = ztilde * S.transpose();
	MatrixXd Vnew = MatrixXd::Constant(q, opts.r3, 0), VR, VQ;
	HouseholderQR<MatrixXd> qr;
	qr.compute(StZ);
	MatrixXd R = qr.matrixQR().triangularView<Upper>();
	MatrixXd Q = qr.householderQ();

	
	kp = StZ.cols();
	MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
	if (pow(condition_numberQRSym(R),2) > 1e10){
	  temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
	  for (j = 0; j < q; j++) Vnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
	}
	else
	  for (j = 0; j < q; j++) Vnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();	
	qr.compute(Vnew);
	VR = qr.matrixQR().triangularView<Upper>();
	VQ = qr.householderQ();
	return List::create(Named("Vnew") = VQ.block(0, 0, q, opts.r3), Named("Snew") = QbyR(S, VR.block(0, 0, opts.r3, opts.r3), 0));
}

//----------------------------------------------------------------**
//***--------------------updateS----------------------------------**
MatrixXd updateS(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd C)
{
	int r1=opts.r1, r2=opts.r2, r3=opts.r3, d = r1*(r1+1)*r3/2, d2=r1*(r1-1)/2, k, k1, j, j1, t,t1;
	MatrixXd ztilde = Z * kroneckerProduct(A, A), V, vectorS, Snew, Snew2;
	VectorXd U, indA, indB;
	U.setZero(d); 
	indA.setZero(r1);
	indB.setZero(d2);
	V = MatrixXd::Constant(d, d, 0);
	Snew2 = MatrixXd::Constant(r2, r1*r3, 0);

	t=0; for(k=1;k<r1;k++) for(k1=1;k1<r1-k+1;k1++) indB[t++]=(k-1)*(r1+1)+k1;
	for(k=0;k<r1;k++) indA[k]=(r1+1)*k;

	for (k = 0; k < r1; k++) {
		for (j = 0; j < r3; j++) {
			U[k*r3 + j] = ztilde.col(indA[k]).transpose()*Y*C.col(j);
			for (k1 = 0; k1 < r1; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(k*r3 + j, k1*r3 + j1) = kroneckerProduct(
						ztilde.col(indA[k1]).array()*ztilde.col(indA[k]).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
			}
		}
	}
	
    t=r1;
	for (k = 0; k < d2; k++) {
		for (j = 0; j < r3; j++) {
			U[t*r3 + j] = 2*ztilde.col(indB[k]).transpose()*Y*C.col(j);
			t1=r1;
			for (k1 = 0; k1 < d2; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(t*r3 + j, t1*r3 + j1) = 4*kroneckerProduct(
						ztilde.col(indB[k1]).array()*ztilde.col(indB[k]).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
				t1++;
			}
		}
		t++;
	}	
	
	for (k = 0; k < r1; k++) {
		for (j = 0; j < r3; j++) {
			t1=r1;
			for (k1 = 0; k1 < d2; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(k*r3 + j, t1*r3 + j1) =V(t1*r3 + j1, k*r3 + j) = 2*kroneckerProduct(
						ztilde.col(indB[k1]).array()*ztilde.col(indA[k]).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
				t1++;
			}
		}
		t++;
	}	
	vectorS = V.colPivHouseholderQr().solve(U);
	vectorS.resize(r3, r1*(r1+1)/2);
	
	for(j=0;j<r3;j++){
		t=0; 		
		for(k=0;k<r1;k++) Snew2(k,k+j*r1) = vectorS(j,t++); 	
		for(k=1;k<r1;k++) 
			for(k1=k;k1<r1;k1++){ 
				Snew2(k-1,k1+j*r1) = vectorS(j,t); 
				Snew2(k1,k-1+j*r1) = vectorS(j,t++);
			}	
	}	
	Snew = TransferModalUnfoldings(Snew2, 2, 3, r1, r2, r3);	
	return Snew;
}

double innerProduct(VectorXd A, VectorXd B, int n)
{
	double out=0.0;
	for(int i=0;i<n;i++) out += A[i]*B[i];
	return out;
}

double innerProductMatrix(MatrixXd A, MatrixXd B)
{
	return (A.array()*B.array()).sum();
}

void derivF(MatrixXd &Gamma, MatrixXd Y, MatrixXd Z, MatrixXd U, MatrixXd alpha)
{
	int k, t, m, n = Y.rows(), r1=U.cols(), p=U.rows();
	VectorXd temp,temp1,temp2;
	Gamma=MatrixXd::Constant(p, r1, 0);
	for(m = 0; m < n; m++){
		for (k = 0; k < p; k++){
			for(t = 0; t < r1; t++){
				temp = temp1 = VectorXd::Constant(r1*r1, 0);
				temp2 = U.transpose()*subvector_extract(Z.row(m),k*p,p,1);
				subvector_updated(temp, temp2, t*r1,r1,1);
				subvector_updated(temp1, temp2, t,r1,r1);
				Gamma(k, t) -= innerProduct(alpha.row(m),temp+temp1,r1*r1);			
			}
		}
	}
}

double lossFunction(MatrixXd Y, MatrixXd Z, MatrixXd VS, MatrixXd U, MatrixXd VSUU, int isPen, double lambda, int isvsuu)
{
	double likhd;
	if(!isvsuu) VSUU = kroneckerProduct(U, U)*VS.transpose();
	if(!isPen) likhd = (Y - Z*VSUU).squaredNorm()/2;
	else{
		int k, p=opts.p;
		double pen_norm=0, Unorm0, gamma_pen=opts_pen.gamma_pen;
		VectorXd Unorm = VectorXd::Constant(p, 0);
		if(opts_pen.isPenU){
			for(k = 1; k < p; k++) Unorm[k] = U.row(k).norm();
		}
		else{
			int j,m;
			MatrixXd djU;
			VectorXi ind1=VectorXi::Constant(2*p-1, 0);			
			for (j = 1; j < p; j++){	        
				for(m = 0; m < j; m++)	ind1[m] = m*p+j;
				for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
				for(m = j+1; m < p; m++)	ind1[m+p-1] = m*p+j;
				djU = extractCols(Z,2*p-1,ind1)*extractRows(VSUU,2*p-1,ind1);	
				Unorm[j] = djU.norm();
			}		
		}
		if(opts_pen.pen==1) for(k = 1; k < p; k++) pen_norm += lambda*Unorm[k];
		if(opts_pen.pen==2) for(k = 1; k < p; k++){
			Unorm0 = Unorm[k];		
			pen_norm += (lambda*Unorm0-Unorm0*Unorm0/gamma_pen/2)*IDEX(gamma_pen*lambda,Unorm0)+(1-IDEX(gamma_pen*lambda,Unorm0))*gamma_pen*lambda*lambda/2;
		}
		if(opts_pen.pen==3) for(k = 1; k < p; k++){ 
			Unorm0 = Unorm[k];
			pen_norm += lambda*Unorm0*IDEX(lambda,Unorm0)+ (1-IDEX(lambda,Unorm0))*IDEX(gamma_pen*lambda,Unorm0)*(gamma_pen*lambda*Unorm0-(Unorm0*Unorm0+lambda*lambda)/2)/(gamma_pen-1);
			pen_norm += (1-IDEX(gamma_pen*lambda,Unorm0))*lambda*lambda*(gamma_pen*gamma_pen-1)/(gamma_pen-1)/2;
		}	
		likhd = (Y - Z*VSUU).squaredNorm()/2 + opts.n*pen_norm;
	}
	return likhd;
}
// ==================================================================================
// ================= Curvilinear searching algorithm ================================
double Curvilinear(MatrixXd &Unew, MatrixXd &VSUU, double &likhd1, double tau, double Ck, double deriv, 
	MatrixXd Y, MatrixXd Z, MatrixXd Gamma, MatrixXd VS, MatrixXd U, int isPen, double lambda)
{
/*
	Input:	
	Ck has formula C_k+1 = (eta*Qk*Ck+F(U^(k+1)))/Q_(k+1)
	deriv is the derivetive of F(U_k)
	W has formula W=0.5*(H*U' - U*H')
	BA has formula BA = B'*A, where A =  cbind(H,U), B = cbind(U,-H), and W = A*B'
	BU has formula BU = B'*U
	VS has formula VS = V*S3
	(p,r1) is the dimension of U	
	opts includes tau,eta,is_LR, where
	     	tau is initial value of tau, default is 1e-3
			eta is a constant multiply of deriv, default is 0.1
			is_LR logical, 1 for W being the outer product of two low-rank matrices; 0 otherwise
	
	Output:
	tau is the optimal tau
*/
	int r1=opts.r1;	
    MatrixXd IDENr1, IDENp, IW,W1,W,A,B,BA,BU;
	IDENr1 = MatrixXd::Identity(2*r1, 2*r1);
	IDENp = MatrixXd::Identity(opts.p, opts.p);
	
	W1 = Gamma*U.transpose();
	W = 0.5*(W1 - W1.transpose());
	A = cbind_rcpp(Gamma,U);
	B = cbind_rcpp(U,-Gamma);		
	BA = B.transpose()*A;
	BU = B.transpose()*U;
	int step=0;
	while(step<5){
		step++;
		if(opts.is_LR){
			IW = IDENr1 + 0.5*tau*BA;			
			Unew = U - tau*A*(IW.colPivHouseholderQr().solve(BU));
		}
		else{
			IW = IDENp + tau*W;
			Unew = IW.colPivHouseholderQr().solve(U-tau*W*U);
		}
		if((Unew.transpose()*Unew-IDENr1.block(0,0,r1,r1)).norm()>opts.tiny)
			Unew = myQR(Unew,r1);
		VSUU = kroneckerProduct(Unew, Unew)*VS.transpose();
		likhd1 = lossFunction(Y,Z,VS,U,VSUU,isPen,lambda,1);
		if(likhd1 <= Ck - tau*deriv) break;
		tau = opts.eta*tau;  
	}
	return tau;	
}

// ==================================================================================
// ======================= updateU ======================================================
List updateU(MatrixXd Y, MatrixXd Z, MatrixXd VS, MatrixXd U0)
{			
	int step=0, p=opts.p;
	double normPi, deriv, SY, Ck, tau=opts.tau_min;
	double  likhd0 = pow(10, 6), likhd1 = 0, Q=1, Qp;

	MatrixXd Gamma, VSUU, alpha;
	MatrixXd U = U0, Old_U = U, dtX, dtXP,dtXdiff, Udiff;
	VectorXd convergence1;
	
	Gamma=MatrixXd::Constant(p, opts.r1, 0);
	VSUU = kroneckerProduct(U, U)*VS.transpose();
	alpha = (Y - Z*VSUU)*VS;		
	derivF(Gamma,Y,Z,U,alpha);
	dtXP = Gamma - U*Gamma.transpose()*U;
	normPi=dtXP.norm();

	Ck = (Y - Z*VSUU).squaredNorm()/2;
	while (step<opts.max_step){
		convergence1 = VectorXd::Constant(2, 1);
		step++;
		deriv=opts.rhols*normPi*normPi;	
		tau = Curvilinear(U, VSUU, likhd1, tau, Ck, deriv, Y, Z, Gamma, VS, Old_U, 0, 0);				
		alpha = (Y - Z*VSUU)*VS;
		derivF(Gamma,Y,Z,U,alpha);
		dtX = Gamma - U*Gamma.transpose()*U;
		normPi=dtX.norm();			
		dtXdiff = dtX - dtXP;     
		Udiff = U-Old_U;
		SY = fabs(innerProductMatrix(Udiff,dtXdiff));
		if(step%2) tau  = SY/dtXdiff.squaredNorm();
		else tau = Udiff.squaredNorm()/SY; 
		tau = MAX(MIN(tau, 1e20), 1e-20);
		if(fabs(likhd1-likhd0)/fabs(likhd0+1)<opts.ftol && Udiff.norm()/sqrt(p) < opts.utol) break;
		if(normPi<opts.Pitol) break;
		if (likhd1>likhd0) break;
		Old_U = U;
		likhd0 = likhd1;
		dtXP = dtX;
        Qp = Q; 
	    Q = opts.gamma*Qp + 1; 
	    Ck = (opts.gamma*Qp*Ck + likhd1)/Q;		
	}
	return List::create(Named("likhd") = likhd1, Named("U") = U, Named("Unorm") = (U-Old_U).norm(), Named("step")=step);
}

//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
// [[Rcpp::export]]
List Estimation(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd U, MatrixXd V, List optsList)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	U is p*r1 matrix
	V is q*r3 matrix
	S is r3*(r1*r2) matrix

	Output:
	Dnew is a estimator of D(3) = V%*%S%*%kronecker(t(U),t(U))
	*/	
	opts.utol = as<double>(optsList[0]);
	opts.ftol = as<double>(optsList[1]);
	opts.Pitol = as<double>(optsList[2]);
	opts.tau_min = as<double>(optsList[3]);
	opts.eta = as<double>(optsList[4]);
	opts.tiny = as<double>(optsList[5]);
	opts.gamma = as<double>(optsList[6]);
	opts.rhols = as<double>(optsList[7]);
	opts.max_step = as<int>(optsList[8]);
	opts.max_step1 = as<int>(optsList[9]);
	opts.is_LR = as<int>(optsList[10]);
	opts.n = as<int>(optsList[11]);
	opts.r1 = as<int>(optsList[12]);
	opts.r2 = as<int>(optsList[13]);
	opts.r3 = as<int>(optsList[14]);
	opts.p = as<int>(optsList[15]);
	opts.q = as<int>(optsList[16]);
	
	
	double  likhd0 = pow(10, 6), likhd1 = 0;
	MatrixXd Dnew, Unew, Vnew, Snew, Z = produceZ(X);
	VectorXd convergence1;	
	List fit;
	int step = 0;
	while (step<opts.max_step1) {
		convergence1 = VectorXd::Constant(3, 1);
		step++;
		
		Snew = updateS(Y, Z, U, V);
		Dnew = V*Snew*kroneckerProduct(U.transpose(), U.transpose());
		likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;	
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;		
		fit = updateV(Y, Z, U, S);
		Vnew = fit[0];
		Snew = fit[1];
		Dnew = Vnew*Snew*kroneckerProduct(U.transpose(), U.transpose());
		likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;
		if (likhd1<likhd0) {
			V = Vnew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;
		fit = updateU(Y, Z, V*S, U);	
        Unew = fit[1];		
		Dnew = V*S*kroneckerProduct(Unew.transpose(), Unew.transpose());
		likhd1 = fit[0];
		if (likhd1<likhd0) {
			U = Unew;
			if ((likhd0 - likhd1) / likhd0<opts.ftol) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[2]=0;
		if(convergence1.sum()==0) break;
	}
	Dnew = V*S*kroneckerProduct(U.transpose(), U.transpose());
	return List::create(Named("likhd") = likhd0, Named("Dnew") = Dnew, Named("S") = S, Named("U")=U, Named("V")=V);
}

//----------------------------------------------------------------**
//***--------------------setup tuning parameters------------------**
// [[Rcpp::export]]
VectorXd setuplambda(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd U, MatrixXd V, int isPenU, int nlam, VectorXd setlam)
{
	int n=Y.rows(), p = U.rows(), r1 = U.cols(), j;
	double lam_max, lam_min, alpha, max_lam, max_tmp=0;
	VectorXd lambda, lambda1;	
	MatrixXd Z = produceZ(X);
	if(isPenU){
		MatrixXd VS = V*S, Gamma=MatrixXd::Constant(p, U.cols(), 0), alpha1;	
		alpha1 = (Y - Z*(kroneckerProduct(U, U)*VS.transpose()))*VS;		
		derivF(Gamma,Y,Z,U,alpha1);			
		for(j=1;j<p;j++) max_tmp = MAX(max_tmp,Gamma.row(j).norm());
	}
	else max_tmp = Y.norm();	
	
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];
	max_tmp/=sqrt(n/r1);
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}
	return lambda;
}

void derivF_hdU(MatrixXd &Gamma, MatrixXd Y, MatrixXd Z, MatrixXd U, MatrixXd V, MatrixXd S, MatrixXd alpha, double lambda)
{
	int k, t, m, r1=opts.r1, p=opts.p;
	double Unorm;
	VectorXd temp,temp1,temp2, derivPen;
	Gamma=MatrixXd::Constant(p, r1, 0);
	derivPen = VectorXd::Constant(p, lambda);
	for(m = 0; m < opts.n; m++){
		for (k = 0; k < p; k++){
			for(t = 0; t < r1; t++){
				temp = temp1 = VectorXd::Constant(r1*r1, 0);
				temp2 = U.transpose()*subvector_extract(Z.row(m),k*p,p,1);
				subvector_updated(temp, temp2, t*r1,r1,1);
				subvector_updated(temp1, temp2, t,r1,r1);
				Gamma(k, t) -= innerProduct(alpha.row(m),temp+temp1,r1*r1);			
			}
		}
	}
	if(opts_pen.pen==2) for(k = 1; k < p; k++) derivPen[k] = MAX(lambda-U.row(k).norm()/opts_pen.gamma_pen,0);
	if(opts_pen.pen==3) for(k = 1; k < p; k++){ 
	    Unorm = U.row(k).norm();
		derivPen[k] = lambda*IDEX(lambda,Unorm)+(1-IDEX(lambda,Unorm))*MAX(opts_pen.gamma_pen*lambda-Unorm,0)/(opts_pen.gamma_pen-1);
	}	
	for(k = 1; k < p; k++) Gamma.row(k) += 2*opts.n*derivPen[k]*U.row(k)/MAX(U.row(k).norm(),opts_pen.thresh);
}

void derivF_hd(MatrixXd &Gamma, MatrixXd Y, MatrixXd Z, MatrixXd U, MatrixXd V, MatrixXd S, MatrixXd alpha, double lambda)
{
	int k, t, j, m, r1=opts.r1, r3=opts.r3, p=opts.p, n=opts.n;
	double Unorm0, gamma_pen=opts_pen.gamma_pen;
	MatrixXd dGU, UUS, djU, US;
	VectorXd temp,temp1,temp2, derivPen;
	VectorXi ind1=VectorXi::Constant(2*p-1, 0);
	Gamma=MatrixXd::Constant(p, r1, 0);
	derivPen = VectorXd::Constant(p, lambda);	
	for(m = 0; m < n; m++){
		for (k = 0; k < p; k++){
			for(t = 0; t < r1; t++){
				temp = temp1 = VectorXd::Constant(r1*r1, 0);
				temp2 = U.transpose()*subvector_extract(Z.row(m),k*p,p,1);
				subvector_updated(temp, temp2, t*r1,r1,1);
				subvector_updated(temp1, temp2, t,r1,r1);
				Gamma(k, t) -= innerProduct(alpha.row(m),temp+temp1,r1*r1);			
			}
		}
	}
	if(opts_pen.isPenColumn){
		VectorXd Unorm = VectorXd::Constant(p, 0);
		VectorXd USj;
		UUS = kroneckerProduct(U, U)*S.transpose();	
		dGU=MatrixXd::Constant(p, r1, 0);
		for (j = 1; j < p; j++){	 		
			for(m = 0; m < j; m++)	ind1[m] = m*p+j;
			for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
			for(m = j+1; m < p; m++)  ind1[m+p-1] = m*p+j;
			djU = extractCols(Z,2*p-1,ind1)*extractRows(UUS,2*p-1,ind1);		
			Unorm[j] = djU.norm()/sqrt(1.0*n);
			for (k = 0; k < p; k++){
				for(t = 0; t < r1; t++){								
					if(k==j){
						US = S.block(0, t*r1, r3, r1)*U.transpose();
						dGU(k,t) += 2*innerProductMatrix(extractCols_by(Z,k*p,p,1), djU*US);				
					}
					else{
						USj = S.block(0, t*r1, r3, r1)*U.row(j).transpose();
						dGU(k,t) += 2*innerProduct(Z.col(j*p+k),djU*USj,n);
					}
				}
			}
		}
		if(opts_pen.pen==2) for(k = 1; k < p; k++){
			Unorm0 = Unorm[k];		
			derivPen[k] = MAX(lambda-Unorm0/gamma_pen,0);
		}
		if(opts_pen.pen==3) for(k = 1; k < p; k++){ 
			Unorm0 = Unorm[k];
			derivPen[k] = lambda*IDEX(lambda,Unorm0)+IDEX(Unorm0,lambda)*MAX(gamma_pen*lambda-Unorm0,0)/(gamma_pen-1);
		}	
		for(k = 1; k < p; k++) Gamma.row(k) += 2*derivPen[k]*dGU.row(k)/MAX(Unorm[k],opts_pen.thresh);		
	}
	else for(int l=0;l<opts.q;l++){
		double USj;
		MatrixXd Unorm = MatrixXd::Constant(opts.q, p, 0);
		UUS = kroneckerProduct(U, U)*(V*S).transpose();	
		dGU=MatrixXd::Constant(p, r1, 0);
		for (j = 1; j < p; j++){	 		
			for(m = 0; m < j; m++)	ind1[m] = m*p+j;
			for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
			for(m = j+1; m < p; m++)  ind1[m+p-1] = m*p+j;
			djU = extractCols(Z,2*p-1,ind1)*extractRows(UUS.col(l),2*p-1,ind1);		
			Unorm(l,j) = djU.norm()/sqrt(1.0*n);
			for (k = 0; k < p; k++){
				for(t = 0; t < r1; t++){								
					if(k==j){
						US = V.row(l)*S.block(0, t*r1, r3, r1)*U.transpose();
						dGU(k,t) += 2*innerProduct(extractCols_by(Z,k*p,p,1)*US.transpose(), djU, n);				
					}
					else{
						USj = V.row(l)*S.block(0, t*r1, r3, r1)*U.row(j).transpose();
						dGU(k,t) += 2*innerProduct(Z.col(j*p+k),djU,n)*USj;
					}
				}
			}
		}	
		if(opts_pen.pen==2) for(k = 1; k < p; k++){
			Unorm0 = Unorm(l,k);		
			derivPen[k] = MAX(lambda-Unorm0/gamma_pen,0);
		}
		if(opts_pen.pen==3) for(k = 1; k < p; k++){ 
			Unorm0 = Unorm(l,k);
			derivPen[k] = lambda*IDEX(lambda,Unorm0)+IDEX(Unorm0,lambda)*MAX(gamma_pen*lambda-Unorm0,0)/(gamma_pen-1);
		}	
		for(k = 1; k < p; k++) Gamma.row(k) += 2*derivPen[k]*dGU.row(k)/MAX(Unorm(l,k),opts_pen.thresh);		
	}
}
//***-------------------------------------------------------------**
// ======================= updateU with penalty======================================
List updateU_Approx_onStief(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd U0, MatrixXd V, double lambda1)
{			
	int step=0, p=opts.p,isPenU=opts_pen.isPenU;
	double normPi, deriv, SY, Ck, tau=opts.tau_min;
	double  likhd0 = pow(10, 6), likhd1 = 0, Q=1, Qp;

	MatrixXd Gamma, VSUU, alpha, VS=V*S;
	MatrixXd U = U0, Old_U = U, dtX, dtXP,dtXdiff, Udiff;
	
	VSUU = kroneckerProduct(U, U)*VS.transpose();
	alpha = (Y - Z*VSUU)*VS;		
	if(isPenU)derivF_hdU(Gamma,Y,Z,U,V,S,alpha,lambda1);
	else derivF_hd(Gamma,Y,Z,U,V,S,alpha,lambda1);
	dtXP = Gamma - U*Gamma.transpose()*U;
	normPi=dtXP.norm();

	Ck = (Y - Z*VSUU).squaredNorm()/2;
	while (step<opts.max_step){
		step++;
		deriv=opts.rhols*normPi*normPi;	
		//tau = Curvilinear(U, VSUU, likhd1, tau, Ck, deriv, Y, Z, Gamma, VS, Old_U, 1, lambda1);		
		tau = Curvilinear(U, VSUU, likhd1, tau, Ck, deriv, Y, Z, Gamma, VS, Old_U, 0, lambda1);		
		alpha = (Y - Z*VSUU)*VS;
		if(isPenU) derivF_hdU(Gamma,Y,Z,U,V,S,alpha,lambda1);
		else derivF_hd(Gamma,Y,Z,U,V,S,alpha,lambda1);
		dtX = Gamma - U*Gamma.transpose()*U;
		normPi=dtX.norm();			
		dtXdiff = dtX - dtXP;     
		Udiff = U-Old_U;
		SY = fabs(innerProductMatrix(Udiff,dtXdiff));
		if(step%2) tau  = SY/dtXdiff.squaredNorm();
		else tau = Udiff.squaredNorm()/SY; 
		tau = MAX(MIN(tau, 1e20), 1e-20);
		if(fabs(likhd1-likhd0)/fabs(likhd0+1)<opts.ftol && Udiff.norm()/sqrt(p) < opts.utol) break;
		if(normPi<opts.Pitol) break;
		if (likhd1>likhd0) break;
		Old_U = U;
		likhd0 = likhd1;
		dtXP = dtX;
        Qp = Q; 
	    Q = opts.gamma*Qp + 1; 
	    Ck = (opts.gamma*Qp*Ck + likhd1)/Q;		
	}
	likhd1 = lossFunction(Y,Z,VS,U,VSUU,0,0,1);
	return List::create(Named("likhd") = likhd1, Named("U") = U, Named("Unorm") = (U-Old_U).norm(), Named("step")=step);
}

//***-------------------------------------------------------------**
// ======================= updateU with penalty======================================
List updateU_Approx_onR(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd U0, MatrixXd V, double lambda1)
{			
	int step=0, p=opts.p,isPenU=opts_pen.isPenU;
	double normPi, deriv, Ck, tau=opts.tau_min;
	double  likhd0 = pow(10, 6), likhd1 = 0, Q=1, Qp;

	MatrixXd Gamma, VSUU, alpha, VS=V*S;
	MatrixXd U = U0, Old_U = U, Udiff;
	
	VSUU = kroneckerProduct(U, U)*VS.transpose();
	alpha = (Y - Z*VSUU)*VS;		
	if(isPenU)derivF_hdU(Gamma,Y,Z,U,V,S,alpha,lambda1);
	else derivF_hd(Gamma,Y,Z,U,V,S,alpha,lambda1);
	normPi=Gamma.norm();

	Ck = (Y - Z*VSUU).squaredNorm()/2;
	while (step<opts.max_step){
		step++;
		deriv=opts.rhols*normPi*normPi;			
		tau = Curvilinear(U, VSUU, likhd1, tau, Ck, deriv, Y, Z, Gamma, VS, Old_U, 0, lambda1);		
		alpha = (Y - Z*VSUU)*VS;
		if(isPenU) derivF_hdU(Gamma,Y,Z,U,V,S,alpha,lambda1);
		else derivF_hd(Gamma,Y,Z,U,V,S,alpha,lambda1);		
		normPi=Gamma.norm();
		
		Udiff = U-Old_U;	
		tau = MAX(MIN(tau, 1e20), 1e-20);
		if(fabs(likhd1-likhd0)/fabs(likhd0+1)<opts.ftol && Udiff.norm()/sqrt(p) < opts.utol) break;
		if (likhd1>likhd0) break;
		Old_U = U;
		likhd0 = likhd1;
        Qp = Q; 
	    Q = opts.gamma*Qp + 1; 
	    Ck = (opts.gamma*Qp*Ck + likhd1)/Q;		
	}
	likhd1 = lossFunction(Y,Z,VS,U,VSUU,0,0,1);
	MatrixXd VR,VQ,R,Snew;
	HouseholderQR<MatrixXd> qr;
	qr.compute(U);
	VR = qr.matrixQR().triangularView<Upper>();
	R = VR.block(0, 0, opts.r1, opts.r1);
	VQ = qr.householderQ();
	U = VQ.block(0, 0, p, opts.r1);
	Snew = S*kroneckerProduct(R, R).transpose();	
	return List::create(Named("likhd") = likhd1, Named("U") = U, Named("Snew") = Snew, Named("step")=step);
}

//***-------------------------------------------------------------**
//***------main function: Estimation with penalizing coefficients in a whole column -----------------**
// [[Rcpp::export]]
List EstPenColumn(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd U, MatrixXd V, VectorXd lambda, List optsList, List optsList_pen)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	U is p*r1 matrix
	V is q*r3 matrix
	S is r3*(r1*r2) matrix

	Output:
	Dnew is a estimator of D(3) = V%*%S%*%kronecker(t(U),t(U))
	*/
	//Rcout << "U= " << U << std::endl;
	
	opts.utol = as<double>(optsList[0]);
	opts.ftol = as<double>(optsList[1]);
	opts.Pitol = as<double>(optsList[2]);
	opts.tau_min = as<double>(optsList[3]);
	opts.eta = as<double>(optsList[4]);
	opts.tiny = as<double>(optsList[5]);
	opts.gamma = as<double>(optsList[6]);
	opts.rhols = as<double>(optsList[7]);
	opts.max_step = as<int>(optsList[8]);
	opts.max_step1 = as<int>(optsList[9]);
	opts.is_LR = as<int>(optsList[10]);
	opts.n = as<int>(optsList[11]);
	opts.r1 = as<int>(optsList[12]);
	opts.r2 = as<int>(optsList[13]);
	opts.r3 = as<int>(optsList[14]);
	opts.p = as<int>(optsList[15]);
	opts.q = as<int>(optsList[16]);
	opts.onStiefel = as<int>(optsList[17]);
	
	opts_pen.pen = as<int>(optsList_pen[0]);
	opts_pen.nlam = as<int>(optsList_pen[1]);
	opts_pen.lam_max = as<double>(optsList_pen[2]);
	opts_pen.lam_min = as<double>(optsList_pen[3]);
	opts_pen.gamma_pen = as<double>(optsList_pen[4]);
	opts_pen.alpha = as<double>(optsList_pen[5]);
	opts_pen.dfmax = as<int>(optsList_pen[6]);
	opts_pen.gamma_tanh = as<double>(optsList_pen[7]);
	opts_pen.thresh = as<double>(optsList_pen[8]);	
	opts_pen.isPenU = as<int>(optsList_pen[9]);
	opts_pen.isPenColumn = as<int>(optsList_pen[10]);	
	opts_pen.delta = as<double>(optsList_pen[11]);
	opts_pen.max_step = as<double>(optsList_pen[12]);
	opts_pen.max_step1 = as<double>(optsList_pen[13]);
	opts_pen.isFISC = as<double>(optsList_pen[14]);		
	
	int l,j, step,nlam=opts_pen.nlam, p=opts.p, q=opts.q, r1=opts.r1, r3=opts.r3;
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0, thresh=0;
	MatrixXd Dnew, Unew, Snew, Vnew, Z = produceZ(X);
	VectorXd activeA, convergence1;
	List fit;
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	VectorXd df = VectorXd::Constant(nlam, 0);
	MatrixXd betapath, Upath, Vpath, Spath, temp;
	betapath = MatrixXd::Constant(p-1, nlam, 0);
	Upath = MatrixXd::Constant(p*r1, nlam, 0);
	Vpath = MatrixXd::Constant(q*r3, nlam, 0);
	Spath = MatrixXd::Constant(r1*r1*r3, nlam, 0);	

	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<opts.max_step) {
		  convergence1 = VectorXd::Constant(3, 1);
			step ++;
			
			Snew = updateS(Y, Z, U, V);
			Dnew = V*Snew*kroneckerProduct(U.transpose(), U.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;	
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;		
			fit = updateV(Y, Z, U, S);
			Vnew = fit[0];
			Snew = fit[1];
			Dnew = Vnew*Snew*kroneckerProduct(U.transpose(), U.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;			
			if (likhd1<likhd0) {
				V = Vnew;
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;			
			if(opts.onStiefel){
				fit = updateU_Approx_onStief(Y, Z, S, U, V, lambda1);	
				Unew = fit[1];		
				Dnew = V*S*kroneckerProduct(Unew.transpose(), Unew.transpose());
			}
			else{
				fit = updateU_Approx_onR(Y, Z, S, U, V, lambda1);
				Unew = fit[1];
				Snew = fit[2];
				Dnew = V*Snew*kroneckerProduct(Unew.transpose(), Unew.transpose());
			}
			likhd1 = fit[0];
			if (likhd1<likhd0) {
				U = Unew;
				if(!opts.onStiefel) S = Snew;
				if ((likhd0 - likhd1) / likhd0<opts.ftol) break;
				else  likhd0 = likhd1;		
			}
			else convergence1[2]=0;
			if(convergence1.sum()==0) break;
		} //end while		
		activeA = VectorXd::Constant(p-1, 0);
		if(opts_pen.isPenU){
			thresh = sqrt(r1)*opts_pen.thresh;
			for(j=1;j<p;j++) if(U.row(j).norm()>thresh) activeA[j-1] = 1;
		}
		else{
			int m;
			MatrixXd UUS, djU;
			VectorXi ind1=VectorXi::Constant(2*p-1, 0);
			UUS = kroneckerProduct(U, U)*S.transpose();
			thresh = sqrt(opts.n*r1)*2*p*opts_pen.thresh;
			for (j = 1; j < p; j++){	        
				for(m = 0; m < j; m++)	ind1[m] = m*p+j;
				for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
				for(m = j+1; m < p; m++)	ind1[m+p-1] = m*p+j;
				djU = extractCols(Z,2*p-1,ind1)*extractRows(UUS,2*p-1,ind1);
				if(djU.norm()>thresh) activeA[j-1] = 1;	
			}	
		}	
		df[l] = activeA.sum();
		likhd[l] = likhd0;		
		betapath.col(l) = activeA;
		temp = U; temp.resize(p*r1, 1);
		Upath.col(l) = temp;
		temp = V; temp.resize(q*r3, 1);
		Vpath.col(l) = temp;
		temp = S; temp.resize(r1*r1*r3, 1);
		Spath.col(l) = temp;
	}// end for
	return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda,Named("Upath")=Upath,Named("Vpath")=Vpath,Named("Spath")=Spath);
}

//***-------------------------------------------------------------**
//***------main function: Estimation with penalizing coefficient in single row -----------------**
// [[Rcpp::export]]
List EstPenSingle(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd U, MatrixXd V, VectorXd lambda, List optsList, List optsList_pen)
{
	/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	U is p*r1 matrix
	V is q*r3 matrix
	S is r3*(r1*r2) matrix

	Output:
	Dnew is a estimator of D(3) = V%*%S%*%kronecker(t(U),t(U))
	*/
	opts.utol = as<double>(optsList[0]);
	opts.ftol = as<double>(optsList[1]);
	opts.Pitol = as<double>(optsList[2]);
	opts.tau_min = as<double>(optsList[3]);
	opts.eta = as<double>(optsList[4]);
	opts.tiny = as<double>(optsList[5]);
	opts.gamma = as<double>(optsList[6]);
	opts.rhols = as<double>(optsList[7]);
	opts.max_step = as<int>(optsList[8]);
	opts.max_step1 = as<int>(optsList[9]);
	opts.is_LR = as<int>(optsList[10]);
	opts.n = as<int>(optsList[11]);
	opts.r1 = as<int>(optsList[12]);
	opts.r2 = as<int>(optsList[13]);
	opts.r3 = as<int>(optsList[14]);
	opts.p = as<int>(optsList[15]);
	opts.q = as<int>(optsList[16]);
	opts.onStiefel = as<int>(optsList[17]);
	
	opts_pen.pen = as<int>(optsList_pen[0]);
	opts_pen.nlam = as<int>(optsList_pen[1]);
	opts_pen.lam_max = as<double>(optsList_pen[2]);
	opts_pen.lam_min = as<double>(optsList_pen[3]);
	opts_pen.gamma_pen = as<double>(optsList_pen[4]);
	opts_pen.alpha = as<double>(optsList_pen[5]);
	opts_pen.dfmax = as<int>(optsList_pen[6]);
	opts_pen.gamma_tanh = as<double>(optsList_pen[7]);
	opts_pen.thresh = as<double>(optsList_pen[8]);	
	opts_pen.isPenU = as<int>(optsList_pen[9]);
	opts_pen.isPenColumn = as<int>(optsList_pen[10]);	
	opts_pen.delta = as<double>(optsList_pen[11]);
	opts_pen.max_step = as<double>(optsList_pen[12]);
	opts_pen.max_step1 = as<double>(optsList_pen[13]);
	opts_pen.isFISC = as<double>(optsList_pen[14]);	
	
	int l,j, step,nlam=opts_pen.nlam, p=opts.p, q=opts.q, r1=opts.r1, r3=opts.r3;
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0, thresh=0;
	MatrixXd Dnew, Unew, Snew, Vnew, Z = produceZ(X);
	VectorXd convergence1;
	List fit;
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	VectorXi df = VectorXi::Constant(nlam, 0), activeX;
	MatrixXd Upath, Vpath, Spath, temp;
	MatrixXi betapath = MatrixXi::Constant(q*(p-1), nlam, 0), activeXpath = MatrixXi::Constant(p-1, nlam, 0);;
	Upath = MatrixXd::Constant(p*r1, nlam, 0);
	Vpath = MatrixXd::Constant(q*r3, nlam, 0);
	Spath = MatrixXd::Constant(r1*r1*r3, nlam, 0);	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<opts.max_step) {
		  convergence1 = VectorXd::Constant(3, 1);
			step ++;
			
			Snew = updateS(Y, Z, U, V);
			Dnew = V*Snew*kroneckerProduct(U.transpose(), U.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;
		
			if (likhd1<likhd0) {
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[0]=0;		
			fit = updateV(Y, Z, U, S);
			Vnew = fit[0];
			Snew = fit[1];
			Dnew = Vnew*Snew*kroneckerProduct(U.transpose(), U.transpose());
			likhd1 = (Y - Z * Dnew.transpose()).squaredNorm()/2;			

			if (likhd1<likhd0) {
				V = Vnew;
				S = Snew;
				likhd0 = likhd1;
			}
			else convergence1[1]=0;			
			if(opts.onStiefel){
				fit = updateU_Approx_onStief(Y, Z, S, U, V, lambda1);	
				Unew = fit[1];		
				Dnew = V*S*kroneckerProduct(Unew.transpose(), Unew.transpose());
			}
			else{
				fit = updateU_Approx_onR(Y, Z, S, U, V, lambda1);
				Unew = fit[1];
				Snew = fit[2];
				Dnew = V*Snew*kroneckerProduct(Unew.transpose(), Unew.transpose());
			}
			likhd1 = fit[0];
			if (likhd1<likhd0) {
				U = Unew;
				if(!opts.onStiefel) S = Snew;
				if ((likhd0 - likhd1) / likhd0<opts.ftol) break;
				else  likhd0 = likhd1;		
			}
			else convergence1[2]=0;
			if(convergence1.sum()==0) break;
		} //end while		
		activeX = VectorXi::Constant(p-1, 0);
		if(opts_pen.isPenU){	
		    thresh = sqrt(r1)*opts_pen.thresh;
			for(j=1;j<p;j++) if(U.row(j).norm()>thresh) activeX[j-1] = 1;
		}
		else{
			int m;		
			VectorXi ind1=VectorXi::Constant(2*p-1, 0);
			MatrixXd UUS, djU;
			
			UUS = kroneckerProduct(U, U)*(V*S).transpose();
			MatrixXi activeA = MatrixXi::Constant(q, p-1, 0);
			thresh = sqrt(opts.n*r1)*2*p*opts_pen.thresh;
			for(int k=0;k<q;k++){
				for (j = 1; j < p; j++){	        
					for(m = 0; m < j; m++)	ind1[m] = m*p+j;
					for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
					for(m = j+1; m < p; m++)	ind1[m+p-1] = m*p+j;
					djU = extractCols(Z,2*p-1,ind1)*extractRows(UUS.col(k),2*p-1,ind1);
					if(djU.norm()>thresh/q) activeA(k,j-1) = 1;	
				}
			}
			for (j = 1; j < p; j++){	        
				for(m = 0; m < j; m++)	ind1[m] = m*p+j;
				for(m = 0; m < p; m++)	ind1[m+j] = j*p+m;
				for(m = j+1; m < p; m++)	ind1[m+p-1] = m*p+j;
				djU = extractCols(Z,2*p-1,ind1)*extractRows(UUS,2*p-1,ind1);
				if(djU.norm()>thresh) activeX[j-1] = 1;	
			}	
			activeA.resize(q*(p-1),1);
			betapath.col(l) = activeA;			
		}			
		df[l] = activeX.sum();
		activeXpath.col(l) = activeX;
		likhd[l] = likhd0;		
		temp = U; temp.resize(p*r1, 1);
		Upath.col(l) = temp;
		temp = V; temp.resize(q*r3, 1);
		Vpath.col(l) = temp;
		temp = S; temp.resize(r1*r1*r3, 1);
		Spath.col(l) = temp;
	}// end for
	return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda,Named("Upath")=Upath,Named("Vpath")=Vpath,Named("Spath")=Spath, Named("activeXpath") = activeXpath);
}
//----------------------------------------------------------------**
//***--------------------EstimationD3 directly--------------------**
// [[Rcpp::export]]
List EstimationD3(MatrixXd Y, MatrixXd X)
{
  int j,q = Y.cols(),p = X.cols(), kp=p*p;
  MatrixXd Dnew = MatrixXd::Constant(q, kp, 0), Z = produceZ(X);
  HouseholderQR<MatrixXd> qr;
  qr.compute(Z);
  MatrixXd R = qr.matrixQR().triangularView<Upper>();
  MatrixXd Q = qr.householderQ();
  MatrixXd RQ = UpTriangularInv(R) * Q.transpose();
  
  MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
  if (pow(condition_numberQRSym(R),2) > 1e10){
    temp = tRbyR(R) + (IDEN.array()*1e-4).matrix();
    for (j = 0; j < q; j++) Dnew.row(j) = (temp.colPivHouseholderQr().solve(Z.transpose()*Y.col(j))).transpose();
  }
  else
    for (j = 0; j < q; j++) Dnew.row(j) = (RQ * Y.col(j)).transpose();
  double likhd = (Y - Z * Dnew.transpose()).squaredNorm()/2;
  return List::create(Named("likhd") = likhd, Named("Dnew") = Dnew);
}





























