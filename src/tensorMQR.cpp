//[[Rcpp::depends(RcppEigen)]]
#include "mqr.h"

//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
// [[Rcpp::export]]
List EstUnconstr(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, VectorXd mu, List optsList)
{	
	opts.n = as<int>(optsList["n"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.q = as<int>(optsList["q"]);
	opts.p = as<int>(optsList["p"]);
	opts.K = as<int>(optsList["p"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.intercept = as<int>(optsList["intercept"]);
	
	double  likhd0 = pow(10, 6), likhd1 = 0;
	int n = Y.rows(),q = opts.q;
	opts.n=n;
	MatrixXd Dnew,Cnew,Anew,Bnew,Snew, Y1=Y, Z = produceZ(X);
	VectorXi convergence1;
	VectorXd Ones;
	Ones.setOnes(n);
	List fit;
	int step = 0;
	while (step<opts.max_step1) {
		convergence1 = VectorXi::Constant(4, 1);
		step = step + 1;
		Snew = updateS(Y1, Z, A, B, C);
		Dnew = C * Snew*kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;
		fit = updateC(Y1, Z, A, B, C, S);
		Cnew = fit[0];
		Snew = fit[1];
		Dnew = Cnew * Snew *kroneckerProduct(B.transpose(), A.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			C = Cnew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;
		fit = updateA(Y1, Z, A, B, C, S);
		Anew = fit[0];
		Snew = fit[1];
		Dnew = C * Snew *kroneckerProduct(B.transpose(), Anew.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			A = Anew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[2]=0;
		fit = updateB(Y1, Z, A, B, C, S);
		Bnew = fit[0];
		Snew = fit[1];		
		Dnew = C * Snew *kroneckerProduct(Bnew.transpose(), A.transpose());
		likhd1 = (Y1  - Z * Dnew.transpose()).squaredNorm();
		if (likhd1<likhd0) {
			B = Bnew;
			S = Snew;
			if(opts.intercept){
				mu = (Y - Z * Dnew.transpose()).colwise().sum()/n;
				Y1 = Y - kroneckerProduct(Ones,mu);
				Y1.resize(n,q);	
			}
			if ((likhd0 - likhd1) / likhd0<opts.eps) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[3]=0;
		if(convergence1.sum()==0) break;
	}
	Dnew = C * S * kroneckerProduct(B.transpose(), A.transpose());
	return List::create(Named("likhd") = likhd0, Named("Dnew") = Dnew, Named("A") = A, Named("B") = B, Named("C") = C,Named("S") = S, Named("mu") = mu);
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
MatrixXd updateS_sym(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd C)
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
List Estimation(MatrixXd Y, MatrixXd X, MatrixXd S, MatrixXd U, MatrixXd V, VectorXd mu, List optsList)
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
	opts.utol = as<double>(optsList["utol"]);
	opts.ftol = as<double>(optsList["ftol"]);
	opts.Pitol = as<double>(optsList["Pitol"]);
	opts.tau_min = as<double>(optsList["tau_min"]);
	opts.eta = as<double>(optsList["eta"]);
	opts.tiny = as<double>(optsList["tiny"]);
	opts.gamma = as<double>(optsList["gamma"]);
	opts.rhols = as<double>(optsList["rhols"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.is_LR = as<int>(optsList["isLR"]);
	opts.n = as<int>(optsList["n"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.intercept = as<int>(optsList["intercept"]);
	
	int n = Y.rows(), q = opts.q;
	double  likhd0 = pow(10, 6), likhd1 = 0;
	opts.n = n;
	MatrixXd Dnew, Unew, Vnew, Snew, Y1=Y, Z = produceZ(X);
	VectorXd convergence1;	
	VectorXd Ones;
	Ones.setOnes(n);
	List fit;
	int step = 0;
	while (step<opts.max_step1) {
		convergence1 = VectorXd::Constant(3, 1);
		step++;
		
		Snew = updateS_sym(Y1, Z, U, V);
		Dnew = V*Snew*kroneckerProduct(U.transpose(), U.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm()/2;	
		if (likhd1<likhd0) {
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[0]=0;		
		fit = updateV(Y1, Z, U, S);
		Vnew = fit[0];
		Snew = fit[1];
		Dnew = Vnew*Snew*kroneckerProduct(U.transpose(), U.transpose());
		likhd1 = (Y1 - Z * Dnew.transpose()).squaredNorm()/2;
		if (likhd1<likhd0) {
			V = Vnew;
			S = Snew;
			likhd0 = likhd1;
		}
		else convergence1[1]=0;
		fit = updateU(Y1, Z, V*S, U);	
        Unew = fit[1];		
		Dnew = V*S*kroneckerProduct(Unew.transpose(), Unew.transpose());
		likhd1 = fit[0];
		if (likhd1<likhd0) {
			U = Unew;
			if(opts.intercept){
				mu = (Y - Z * Dnew.transpose()).colwise().sum()/n;
				Y1 = Y - kroneckerProduct(Ones,mu);
				Y1.resize(n,q);	
			}
			if ((likhd0 - likhd1) / likhd0<opts.ftol) break;
			else  likhd0 = likhd1;		
		}
		else convergence1[2]=0;
		if(convergence1.sum()==0) break;
	}
	Dnew = V*S*kroneckerProduct(U.transpose(), U.transpose());
	return List::create(Named("likhd") = 2*likhd0, Named("Dnew") = Dnew, Named("S") = S, Named("U")=U, Named("V")=V, Named("mu") = mu);
}

void derivF_hdU(MatrixXd &Gamma, MatrixXd Y, MatrixXd Z, MatrixXd U, MatrixXd V, MatrixXd S, MatrixXd alpha, double lambda)
{
	int k, t, m, r1=opts.r1, p=opts.p, n=opts.n;
	double Unorm;
	VectorXd temp,temp1,temp2, derivPen;
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
	if(opts_pen.pen==2) for(k = 1; k < p; k++) derivPen[k] = MAX(lambda-U.row(k).norm()/opts_pen.gamma_pen,0);
	if(opts_pen.pen==3) for(k = 1; k < p; k++){ 
	    Unorm = U.row(k).norm();
		derivPen[k] = lambda*IDEX(lambda,Unorm)+(1-IDEX(lambda,Unorm))*MAX(opts_pen.gamma_pen*lambda-Unorm,0)/(opts_pen.gamma_pen-1);
	}	
	for(k = 1; k < p; k++) Gamma.row(k) += 2*n*derivPen[k]*U.row(k)/MAX(U.row(k).norm(),opts_pen.thresh);
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
	
	opts.utol = as<double>(optsList["utol"]);
	opts.ftol = as<double>(optsList["ftol"]);
	opts.Pitol = as<double>(optsList["Pitol"]);
	opts.tau_min = as<double>(optsList["tau_min"]);
	opts.eta = as<double>(optsList["eta"]);
	opts.tiny = as<double>(optsList["tiny"]);
	opts.gamma = as<double>(optsList["gamma"]);
	opts.rhols = as<double>(optsList["rhols"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.is_LR = as<int>(optsList["isLR"]);
	opts.n = as<int>(optsList["n"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.intercept = as<int>(optsList["intercept"]);
	opts.onStiefel = as<int>(optsList["onStiefel"]);	

	opts_pen.gamma_tanh = as<double>(optsList_pen["gamma_tanh"]);
	opts_pen.thresh = as<double>(optsList_pen["thresh"]);	
	opts_pen.isPenU = as<int>(optsList_pen["isPenU"]);
	opts_pen.delta = as<double>(optsList_pen["delta"]);
	opts_pen.isFISC = as<double>(optsList_pen["isFISC"]);
	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();
	
	int l,j, step,nlam=opts_pen.nlam, p=opts.p, q=opts.q, r1=opts.r1, r3=opts.r3,n=Y.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0, thresh=0;
	opts.n=n;
	MatrixXd Dnew, Unew, Snew, Vnew, Z = produceZ(X);
	VectorXi activeA, convergence1, df;	
	df = VectorXi::Constant(nlam, 0);
	MatrixXi betapath = MatrixXi::Constant(p-1, nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);	
	MatrixXd Upath, Vpath, Spath, temp;	
	Upath = MatrixXd::Constant(p*r1, nlam, 0);
	Vpath = MatrixXd::Constant(q*r3, nlam, 0);
	Spath = MatrixXd::Constant(r1*r1*r3, nlam, 0);	
	List fit;

	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<opts.max_step1) {
		  convergence1 = VectorXi::Constant(3, 1);
			step ++;
			
			Snew = updateS_sym(Y, Z, U, V);
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
		activeA = VectorXi::Constant(p-1, 0);
		if(opts_pen.isPenU){
			thresh = sqrt(r1)*opts_pen.thresh;
			for(j=1;j<p;j++) if(U.row(j).norm()>thresh) activeA[j-1] = 1;
		}
		else{
			int m;
			MatrixXd UUS, djU;
			VectorXi ind1=VectorXi::Constant(2*p-1, 0);
			UUS = kroneckerProduct(U, U)*S.transpose();
			thresh = sqrt(n*r1)*2*p*opts_pen.thresh;
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
	return List::create(Named("likhd") = 2*likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda,Named("Upath")=Upath,Named("Vpath")=Vpath,Named("Spath")=Spath);
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
	opts.utol = as<double>(optsList["utol"]);
	opts.ftol = as<double>(optsList["ftol"]);
	opts.Pitol = as<double>(optsList["Pitol"]);
	opts.tau_min = as<double>(optsList["tau_min"]);
	opts.eta = as<double>(optsList["eta"]);
	opts.tiny = as<double>(optsList["tiny"]);
	opts.gamma = as<double>(optsList["gamma"]);
	opts.rhols = as<double>(optsList["rhols"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);
	opts.max_step1 = as<int>(optsList["max_step1"]);
	opts.is_LR = as<int>(optsList["isLR"]);
	opts.n = as<int>(optsList["n"]);
	opts.r1 = as<int>(optsList["r1"]);
	opts.r2 = as<int>(optsList["r2"]);
	opts.r3 = as<int>(optsList["r3"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.intercept = as<int>(optsList["intercept"]);
	opts.onStiefel = as<int>(optsList["onStiefel"]);	

	opts_pen.gamma_tanh = as<double>(optsList_pen["gamma_tanh"]);
	opts_pen.thresh = as<double>(optsList_pen["thresh"]);	
	opts_pen.isPenU = as<int>(optsList_pen["isPenU"]);
	opts_pen.delta = as<double>(optsList_pen["delta"]);
	opts_pen.isFISC = as<double>(optsList_pen["isFISC"]);
	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();	
	
	int l,j, step,nlam=opts_pen.nlam, p=opts.p, q=opts.q, r1=opts.r1, r3=opts.r3,n=Y.rows();
	double  likhd0 = pow(10, 6), lambda1, likhd1 = 0, thresh=0;
	opts.n=n;
	MatrixXd Dnew, Unew, Snew, Vnew, Z = produceZ(X);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	VectorXi df = VectorXi::Constant(nlam, 0), activeX, convergence1;
	MatrixXd Upath, Vpath, Spath, temp;
	MatrixXi betapath = MatrixXi::Constant(q*(p-1), nlam, 0), activeXpath = MatrixXi::Constant(p-1, nlam, 0);;
	Upath = MatrixXd::Constant(p*r1, nlam, 0);
	Vpath = MatrixXd::Constant(q*r3, nlam, 0);
	Spath = MatrixXd::Constant(r1*r1*r3, nlam, 0);	
	List fit;
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<opts.max_step1) {
		  convergence1 = VectorXi::Constant(3, 1);
			step ++;
			
			Snew = updateS_sym(Y, Z, U, V);
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
	return List::create(Named("likhd") = 2*likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda,Named("Upath")=Upath,Named("Vpath")=Vpath,Named("Spath")=Spath, Named("activeXpath") = activeXpath);
}