//[[Rcpp::depends(RcppEigen)]]
#include "mvrhd.h"

//----------------------------------------------------------------**
//***----update the jth row of matrix A with penalty--------------**
VectorXd update_colj(VectorXd z, double lambda, double alpha, double gamma, int penalty)
{
	double znorm = z.norm();
	znorm = penalties(znorm, 1, lambda, alpha, gamma, penalty) / znorm;
	return znorm * z;
}
//----------------------------------------------------------------**
//***----update the jth row of matrix A with penalty--------------**
MatrixXd update_blockj(MatrixXd z, double lambda, double alpha, double gamma, int penalty)
{
	double znorm = z.norm();
	znorm = penalties(znorm, 1, lambda, alpha, gamma, penalty) / znorm;
	return znorm * z;
}
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
VectorXd setuplambdaMVR_colwise(MatrixXd Y, MatrixXd Z, int nlam, VectorXd setlam)
{
	int n=Y.rows(), p = Z.cols(), q = Y.cols(), j;
	double lam_max, lam_min, alpha, max_lam, max_tmp=0;
	VectorXd lambda, lambda1, znorm;	
	znorm = Z.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) max_tmp = MAX(max_tmp,(Y.transpose()*Z.col(j)).norm()/znorm[j]);	
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];
	max_tmp/=n*sqrt(q);
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
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
MatrixXd setuplambdaMVR_lasso(MatrixXd Y, MatrixXd Z, int nlam, VectorXd setlam)
{
	int n = Y.rows(), p = Z.cols(), q = Y.cols(), j, k;
	double lam_max, lam_min, alpha, max_lam, max_tmp;
	VectorXd lambda, lambda1, znorm;	
	MatrixXd lambda_all = MatrixXd::Constant(nlam, q, 0);
	
	znorm = Z.colwise().norm()/sqrt(n);	
	for(k=0; k<q; k++){	
    	max_tmp=0;
		for(j=0; j<p; j++) max_tmp = MAX(max_tmp,fabs(Y.col(k).transpose()*Z.col(j))/znorm[j]);	
		lam_max = setlam[0];
		lam_min = setlam[1];
		alpha = setlam[2];
		max_tmp/=n;
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
		lambda_all.col(k) = lambda;
	}
	return lambda_all;
}
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
VectorXd setuplambdaMVR_blockwise(MatrixXd Y, MatrixXd Z, int nlam, VectorXd setlam, VectorXi lengths)
{
	int n = Y.rows(), p = lengths.size(), q = Y.cols(), j, d = 0, dj;
	double lam_max, lam_min, alpha, max_lam, max_tmp=0;
	VectorXd lambda, lambda1;	
	MatrixXd V, L, Vnorm;
	for(j=0;j<p;j++){		 
	    dj = lengths[j];
		V = Z.middleCols(d,dj);
		L = (V.transpose()*V/n).llt().matrixL();
		Vnorm = QbyR(V, UpTriangularInv(L.transpose()),1);
		max_tmp = MAX(max_tmp,(Y.transpose()*Vnorm).norm()/(n*sqrt(q*dj)));
        d += dj;		
	}
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];
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
//----------------------------------------------------------------**
//***-------------setup tuning parameters for MVR-----------------**
// [[Rcpp::export]]
MatrixXd setuplambdaMVR_glasso(MatrixXd Y, MatrixXd Z, int nlam, VectorXd setlam, VectorXi lengths)
{
	int n=Y.rows(), p = lengths.size(), q = Y.cols(), j, k, d = 0, dj;
	double lam_max, lam_min, alpha, max_lam, max_tmp;
	VectorXd lambda, lambda1, znorm;	
	MatrixXd lambda_all = MatrixXd::Constant(nlam, q, 0);
	MatrixXd V, L, Vnorm, zynorm;
	zynorm.setZero(p, q);

	for(j=0;j<p;j++){		 
	    dj = lengths[j];
		V = Z.middleCols(d,dj);
		L = (V.transpose()*V/n).llt().matrixL();
		Vnorm = QbyR(V, UpTriangularInv(L.transpose()),1);
		zynorm.row(j) = (Vnorm.transpose()*Y).colwise().norm()/(n*sqrt(dj));
        d += dj;		
	}	
	for(k=0; k<q; k++){	
    	max_tmp=0;
		for(j=0; j<p; j++)	max_tmp = MAX(max_tmp,zynorm(j,k));	       
		lam_max = setlam[0];
		lam_min = setlam[1];
		alpha = setlam[2];
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
		lambda_all.col(k) = lambda;
	}
	return lambda_all;
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
List MVR_colwise(MatrixXd Y, MatrixXd Z1, MatrixXd W, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd)
{
/*
	Input:
	Y is n*q matrix
	Z is n*p matrix, which is standardized
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is q*p coefficient metrix
*/  
	int l,j, active, step, nlam = lambda.size(), n = Y.rows(), q = Y.cols(), p = Z1.cols();
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, gamma = opts_pen.gamma_pen, penalty = opts_pen.pen, pz=opts.pz;
	double alpha = opts_pen.alpha, eps = opts.eps;
	double lambda1,diffmax=0,diffnorm;

    MatrixXd Anew=MatrixXd::Constant(q, p, 0), Bnew, Z = Z1;
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0), Cbeta, Cnew0, Cnew;
	VectorXd ajnew,gj, ajnorm, znorm, diff;
	if(pz){
		Cbeta.setZero(pz*q,nlam);
		Cnew0.setZero(pz,q);
	}	
	else{
		Cnew0.setZero(1,q);
		Cbeta.setZero(q,nlam);
	}
	MatrixXd r = Y;
	ajnorm.setZero(p);
	znorm = Z1.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) Z.col(j) = Z1.col(j)/znorm[j];	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<max_step) {
			step++;
			active = 0;
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active++;
			if (active>dfmax)
				return List::create(Named("betapath") = NULL,Named("Cbeta") = NULL);
			if(pz){
				Cnew = W.transpose()*r/n;
				r -= W*Cnew; 
				Cnew += Cnew0;
				Cnew0 = Cnew;
			}		
			for (j = 0; j < p; j++) {
				gj = r.transpose()*Z.col(j)/n + Anew.col(j);			
				ajnew = update_colj(gj, lambda1, alpha, gamma, penalty);
				diff = ajnew-Anew.col(j);
				diffnorm = diff.norm();
				if(diffnorm!=0){
					r -= kroneckerProduct(Z.col(j),diff.transpose());
					Anew.col(j) = ajnew;
					ajnorm[j] = ajnew.norm();
					if(diffnorm>diffmax) diffmax = diffnorm;
				}
			}
			if(diffmax<eps) break;
		}//end while
		for (j = 0; j<p; j++) 	if (ajnorm[j]) activeA(j,l) = 1;
		if(pz){	
			likhd[l] = (Y - W*Cnew- Z * Anew.transpose()).squaredNorm();
			Cnew.resize(pz*q,1);
			Cbeta.col(l) = Cnew;				
		}	
		else  likhd[l] = (Y - Z * Anew.transpose()).squaredNorm();					
		Bnew = Anew.transpose();
		for(j=0; j<p; j++) Bnew.row(j) /= znorm[j];
		Bnew.resize(p*q,1);
		betapath.col(l) = Bnew;
	}//end for	
	return List::create(Named("betapath") = betapath, Named("Cbeta") = Cbeta);
}
//***-------------------------------------------------------------** 
//***-------update the jth row of matrix A with penalty-----------** 
List MVR_lasso(VectorXd Y, MatrixXd Z, MatrixXd W, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd)
{
/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is p-dimensional coefficient
*/  
	int l,j, active, step, n = opts.n, p = opts.p, nlam = lambda.size();
    double ajnew,gj, lambda1,diffmax, diff;
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, gamma = opts_pen.gamma_pen, penalty = opts_pen.pen, pz=opts.pz;
	double alpha = opts_pen.alpha, eps = opts.eps;
	
	MatrixXd beta = MatrixXd::Constant(p,nlam,0), Cbeta;
	VectorXd Anew, ajnorm, r = Y, Cnew0, Cnew;
	Anew.setZero(p);
	ajnorm = Anew;
	if(pz){
		Cbeta.setZero(pz,nlam);
		Cnew0.setZero(pz);
	}
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step=0;
		while (step<max_step){
			step++;
			active = 0;
			diffmax = 0;
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active++;	
			if (active>dfmax)
				return List::create(Named("beta") = NULL,Named("Cbeta") = NULL);					
			if(pz){
				Cnew = W.transpose()*r/n;
				r -= W*Cnew;
				Cnew += Cnew0;	
				Cnew0 = Cnew;				
			}				
			for (j = 0; j < p; j++) {
				gj = r.dot(Z.col(j))/n + Anew[j];
				ajnew = penalties(gj, 1, lambda1, alpha, gamma, penalty);
				diff = ajnew-Anew[j];
				if(diff!=0){ 
             		r -= Z.col(j)*diff;	
					Anew[j] = ajnew;
					ajnorm[j] = fabs(ajnew);
					diff = fabs(diff); 
					if(diff>diffmax) diffmax = diff;
				}				
			}			
			if(diffmax<eps) break;
		}//end while
		if(pz){
			likhd[l] = (Y - W*Cnew- Z * Anew).squaredNorm();
			Cbeta.col(l) = Cnew;
		}
		else  likhd[l] = (Y - Z * Anew).squaredNorm();
		for (j = 0; j<p; j++) if (ajnorm[j]) activeA(j,l) = 1;
		beta.col(l) = Anew;			
	}//end for	
	return List::create(Named("beta") = beta,Named("Cbeta") = Cbeta);
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
List MVR_blockwise(MatrixXd Y, MatrixXd Z, MatrixXd W, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd, VectorXi lengths)
{
/*
	Input:
	Y is n*q matrix
	Z is n*p matrix, which is standardized
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is q*p coefficient metrix
*/  
	int l,j,d,dj, active, step, nlam = opts_pen.nlam, n = opts.n, q = opts.q, p = opts.p, p0 = lengths.sum();
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, penalty = opts_pen.pen, pz=opts.pz;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma_pen;
	double lambda1,diffmax,diffnorm;

    MatrixXd Anew=MatrixXd::Constant(p0, q, 0), Bnew, Cbeta, Cnew0, Cnew;
	MatrixXd betapath = MatrixXd::Constant(p0*q, nlam, 0), gj, diff, ajnew;
	VectorXd ajnorm;
	ajnorm.setZero(p);
	if(pz){
		Cbeta.setZero(pz*q,nlam);
		Cnew0.setZero(pz,q);
	}	
	else{
		Cnew0.setZero(1,q);
		Cbeta.setZero(q,nlam);
	}	
	MatrixXd r = Y;
	MatrixXd V, L, Vnorm, Vnew;
	Vnew = MatrixXd::Constant(n, p0, 0);
	List Gamma_sqrtn(p);
    	
	for(d = 0, j=0;j<p;j++){	
        dj = lengths[j];	
		V = Z.middleCols(d,dj);
		L = (V.transpose()*V/n).llt().matrixL();
		Gamma_sqrtn[j] = UpTriangularInv(L.transpose());
		Vnew.middleCols(d,dj) = QbyR(V, Gamma_sqrtn[j],1);
        d += dj;		
	}	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step = 0;
		while (step<max_step) {
			step++;
			active = 0;			
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active++;
			if (active>dfmax) 
				return List::create(Named("betapath") = NULL,Named("Cbeta") = NULL);
			if(pz){
				Cnew = W.transpose()*r/n;
				r -= W*Cnew; 
				Cnew += Cnew0;
				Cnew0 = Cnew;
			}			
			diffmax=0;
			for (d = 0, j = 0; j < p; j++) {
				dj = lengths[j];
				Vnorm = Vnew.middleCols(d,dj);
				gj = Vnorm.transpose()*r/n + Anew.middleRows(d,dj);					
				ajnew = update_blockj(gj, lambda1, alpha, gamma, penalty);
				diff = ajnew - Anew.middleRows(d,dj);
				diffnorm = diff.norm();
				if(diffnorm!=0){
					r -= Vnorm*diff;
					Anew.middleRows(d,dj) = ajnew;
					ajnorm[j] = ajnew.norm();
					if(diffnorm>diffmax) diffmax = diffnorm;
				}
				d+=dj;
			}
			if(diffmax<eps) break;
		}//end while
		for (j = 0; j<p; j++) 	if (ajnorm[j]) activeA(j,l) = 1;
		if(pz){	
			likhd[l] = (Y - W*Cnew- Z * Anew).squaredNorm();
			Cnew.resize(pz*q,1);
			Cbeta.col(l) = Cnew;				
		}	
		else  likhd[l] = (Y - Z * Anew).squaredNorm();		
		
		for(d = 0, j=0; j<p; j++){
			dj = lengths[j];
			Anew.middleRows(d,dj) = (QbyR(Anew.middleRows(d,dj),Gamma_sqrtn[j],0));
			d+=dj;
		}
		Bnew = Anew;
		Bnew.resize(p0*q,1);
		betapath.col(l) = Bnew;
	}//end for	
	return List::create(Named("betapath") = betapath, Named("Cbeta") = Cbeta);
}
//***-------------------------------------------------------------**
//***-------update the jth row of matrix A with penalty-----------**
List MVR_glasso(VectorXd Y, MatrixXd Z, MatrixXd W, MatrixXi &activeA, VectorXd lambda, VectorXd &likhd, VectorXi lengths)
{
/*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
	
	Output:
	Anew is p-dimensional coefficient
*/  
	int l,j,d,dj, active, step, nlam = opts_pen.nlam, n = opts.n, p = opts.p, p0 = lengths.sum();
    double lambda1,diffmax, diffnorm;
	int dfmax = opts_pen.dfmax, max_step = opts.max_step, penalty = opts_pen.pen, pz=opts.pz;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma_pen;
	
	MatrixXd beta = MatrixXd::Constant(p0,nlam,0), Vnorm, Cbeta;
	VectorXd Anew, ajnorm, r = Y, ajnew, gj, diff, Cnew0, Cnew;
	Anew.setZero(p0);
	ajnorm.setZero(p);
	if(pz){
		Cbeta.setZero(pz,nlam);
		Cnew0.setZero(pz);
	}	
	
	for (l = 0; l < nlam; l++) {
		lambda1 = lambda[l];
		step=0;
		while (step<max_step){
			step++;
			active = 0;
			diffmax = 0;
			for (j = 0; j < p; j++)
				if (ajnorm[j] != 0) active++;
			if (active>dfmax)	
                return List::create(Named("beta") = NULL,Named("Cbeta") = NULL);
			if(pz){
				Cnew = W.transpose()*r/n;
				r -= W*Cnew;
				Cnew += Cnew0;	
				Cnew0 = Cnew;				
			}			
			for (d = 0,j = 0; j < p; j++) {
				dj = lengths[j];
                Vnorm = Z.middleCols(d,dj);				
				gj = Vnorm.transpose()*r/n + Anew.segment(d,dj);
				ajnew = update_colj(gj, lambda1, alpha, gamma, penalty);
				diff = ajnew - Anew.segment(d,dj);
				diffnorm = diff.norm();			
				if(diffnorm!=0){ 
             		r -= Vnorm*diff;	
					Anew.segment(d,dj) = ajnew;
					ajnorm[j] = ajnew.norm();		
					if(diffnorm>diffmax) diffmax = diffnorm;
				}	
                d+=dj;					
			}						
			if(diffmax<eps) break;
		}//end while
		if(pz){
			likhd[l] = (Y - W*Cnew- Z * Anew).squaredNorm();
			Cbeta.col(l) = Cnew;
		}
		else  likhd[l] = (Y - Z * Anew).squaredNorm();		
		for (j = 0; j<p; j++) if (ajnorm[j]) activeA(j,l) = 1;
		beta.col(l) = Anew;
	}//end for	
	return List::create(Named("beta") = beta,Named("Cbeta") = Cbeta);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_colwise(MatrixXd Y, MatrixXd Z, MatrixXd W, VectorXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();

	int nlam=opts_pen.nlam, p=opts.p;
	opts.n=Y.rows();
	MatrixXi df = MatrixXi::Constant(p,nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	List fit;
	fit = MVR_colwise(Y, Z, W, df, lambda, likhd);
	return List::create(Named("betapath") = fit[0], Named("df") = df, Named("likhd") = likhd, Named("Cpath") = fit[1]);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_lasso(MatrixXd Y, MatrixXd Z1, MatrixXd W, MatrixXd lambda, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.rows();

	int j,k,nlam=opts_pen.nlam, p=opts.p, q=opts.q, n=Y.rows(), pz=opts.pz;
	opts.n=n;
	MatrixXd betapath = MatrixXd::Constant(p*q, nlam, 0), Z = Z1, beta, Cnew, Cpath;
	MatrixXi df = MatrixXi::Constant(p*q,nlam, 0), activeA;
	VectorXd znorm, likhd0;
	MatrixXd likhd = MatrixXd::Constant(q, nlam, 0);
	likhd0.setZero(nlam);
	List fit;
	if(pz)	Cpath.setZero(pz*q,nlam);
	else Cpath.setZero(q,nlam);

    znorm = Z1.colwise().norm()/sqrt(n);
	for(j=0;j<p;j++) Z.col(j) = Z1.col(j)/znorm[j];
	for(k=0; k<q; k++){
        activeA.setZero(p,nlam);		
		fit = MVR_lasso(Y.col(k), Z, W, activeA, lambda.col(k), likhd0);
		beta = fit[0];
		Cnew = fit[1];
		df.block(k*p, 0, p, nlam) = activeA;
		for(j=0;j<p;j++) beta.row(j) /= znorm[j];
		betapath.middleRows(k*p, p) = beta;
		likhd.row(k) = likhd0;
		if(pz) Cpath.middleRows(k*pz,pz) = Cnew;
	}
	return List::create(Named("betapath") = betapath, Named("df") = df, Named("likhd") = likhd, Named("Cpath") = Cpath);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_blockwise(MatrixXd Y, MatrixXd Z, MatrixXd W, VectorXd lambda, VectorXi lengths, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();

	int nlam=opts_pen.nlam, p=opts.p;
	opts.n=Y.rows();
	MatrixXi df = MatrixXi::Constant(p,nlam, 0);
	VectorXd likhd = VectorXd::Constant(nlam, 0);
	List fit;
   
	fit = MVR_blockwise(Y, Z, W, df, lambda, likhd, lengths);
	return List::create(Named("betapath") = fit[0], Named("df") = df, Named("likhd") = likhd, Named("Cpath") = fit[1]);
}
//----------------------------------------------------------------**
//***----------Estimation Multivariate with  penalty--------------**
// [[Rcpp::export]]
List EstMVR_glasso(MatrixXd Y, MatrixXd Z, MatrixXd W, MatrixXd lambda, VectorXi lengths, List optsList, List optsList_pen){
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.n = as<int>(optsList["n"]);
	opts.p = as<int>(optsList["p"]);
	opts.q = as<int>(optsList["q"]);
	opts.pz = as<int>(optsList["pz"]);

	opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma_pen = as<double>(optsList_pen["gamma_pen"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.rows();

	int j,k,d,dj,nlam=opts_pen.nlam, p=opts.p, q=opts.q, n=Y.rows(), p0 = lengths.sum(), pz=opts.pz;
	opts.n=n;
	MatrixXd betapath = MatrixXd::Constant(p0*q, nlam, 0), beta, Cbeta, Cpath;
	MatrixXi df = MatrixXi::Constant(p*q,nlam, 0), activeA;
	MatrixXd likhd = MatrixXd::Constant(q, nlam, 0);
	VectorXd likhd0;
	likhd0.setZero(nlam);
	List fit;
	if(pz)	Cpath.setZero(pz*q,nlam);
	else Cpath.setZero(q,nlam);	

	MatrixXd V, L, Vnew = MatrixXd::Constant(n, p0, 0);
	List Gamma_sqrtn(p);
	for(d=0,j=0;j<p;j++){	
        dj = lengths[j];
		V = Z.middleCols(d,dj);
		L = (V.transpose()*V/n).llt().matrixL();
		Gamma_sqrtn[j] = UpTriangularInv(L.transpose());
		Vnew.middleCols(d,dj) = QbyR(V, Gamma_sqrtn[j],1);
        d += dj;		
	}
	for(k=0; k<q; k++){
        activeA.setZero(p,nlam);		
		fit = MVR_glasso(Y.col(k), Vnew, W, activeA, lambda.col(k), likhd0, lengths);
		beta = fit[0];
		Cbeta = fit[1];		
		for(d=0,j=0; j<p; j++){
			dj = lengths[j];
			beta.middleRows(d,dj) = QbyR(beta.middleRows(d,dj),Gamma_sqrtn[j],0);
			d += dj;
		}			
		betapath.middleRows(k*p0,p0) = beta;
		likhd.row(k) = likhd0.transpose();
		df.middleRows(k*p,p) = activeA;	
		if(pz) Cpath.middleRows(k*pz,pz) = Cbeta;
	}
	return List::create(Named("betapath") = betapath, Named("df") = df, Named("likhd") = likhd, Named("Cpath") = Cpath);
}
