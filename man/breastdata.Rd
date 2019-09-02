\name{breastData}
\alias{breastData}
\docType{data}
\title{
Breast cancer gene expression and DNA copy number dataset 
}
\description{
The breast cancer dataset includes gene expressions and comparative genomic
hybridization measurements for 89 subjects, which is from Chin et al. (2006). This dataset has been considered by  Witten et al. (2009) and Chen et al. (2013). In our paper, we selected chromosome 21, including \eqn{q = 44} variables for copy-number variations and \eqn{p = 227} variables for gene expression. As in Chen et al. (2013), we consider copy-number variations
as the responses and gene expressions as the predictors.
}
\usage{
data(breastData)
}
\details{
The "breastData" is formated as a list with elements:

\tabular{ll}{
dna: \tab the CGH spots , a matrix with size \eqn{2149\times 89} and the smaple size 89\cr

rna: \tab genes,  a matrix with size \eqn{19672\times 89} and the smaple size 89\cr

chrom: \tab chromosomal location of each CGH spot, a 2149-vector\cr

nuc: \tab nucleotide position of each CGH spot, a 2149-vector\cr

gene: \tab an accession number for each gene, a 19672-vector\cr

genenames: \tab gene name, a 19672-vector\cr

genechr: \tab chromosomal location of each gene, a 19672-vector\cr

genedesc: \tab description of each gene, a 19672-vector\cr

genepos: \tab nucleotide position of each gene, a 19672-vector\cr
}
}

\references{
Chin, K., DeVries, S., Fridlyand, J., Spellman, P., Roydasgupta, R., Kuo,W.-L., Lapuk, A., Neve,
R., Qian, Z., Ryder, T., Chen, F., Feiler, H., Tokuyasu, T., Kingsley, C., Dairkee, S., Meng, Z., Chew,
K., Pinkel, D., Jain, A., Ljung, B., Esserman, L., Albertson, D.,Waldman, F. & Gray, J. (2006). Genomic and transcriptional aberrations linked to breast cancer pathophysiologies. Cancer cell {\bold{10 (6)}}, 429-541.


Witten, D. M., Tibshirani, R. and Hastie, T. (2009). A penalized matrix decomposition, with applications to sparse
principal components and canonical correlation analysis. Biostatistics {\bold{10 (3)}}, 515-534.


Chen, K., Dong, H., and Chan, K. S. (2013). Reduced rank regression via adaptive
nuclear norm penalization. Biometrika, {\bold{100 (4)}}, 901-920.
}
\examples{
  data(breastData)
  attach(breastData)
  Y = t(dna[chrom==21,])
  Xt = t(rna[which(genechr==21),])
  n = nrow(Y)

  minX = apply(Xt,2,min)
  maxX = apply(Xt,2,max)
  X = (Xt - matrix(rep(minX,each = n),n))/matrix(rep(maxX-minX,each = n),n)
  Y = scale(Y)
  fit <- mqr_sparse_dr(Y[,1:5], X[,1:10])
  D3hat <- fit$Dnew
  opt <- fit$rk_opt
  detach(breastData)
}
\keyword{datasets}
