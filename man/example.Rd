\name{example}
\docType{data}
\alias{example}
\title{matrices and simulated data for small example}
\description{
This data set contains the 3 by 2 matrix of contrasts C, the 150 by 3 design matrix D,
the 150 by 5 data matrix X and the 5 by 5 inverse covariace matrix sinv (as a sparse symmetricMatrix)
Library(matrix) is required to work with sinv. The 150 by 5 matrix M is the true mean structure of the data.
The simulated data X comes from a normal distribution wth mean M and the covariance matrix of each row is sinv^-1.
For variables 2,4 and 5 the mean of observations 1-50 is 1, 51-100 is 2 and 101:150 is 3. The
remaining variables have mean zero.
}
\usage{example}
\format{matrices C,D,X,M,sinv}
\examples{
data(example)
dim(C)
}
\source{Simulated data with known structure}
\references{
Kiiveri (2001)) Multivariate analysis of microarray data: differential expression and differential connection.
}
\keyword{datasets}