\name{ipr}
\alias{ipr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterative Proportional Repartition (IPR) algorithm
}
\description{
Estimating the health cost repartition among diseases in the presence of multimorbidity, i.e. when some patients have multiple diseases. Using the Iterative Proportional Repartition algorithm (see reference below), the goal is to estimate the average cost for each disease, starting from the global health costs available for each patient.
}
\usage{
ipr(X, y, print.it=FALSE, start=rep(1,dim(X)[2]), cutup=Inf, cutlow=cutup,
epsrel=0.001, epsabs=0.1, maxiter=1000, det=FALSE)
}
\arguments{
  \item{X}{Matrix with \eqn{x_{ij}=1} if patient \eqn{i} suffers from disease \eqn{j} and \eqn{x_{ij}=0} otherwise. Each row thus refers to one patient and each column to one disease. The number of columns of \code{X} corresponds to the number of diseases considered.}
  \item{y}{Vector where \eqn{y_i} is the global health cost of patient \eqn{i}. The length of \code{y} must be equal to the number of rows of \code{X}.}
  \item{print.it}{Logical. If \code{TRUE}, the number of the current iteration and the current estimates are printed.}
  \item{start}{Vector of initial estimates of the average cost for each disease to start IPR algorithm. Default is an initial average cost of 1 for all diseases. The length of \code{start} must be equal to the number of columns of \code{X}.}
  \item{cutup, cutlow}{Options which can be used to get a robust version of IPR. If the current allocated cost of disease \eqn{j} for patient \eqn{i} is more than \code{cutup} times more expansive (or less then \code{cutlow} times less expansive) than the current average cost estimate of that disease \eqn{j}, then this outlying allocated cost is not taken into account in the next iteration to compute the average cost of disease \eqn{j}. By default, \code{cutup} and \code{cutlow} are set to \code{Inf}.}
  \item{epsrel}{Stopping criterion such that the IPR algorithm stops if for all diseases, the current estimated average cost differs by less than 100*\code{epsrel} percent from what it was at the previous iteration. The default value is 0.001. Should be set to 0 to ignore that criterion.}
  \item{epsabs}{Stopping criterion such that the IPR algorithm stops if for all diseases, the current estimated average cost differs (in absolute value) by less than \code{epsabs} percent from what it was at the previous iteration. The default value is 0.1. Should be set to 0 to ignore that criterion.}
  \item{maxiter}{Maximal number of iterations of IPR algorithm. The default value is 1000.}
  \item{det}{Logical. If \code{TRUE}, the allocated costs of each disease for each patient are given, by returning a matrix \eqn{Y} where \eqn{y_{ij}} is the estimated cost of disease \eqn{j} for patient \eqn{i}.}
}
\details{
Let us consider \eqn{n} patients and \eqn{p} diseases. We are given a matrix \eqn{X} such that \eqn{x_{ij}=1} if the patient \eqn{i} suffers from disease \eqn{j} and \eqn{x_{ij}=0} otherwise. We are also given a vector \eqn{y}, where \eqn{y_i} is the global health cost of patient \eqn{i}. In order to estimate the average cost of each disease, the IPR algorithm works as follows:

1. Start with some initial estimates \eqn{mu_j}, e.g. \eqn{mu_j}=1 for all \eqn{j=1,\dots,p}. Those initial estimates are stored in the vector \code{start}.

2. Allocate the cost \eqn{y_i} among the diseases diagnosed for patient \eqn{i}, proportionally to the current estimates \eqn{mu_j}.

3. Update the current estimate of \eqn{mu_j} by averaging the specific costs obtained in step 2 for the disease \eqn{j} over the patients having that disease.

4. Repeat steps 2 and 3 until a stopping criterion, based on relative or absolute distance between two consecutive iterations. The stopping criterion can be defined with \code{epsabs} or \code{epsrel}.


By construction, the IPR algorithm satisfies two properties. First, it allows to obtain positive estimates for each average disease cost. Secondly, it allows to retrieve the total health costs. In other words, the sum of the estimates \eqn{mu_j} multiplied by the number of patients suffering from \eqn{j} is equal to the sum of the costs \eqn{y_i}.

The estimate of total cost \eqn{tau_j} spent for disease \eqn{j} as well as the estimated proportion \eqn{pi_j} of the total costs which is allocated to disease \eqn{j} are also returned by our function.

Mathematically, \eqn{tau_j} is the sum over \eqn{i=1} to \eqn{i=n} of \eqn{X_{ij}*mu_j}, while \eqn{pi_j} is defined by \eqn{tau_j} divided by the sum of all \eqn{tau_k}.
}
\value{
\item{coef}{A vector with the estimated average cost of each disease.}
\item{total}{A vector with the estimated total cost spent for each disease.}
\item{proportions}{A vector with the estimated proportion of total cost spent for each disease.}
\item{niter}{The number of iterations of IPR algorithm until the stopping criterion is achieved.}
\item{esprel}{The stopping criterion based on a relative distance between two consecutive iterations which has been used.}
\item{epsabs}{The stopping criterion based on an absolute distance between two consecutive iterations which has been used.}
\item{detail}{A matrix with the allocated costs of each disease for each patient, if \code{det} is set to \code{TRUE}.}
}
\references{
Rousson, V., Rossel, J.-B. & Eggli, Y. (2019). Estimating Health Cost Repartition Among Diseases in the Presence of Multimorbidity. \emph{Health Services Research and Managerial Epidemiology, 6.}

Rossel, J.-B., Rousson, V. & Eggli, Y. A comparison of statistical methods for allocating disease costs in the presence of interactions. In preparation.
}
\author{
Dr. Jean-Benoit Rossel (\email{jean-benoit.rossel@unisante.ch}), Prof. Valentin Rousson and Dr. Yves Eggli.
}
\examples{
# Here is a first example with 10 patients and 4 diseases:
X <- matrix(c(1,0,0,0,
0,1,1,0,
0,1,0,1,
1,0,0,1,
1,1,1,0,
0,0,1,1,
0,1,0,0,
1,1,0,0,
0,1,1,1,
0,0,0,1),ncol=4,byrow=TRUE)

y <- c(500,200,100,400,1000,500,100,300,800,2000)

# If we would use a linear model without intercept to estimate the average
# disease costs, we would obtain a negative value for disease 2.
lm(y~X-1)

# The IPR algorithm provides only positive estimates
ipr(X,y)


# Here is a second example:
X <- matrix(c(1,0,0,1,1,1),nrow=3,byrow=TRUE)
y <- c(5000,500,6600)

# We have three patients. The first one has only disease 1 with a cost of 5000.
# The second one has only disease 2 with a cost of 500 (i.e. ten times less
# expansive than disease 1). The third patient has both diseases with
# a cost of 6600 (i.e. 5000 + 500 + an extra cost of 1100).

# Using a linear model, one would allocate the extra cost equally between
# the three patients. The estimated average cost would thus be 5000+(1100/3)
# for disease 1 and 500+(1100/3) for disease 2.
lm(y~X-1)

# Using IPR algorithm, one allocates the extra cost taking into account that
# disease 1 is ten times more expansive than disease 2 when occuring alone.
# One thus gets an estimated average cost of 5500 for disease 1 and
# of 550 for disease 2.
ipr(X,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{methods}% __ONLY ONE__ keyword per line
