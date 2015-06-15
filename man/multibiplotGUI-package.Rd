\name{multibiplotGUI-package}
\alias{multibiplotGUI-package}
\alias{multibiplotGUI}
\docType{package}
\title{Multibiplot Analysis package}
\description{A GUI with which users can construct and interact
        with Multibiplot Analysis and provides inferential results by using Bootstrap Methods.}
        
\details{
\tabular{ll}{
Package: \tab multibiplotGUI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-15\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Ana Belen Nieto Librero \email{ananieto@usal.es}, 
Nora Baccala \email{norabaccala@gmail.com}, 
Purificacion Vicente Galindo \email{purivg@usal.es}, 
Purificacion Galindo Villardon \email{pgalindo@usal.es}

Maintainer: Ana Belen Nieto Librero <ananieto@usal.es>
}
\references{
\itemize{

\item Gabriel, K. R. (1971). The Biplot graphic display of matrices with applications to principal components analysis. Biometrika, 58(3), 453-467.

\item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.

\item Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Annals of Statistics, 7, 1-26.

\item Efron, B. (1987). Better bootstrap confidence intervals. Journal of the American Statistical Association, 82, 171-185.

\item Efron, B., & Tibshirani, R. J. (1993). An introduction into the bootstrap. New York: Chapman and Hall.

\item Baccala, N. (2004). \emph{Contribuciones al analisis de matrices de datos multivia: tipologia de las variables}. Phd Thesis.
}}
\keyword{ package }

\examples{
data(iris)
multibiplot(iris[,-5], c(50,50,50))
}
