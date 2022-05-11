\name{multibiplot}
\alias{multibiplot}
\title{Multibiplot Analysis in R}
\description{Provides a GUI with which users can construct and interact
        with Multibiplot Analysis}
\usage{
multibiplot(x, ni)
}

\arguments{
\item{x}{A data frame in which the different sets data have been juxtaposed}
\item{ni}{A vector specifying the length of each set data}
}

\details{
see documentation
When the function is launched,
firstly, it is necessary to enter the number of tables to intervent in the analysis
and choose a option in the dropdown menu. Then, an options window is displayed where
you can change the color, the size, the label and/or the symbol of an element or
of a set of elements; to select the
transformation data and to tick the checkbox to show the axes in the graph.
Press the Graph button and then choose the number of axes to be retained.
When the graph will be shown the function will allow you
 to change characteristics of the points with the mouse.
Press the right mouse button and a window will be displayed to change the color, the size,
the label and/or the symbol of the nearest point of position clicked.
Press the left mouse button and a window will be displayed to select one option:
Change the position label, Remove label or Do nothing.
Press the 3D button and a window will be displayed with the 3D-graph.

}
\value{A graph showing the data representation and an output file containing the
contributions, qualities of representation, goodness of fit, coordinates and eigen values}

\references{Baccala N (2004). \dQuote{Contribuciones al analisis de matrices de datos multivia:
Tipologia de las variables.} Master's thesis, Universidad de Salamanca.}

\author{Ana Belen Nieto Librero \email{ananieto@usal.es},
Nora Baccala \email{norabaccala@gmail.com},
Purificacion Vicente Galindo \email{purivg@usal.es},
Purificacion Galindo Villardon \email{pgalindo@usal.es}}

\examples{
data(iris)
multibiplot(iris[,-5], rep(50,3))
}
\keyword{multivariate}
