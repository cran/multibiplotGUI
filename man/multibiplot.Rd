\name{multibiplot}
\alias{multibiplot}
\title{Multibiplot Analysis in R}
\description{A GUI with which users can construct and interact
        with Multibiplot Analysis and provides inferential results by using Bootstrap Methods.}
\usage{
multibiplot(x, ni)
}

\arguments{
\item{x}{A data frame in which the different datasets have been juxtaposed}
\item{ni}{A vector specifying the length of each dataset}
}

\details{
When the function is launched with the appropriate arguments,  
firstly, it is necessary to select the kind of data to be analyzed. Then, an option window is displayed where 
you can change the color, the size, the label and/or the symbol of an element or 
of a set of elements; to select the 
transformation data, to select the kind of Biplot factorization to be applied, to change the window size containing the graphs and to tick the checkbox to show the axes in the graph.
Press the Graph button and then choose the number of axes to be retained.
When the graph will be shown, the function will allow you
 to change characteristics of the points with the mouse.
Press the right mouse button and a window will be displayed to change the color, the size,
the label and/or the symbol of the nearest point of position clicked.
Press the left mouse button and a window will be displayed to select one option:
Change the position label, Remove label or Do nothing.
It is also possible to select the dimensions shown in the graph and to change the limits of the axes.
In the window there are six menus with their corresponding submenus:

\itemize{
        \item File
	\itemize{
		\item Copy image
		\item Save image
		\itemize{
			\item PDF file
			\item Eps file
			\item Png file
			\item Jpg/Jpeg file
		}
		\item Exit
	}
	\item 3D
	\itemize{
		\item 3D
	}
	\item Projections
	\itemize{
		\item Variables
		\item Back to original data
	}
	\item Options
	\itemize{
		\item Change title
		\item Show/Hide axes
	}	
	\item Cluster
	\itemize{
		\item Hierarchical cluster with biplot coordinates
		\item K-means with biplot coordinates
		\item K-medoids with biplot coordinates
		\item Back to original graph
	}
	\item Bootstrap
	\itemize{
		\item Bootstrap
		}
}

The File menu provides different options to save the graph and permits to exit the program. The second menu shows the graph in 3 dimensions. The third menu allows the user to project the individuals onto the direction representing one variable selected from a listbox. This menu permits to go back to original graph. The following menu permits to change the title and to show/hide the axes in the graph. The fifth menu allows the user to analyze the biplot coordinates with clustering techniques. The last menu provides inferential results using Bootstrap methods. It is possible to select the number of resamples to be extracted and the conficende level to calculate the intervals presented in the results. This menu saves a file containing the results and graphs with histograms and QQ-plots generated with the bootstrap replications. 
}

\value{A graph showing the data representation, an output file containing the 
contributions, qualities of representation, goodness of fit, coordinates and eigen values and output files containing these results in a classical and inferential form.}

\references{
\itemize{

\item Gabriel, K. R. (1971). The Biplot graphic display of matrices with applications to principal components analysis. Biometrika, 58(3), 453-467.

\item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.

\item Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Annals of Statistics, 7, 1-26.

\item Efron, B. (1987). Better bootstrap confidence intervals. Journal of the American Statistical Association, 82, 171-185.

\item Efron, B., & Tibshirani, R. J. (1993). An introduction into the bootstrap. New York: Chapman and Hall.

\item Baccala, N. (2004). \emph{Contribuciones al analisis de matrices de datos multivia: tipologia de las variables}. Phd Thesis.}
}

\author{Ana Belen Nieto Librero \email{ananieto@usal.es}, 
Nora Baccala \email{norabaccala@gmail.com}, 
Purificacion Vicente Galindo \email{purivg@usal.es}, 
Purificacion Galindo Villardon \email{pgalindo@usal.es}}

\examples{
data(iris)
multibiplot(iris[,-5], c(50,50,50))
}
\keyword{multivariate}

