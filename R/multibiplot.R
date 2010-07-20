multibiplot <-
function()
{
	
#############################################################################
#########	libraries
#############################################################################

	library(RODBC)
	require(tcltk)
	library(tkrplot)
	library(rgl)
	tclRequire("BWidget")


#############################################################################
### Informative window
#############################################################################
symbols <- c("*",".", "o","O","0","+","-","|","%","#")
tipo<<-"RCMP"
winfor<-tktoplevel()
tkwm.title(winfor,"Multibiplot")
OnOKinf <- function()
			{

#############################################################################
#########	Window to enter the number of matrices to analyze
#############################################################################
	
	tkdestroy(winfor)
	wnummat<-tktoplevel()


	OnOK <- function()
	{
		NameVal <<- tclvalue(nummat)
		tkdestroy(wnummat)
		NameVal<<-as.integer(NameVal)
		msg <- paste(gettext("You have selected ",domain="R-multibiplotGUI"),NameVal,gettext("matrices to analyze",domain="R-multibiplotGUI"))
		tkmessageBox(message=msg)



#############################################################################
#########	We create vectors with the names of the matrices and of the  
#########	eigen values
#############################################################################

		nombres<<-rep("X",NameVal)
		nombresslab<<-rep("X",NameVal)
		nombresst<<-rep("Xstand",NameVal)
		valorespropios<<-rep(1,NameVal)
		dimensiones<<-rep(1,NameVal)
		diagonales<<-rep("N",NameVal)
		invdiagonales<<-rep("N",NameVal)


		for (i in 1:NameVal)
		{
			nombres[i]<<-paste("X",i, sep = "")
			nombresslab[i]<<-paste("Xs",i, sep = "")
			nombresst[i]<<-paste("Xstand",i, sep = "")
			valorespropios[i]<<-paste("lambda",i, sep = "")
			dimensiones[i]<<-paste("d",i, sep = "")
			diagonales[i]<<-paste("N",i, sep = "")
		}


#############################################################################
#########	Data
#############################################################################

		for (j in 1:NameVal)
		{
			file <- tclvalue(tkgetOpenFile())
			if (!length(file)) return()
			channel <- odbcConnectExcel(file)
			hojastabla<-sqlTables(channel)$TABLE_NAME
	
			whoja<-tktoplevel()
			tkwm.title(whoja,gettext("Select one table",domain="R-multibiplotGUI"))
			tl<-tklistbox(whoja,height=4,selectmode="single",background="white")
			tkgrid(tklabel(whoja,text=gettext("Select one table:",domain="R-multibiplotGUI")))
			tkgrid(tl)
	
			for (k in (1:(dim(sqlTables(channel))[1])))
			{
    				tkinsert(tl,"end",hojastabla[k])
			}

			tkselection.set(tl,0)  #Indexing starts at zero.

			OnOK <- function()
			{
				hojaChoice <<- hojastabla[as.numeric(tkcurselection(tl))+1]
				matrizd<<-sqlFetch(channel,hojaChoice)
				assign(nombres[j],matrizd,envir = .GlobalEnv)
				odbcClose(channel)
    				tkdestroy(whoja)
    			}
		
			OK.but <-tkbutton(whoja,text="   OK   ",command=OnOK)
			tkgrid(OK.but)
			tkfocus(whoja)
			tkwait.window(whoja)

		}



###############################################################################				
####	We create vectors of matrices without labels
###############################################################################


		for (i in 1: NameVal)
		{
		
			matriz<-array(dim=c(dim(get(nombres[i]))[1],(dim(get(nombres[i]))[2]-1)))
			for(k in 2:dim(get(nombres[i]))[2])
			{
				for(j in 1:dim(get(nombres[i]))[1])
				{
					matriz[j,k-1]<-(get(nombres[i]))[j,k]
				}
			}

			assign(nombresslab[i],matriz,envir= .GlobalEnv)
			
		}

###############################################################################				
####	If the selected option is Some sets of individuals observed in a single 
####	 set of variables
###############################################################################
	
		if ((filas == 1) && (columnas == 0)){
			for (z in 1:NameVal)
			{
			
				assign(dimensiones[z],(dim(get(nombresslab[z])))[1],envir= .GlobalEnv)

				if (z == 1) {
					X<<-get(nombresslab[z])
				} else {
					X<<-rbind(X,get(nombresslab[z]))
				}	
			}

###############################################################################				
####	We standardize the matrices
###############################################################################

			media<<-colMeans(X)
			desvt<<-sqrt(diag(var(X)))
			
			for(i in 1:NameVal)
			{
				matrizst<<-get(nombresslab[i])
				for(k in 1:(dim(get(nombresslab[i]))[2]))
				{	
					for(z in 1:dim(get(nombresslab[i]))[1])
					{
						matrizst[z,k]<<-(get(nombresslab[i])[z,k]-media[k])/desvt[k]
					}
				}

				assign(nombresst[i],matrizst,envir = .GlobalEnv)
				matrizcor<<-cor(matrizst)
				vvpropios<<-eigen(matrizcor)
				assign(valorespropios[i],vvpropios$values[1],envir = .GlobalEnv)
				

				diagonal<<-rep(vvpropios$values[1],dim(get(nombresslab[i]))[1])
				assign(diagonales[i],diagonal,envir = .GlobalEnv)


				if (i == 1) {
					textindividuos<<-t((get(nombres[i]))[1])
					textvariables<<-names(get(nombres[i]))[2:dim(get(nombres[i]))[2]]
					Xpon<<-(get(nombresst[i]))/get(valorespropios[i])
				} else {
					textindividuos<<-cbind(textindividuos,t((get(nombres[i]))[1]))
					Xpon<<-rbind(Xpon,(get(nombresst[i]))/get(valorespropios[i]))
				}	
			}

		
		} else{
			for(i in 1:NameVal)
			{
				matrizcor<<-cor(get(nombresslab[i]))
				vvpropios<<-eigen(matrizcor)
				assign(valorespropios[i],vvpropios$values[1],envir = .GlobalEnv)
				assign(dimensiones[i],(dim(get(nombresslab[i])))[2],envir= .GlobalEnv)

				if (i == 1) {
					textindividuos<<-t((get(nombres[i]))[1])
					textvariables<<-names((get(nombres[i])))[2:dim(get(nombres[i]))[2]]
					Xpon<<-(get(nombresslab[i]))/get(valorespropios[i])
				} else {
					textvariables<<-cbind(t(textvariables),t(names(get(nombres[i]))[2:dim(get(nombres[i]))[2]]))
					Xpon<<-cbind(Xpon,(get(nombresslab[i]))/get(valorespropios[i]))
					
				}	
			}	
		}
		
##############################################################################
####	We create vectors of the colors
##############################################################################

		colvariables<<-rep("blue",times = dim(Xpon)[2])		
		colindividuos<<-rep("green",times = dim(Xpon)[1])


##############################################################################
####	We create vectors of the symbols
##############################################################################

		simvariables<<-rep(" ",times = dim(Xpon)[2])		
		simindividuos<<-rep("+",times = dim(Xpon)[1])



	
##############################################################################
#####	Window to change labels and colors and select the biplot 
##############################################################################

		tt<-tktoplevel()
		tkwm.title(tt,gettext("Options",domain="R-multibiplotGUI"))

#####Dropdown menu#############################

		topMenutt <- tkmenu(tt)
		tkconfigure(tt,menu=topMenutt)
		fileMenutt <- tkmenu(topMenutt,tearoff=FALSE)
		
		if ((filas == 1) && (columnas == 0)){
			tkadd(fileMenutt,"command",label="RMP-biplot",command=function() tipo<<-"RMP")
			tkadd(fileMenutt,"command",label="RCMP-biplot",command=function() tipo<<-"RCMP")
		} else {
			tkadd(fileMenutt,"command",label="RMP-biplot",command=function() tipo<<-"RMP")
			tkadd(fileMenutt,"command",label="RCMP-biplot",command=function() tipo<<-"RCMP")
		}
		tkadd(topMenutt,"cascade",label="Biplot",menu=fileMenutt)


#### Frames

	framet<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")

	frametext<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")


	frameok<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")

	framecol<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")

	framename<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")

	frames<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")

	framegraphic<-tkframe(tt, relief = "flat", borderwidth = 2, 
        background = "white")




#### List of transformations

	scrtrans <- tkscrollbar(framet, repeatinterval=5,
				   command=function(...)tkyview(tltrans,...))
	tltrans<-tklistbox(framet,height=6,width=42,selectmode="single",yscrollcommand=function(...)tkset(scrtrans,...),background="white")

	tipostd <- c(gettext("Subtract the global mean",domain="R-multibiplotGUI"),gettext("Column centering",domain="R-multibiplotGUI"),
  gettext("Standardize columns",domain="R-multibiplotGUI"),gettext("Row centering",domain="R-multibiplotGUI"),
  gettext("Standardize rows",domain="R-multibiplotGUI"),gettext("Raw data",domain="R-multibiplotGUI"))
	for (i in (1:6))
	{
	    tkinsert(tltrans,"end",tipostd[i])
	}
	tkselection.set(tltrans,2)  # Default Standardize Columns.
	OnOKtrans <- function()
			{
				tChoice <<- tipostd[as.numeric(tkcurselection(tltrans))+1]
	
### We standardize the matrix 
		if (tChoice==gettext("Subtract the global mean",domain="R-multibiplotGUI")){
			media<<-mean(Xpon)
			Xstd<<-Xpon




			for(j in 1:dim(Xpon)[2])
				{for(i in 1:dim(Xpon)[1])
					{
						Xstd[i,j]<<-(Xpon[i,j]-media)		
					}
				}
		}


		if (tChoice==gettext("Column centering",domain="R-multibiplotGUI")){
			mediav<<-colMeans(Xpon)
			Xstd<<-Xpon




			for(j in 1:dim(Xpon)[2])
				{for(i in 1:dim(Xpon)[1])
					{
						Xstd[i,j]<<-(Xpon[i,j]-mediav[j])	
					}
				}
		}


		if (tChoice==gettext("Standardize columns",domain="R-multibiplotGUI")){
			mediav<<-colMeans(Xpon)
			desvvar<<-sqrt(diag(var(Xpon)))
			Xstd<<-Xpon




			for(j in 1:dim(Xpon)[2])
				{for(i in 1:dim(Xpon)[1])
					{
						Xstd[i,j]<<-(Xpon[i,j]-mediav[j])/desvvar[j]
					}
				}

		}


		if (tChoice==gettext("Row centering",domain="R-multibiplotGUI")){
		
			mediav<<-rowMeans(Xpon)
			Xstd<<-Xpon


			for(j in 1:dim(Xpon)[2])
				{for(i in 1:dim(Xpon)[1])
					{
						Xstd[i,j]<<-(Xpon[i,j]-mediav[i])
					}
				}
		}


		if (tChoice==gettext("Standardize rows",domain="R-multibiplotGUI")){
			mediav<<-rowMeans(Xpon)
			desvvar<<-mediav
			for (i in 1:dim(Xpon)[1])
			{
				desvvar[i]<<-sqrt(var(Xpon[i,]))
			}


			Xstd<<-Xpon


			for(j in 1:dim(Xpon)[2])
				{for(i in 1:dim(Xpon)[1])
					{
						Xstd[i,j]<<-(Xpon[i,j]-mediav[i])/desvvar[i]
					}
				}
		}
	
		if (tChoice==gettext("Raw data",domain="R-multibiplotGUI")){
		
		Xstd<<-Xpon

		}
	}

		OK.buttrans <-tkbutton(frameok,text="    OK    ",command=OnOKtrans)
		tkpack(OK.buttrans,expand = "TRUE", side="right", fill = "both")



##### Checkbox to show the axes or not #######

	cb <- tkcheckbutton(framecol)
	cbValue <- tclVar("0")
	tkconfigure(cb,variable=cbValue)


##############################################################################
##### 	List of individuals
##############################################################################

		indicei<-NULL
		Namei<-NULL
		NameVali<-NULL


		scri <- tkscrollbar(framet, repeatinterval=5,
				   command=function(...)tkyview(tli,...))

		tli<-tklistbox(framet,height=4,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scri,...),background="white")
	
		tkpack(tklabel(frametext,text=gettext("Individuals",domain="R-multibiplotGUI")),side="left",expand = "TRUE",fill="both")



		for (i in 1:(dim(Xpon)[1]))
		{
    			tkinsert(tli,"end",textindividuos[i])
		}
		tkselection.set(tli,0) #  Indexing starts at zero.

		OnOKi <- function()
		{
			Choicei <<- textindividuos[as.numeric(tkcurselection(tli))+1]
	

##### Color of the selected variable  #############

			indicei<<-as.numeric(tkcurselection(tli))+1
			colori <- colindividuos[indicei[1]]
	 		tkconfigure(canvasi,bg=colori)


##### Text of the selected variable  #############

			Namei <<- tclVar(textindividuos[indicei[1]])
			tkconfigure(entry.Namei,textvariable=Namei)
    
		}
	
		OK.buti <-tkbutton(frameok,text="    OK    ",command=OnOKi)
		
		tkpack(tli,scri,expand = "TRUE", side="left", fill = "both")
		tkpack.configure(scri,side="left")

		tkpack(OK.buti,expand = "TRUE", side="left", fill = "both")
		



#######Color#######################################
	
		indicei<-as.numeric(tkcurselection(tli))+1
		colori <- colindividuos[indicei[1]]
		canvasi <- tkcanvas(framecol,width="57",height="20",bg=colori)

		ChangeColori <- function()
		{
 
			colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colindividuos[indicei[1]],title=gettext("Choose a color",domain="R-multibiplotGUI")))
 
			 if (nchar(colori)>0)
    				{
					tkconfigure(canvasi,bg=colori)
		 			colindividuos[indicei]<<-colori
				}
		}

		ChangeColor.buttoni<- tkbutton(framecol,text=gettext("Change Color",domain="R-multibiplotGUI"),command=ChangeColori,width=2)
		tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")

 

######Labels   ###################################
		
		Namei <- textindividuos[indicei[1]]
		entry.Namei <-tkentry(framename,width="10",textvariable=Namei)

		OnOKli <- function()
		{
			NameVali <<- tclvalue(Namei)
			textindividuos[indicei]<<-NameVali

#####Values of listbox###############################
		
			for (i in 1:dim(Xpon)[1])
			{
				tkdelete(tli,0)
			}

			for (i in 1:(dim(Xpon)[1]))
			{
   				tkinsert(tli,"end",textindividuos[i])
			}


		}

		OK.butli <-tkbutton(framename,text=gettext(" Change label",domain="R-multibiplotGUI"),command=OnOKli,width=2)
		tkbind(entry.Namei, "<Return>",OnOKli)
		tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
	
######Symbols  ###################################

		comboBoxi <- tkwidget(frames,"ComboBox",editable=FALSE,values=symbols,width=10)

		chang.symi <- function()
		{
    		simChoicei <<- symbols[as.numeric(tclvalue(tcl(comboBoxi,"getvalue")))+1]
		simindividuos[indicei]<<-simChoicei
		}
		Change.symboli <-tkbutton(frames,text=gettext("   Change symbol   ",domain="R-multibiplotGUI"),command=chang.symi,width=6)
		tkpack(comboBoxi,Change.symboli,side="left",expand="TRUE", fill="both")




##### List of variables ###########################

		indicev<-NULL
		Namev<-NULL
		NameValv<-NULL


		scrv <- tkscrollbar(framet, repeatinterval=5,
				   command=function(...)tkyview(tlv,...))

		tlv<-tklistbox(framet,height=4,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scrv,...),background="white")
	
		tkpack(tklabel(frametext,text=gettext("Variables",domain="R-multibiplotGUI")),side="left",expand = "TRUE",fill="both")


		for (i in 1:dim(Xpon)[2])
		{
    			tkinsert(tlv,"end",textvariables[i])
		}
		tkselection.set(tlv,0) #  Indexing starts at zero.

		OnOKv <- function()
		{
			Choicev <<- textvariables[as.numeric(tkcurselection(tlv))+1]
	
##### Color of the selected variable  #############
		
			indicev<<-as.numeric(tkcurselection(tlv))+1
			colorv <- colvariables[indicev[1]]
	 		tkconfigure(canvasv,bg=colorv)


##### Text of the selected variable  #############

			Namev <<- tclVar(textvariables[indicev[1]])
			tkconfigure(entry.Namev,textvariable=Namev)
    
		}
	
		OK.butv <-tkbutton(frameok,text="    OK    ",command=OnOKv)
		tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")


		tkpack(tlv,scrv,expand = "TRUE", side="left", fill = "both")
		tkpack.configure(scrv,side="left")

		tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")
		tkfocus(tt)



#######Color#######################################
	
		indicev<-as.numeric(tkcurselection(tlv))+1
		colorv <- colvariables[indicev[1]]
		canvasv <- tkcanvas(framecol,width="57",height="20",bg=colorv)

		ChangeColorv <- function()
		{
 
			colorv <<- tclvalue(tcl("tk_chooseColor",initialcolor=colvariables[indicev[1]],title=gettext("Choose a color",domain="R-multibiplotGUI")))
 
		 	if (nchar(colorv)>0)
    				{
					tkconfigure(canvasv,bg=colorv)
		 			colvariables[indicev]<<-colorv
				}
		}

		ChangeColor.buttonv<- tkbutton(framecol,text=gettext("Change Color",domain="R-multibiplotGUI"),command=ChangeColorv,width=2)
		tkpack(canvasv,ChangeColor.buttonv,expand = "TRUE", side="left", fill = "both")
 


######Labels  ###################################
	
		Namev <- textvariables[indicev[1]]
		entry.Namev <-tkentry(framename,width="10",textvariable=Namev)

		OnOKlv <- function()
		{
			NameValv <<- tclvalue(Namev)
			textvariables[indicev]<<-NameValv

#####Values of listbox###############################
	
			for (i in 1:(dim(Xpon)[2]))
			{
				tkdelete(tlv,0)
			}

			for (i in 1:(dim(Xpon)[2]))
			{
   				tkinsert(tlv,"end",textvariables[i])
			}


		}

		OK.butlv <-tkbutton(framename,text=gettext(" Change label",domain="R-multibiplotGUI"),command=OnOKlv,width=2)
		tkbind(entry.Namev, "<Return>",OnOKlv)
		tkpack(entry.Namev,OK.butlv,expand = "TRUE", side="left", fill = "both")

		tkpack(tklabel(frames,text="      ",width=28),expand = "TRUE", side="left", fill = "both")

		tkpack(tklabel(framecol,text=gettext("Show axes",domain="R-multibiplotGUI")),cb,expand = "TRUE", side="left",expand="TRUE", fill = "both")




		Graphics <- function()
		{
			
			Myhscale <<- 2    # Horizontal scaling
			Myvscale <<- 2    # Vertical scaling

	
##############################################################################
#####		Coordinates
##############################################################################
			ejes<<-c(gettext("Axis 1",domain="R-multibiplotGUI"),gettext("Axis 2",domain="R-multibiplotGUI"),gettext("Axis 3", domain="R-multibiplotGUI"))
			if (tipo == "RMP"){

				descom<<-La.svd(Xstd,nu=3,nv=3)

				coindividuos<<-descom$u%*%diag(descom$d[1:3])
				covariables<<-t(descom$v)

				

##############################################################################
#####		Contributions, goodness of fit and qualities of representation
##############################################################################
	
				suma2valprop<<-sum((descom$d[1:3])^2)
				sumaRvalprop<<-sum((descom$d)^2)

				suma2valprop4<<-sum((descom$d[1:3])^4)
				sumaRvalprop4<<-sum((descom$d)^4)
					

				bonajuste<<-(suma2valprop/sumaRvalprop)*100
					
				calcol<<-2/(dim(diag(descom$d))[1])*100

			
			} else {

				descom<<-La.svd(Xstd,nu=3,nv=3)
				
				coindividuos<<-descom$u%*%diag(descom$d[1:3])
				covariables<<-t(descom$v)%*%diag(descom$d[1:3])

##############################################################################
#####		Contributions, goodness of fit and qualities of representation
##############################################################################
	
				suma2valprop<<-sum((descom$d[1:3])^2)
				sumaRvalprop<<-sum((descom$d)^2)

				suma2valprop4<<-sum((descom$d[1:3])^4)
				sumaRvalprop4<<-sum((descom$d)^4)
					
					
				calcol<<-(suma2valprop4/sumaRvalprop4)*100

			}
	
				
	
				coindividuosnam<<-as.data.frame(coindividuos)
				rownames(coindividuosnam)<<-textindividuos
				colnames(coindividuosnam)<<-ejes

				covariablesnam<<-as.data.frame(covariables)
				rownames(covariablesnam)<<-textvariables
				colnames(covariablesnam)<<-ejes

				coindivcuad<<-coindividuos^2
				CRTi<<-rowSums(coindivcuad)
				CRTi<<-(CRTi*1000)/suma2valprop
				CRTi<<-as.data.frame(CRTi)
				rownames(CRTi)<<-textindividuos

				covarcuad<<-covariables^2
				CRTj<<-rowSums(covarcuad)
				CRTj<<-(CRTj*1000)/suma2valprop
				CRTj<<-as.data.frame(CRTj)
				rownames(CRTj)<<-textvariables


				CREiFq<<-array(dim=dim(coindividuos))
				CREiFq[,1]<<-((coindivcuad)[,1]*1000)/((descom$d[1])^2)
				CREiFq[,2]<<-((coindivcuad)[,2]*1000)/((descom$d[2])^2)
				CREiFq[,3]<<-((coindivcuad)[,3]*1000)/((descom$d[3])^2)
				CREiFq<<-as.data.frame(CREiFq)
				rownames(CREiFq)<<-textindividuos
				colnames(CREiFq)<<-ejes
					
				CREjFq<<-array(dim=dim(covariables))
				sumavar<<-rowSums(covarcuad)
				CREjFq[,1]<<-((covarcuad)[,1]*1000)/(sumavar)
				CREjFq[,2]<<-((covarcuad)[,2]*1000)/(sumavar)
				CREjFq[,3]<<-((covarcuad)[,3]*1000)/(sumavar)
				CREjFq<<-as.data.frame(CREjFq)
				rownames(CREjFq)<<-textvariables
				colnames(CREjFq)<<-ejes


				CRFqEi<<-coindivcuad
				sumaindi<<-rowSums(coindivcuad)
				CRFqEi[,1]<<-((coindivcuad)[,1]*1000)/(sumaindi)
				CRFqEi[,2]<<-((coindivcuad)[,2]*1000)/(sumaindi)
				CRFqEi[,3]<<-((coindivcuad)[,3]*1000)/(sumaindi)
				CRFqEi<<-as.data.frame(CRFqEi)
				rownames(CRFqEi)<<-textindividuos
				colnames(CRFqEi)<<-ejes


				CRFqEj<<-covarcuad
				sumavar<<-rowSums(covarcuad)
				CRFqEj[,1]<<-((covarcuad)[,1]*1000)/(sumavar)
				CRFqEj[,2]<<-((covarcuad)[,2]*1000)/(sumavar)
				CRFqEj[,3]<<-((covarcuad)[,3]*1000)/(sumavar)
				CRFqEj<<-as.data.frame(CRFqEj)
				rownames(CRFqEj)<<-textvariables
				colnames(CRFqEj)<<-ejes

				calfilas<<-(suma2valprop4/sumaRvalprop4)*100


				rangos<<-array(dim=c(NameVal,2))

				for(z in 1:NameVal)
				{
					if (z==1){
						rangos[z,1]<<-1
						rangos[z,2]<<-get(dimensiones[z])
					}else{
						rangos[z,1]<<-rangos[z-1,2]+1
						rangos[z,2]<<-rangos[z-1,2]+get(dimensiones[z])
					}
					
				}		


				sumacoordindgrupo<<-array(dim=c(NameVal,3))
				sumacoordvargrupo<<-array(dim=c(NameVal,3))
				
				CRTt<<-c(1:NameVal)
				CRGtFq<<-array(dim=c(NameVal,3))
				CRFqGt<<-array(dim=c(NameVal,3))

				Xaproxestim<<-descom$u%*%diag(descom$d[1:3])%*%descom$v


				for(z in 1:NameVal){
					if(z==1){
						lambda<<-rep(get(valorespropios[z]),get(dimensiones[z]))
					}else{
						lambda<<-c(lambda,rep(get(valorespropios[z]),get(dimensiones[z])))
					}
				}

				if ((filas==1)&&(columnas==0)){
					
					for(z in 1:NameVal)
					{
						sumacoordindgrupo[z,1]<<-sum(coindivcuad[rangos[z,1]:rangos[z,2],1])
						sumacoordindgrupo[z,2]<<-sum(coindivcuad[rangos[z,1]:rangos[z,2],2])
						sumacoordindgrupo[z,3]<<-sum(coindivcuad[rangos[z,1]:rangos[z,2],3])
						CRTt[z]<<-((sumacoordindgrupo[z,1]+sumacoordindgrupo[z,2]+sumacoordindgrupo[z,3])*1000)/suma2valprop
						
						CRGtFq[z,1]<<-(sumacoordindgrupo[z,1]*1000)/((descom$d[1])^2)
						CRGtFq[z,2]<<-(sumacoordindgrupo[z,2]*1000)/((descom$d[2])^2)
            CRGtFq[z,3]<<-(sumacoordindgrupo[z,3]*1000)/((descom$d[3])^2)

						CRFqGt[z,1]<<-(sumacoordindgrupo[z,1]*1000)/(sumacoordindgrupo[z,1]+sumacoordindgrupo[z,2]+sumacoordindgrupo[z,3])
						CRFqGt[z,2]<<-(sumacoordindgrupo[z,2]*1000)/(sumacoordindgrupo[z,1]+sumacoordindgrupo[z,2]+sumacoordindgrupo[z,3])
            CRFqGt[z,3]<<-(sumacoordindgrupo[z,3]*1000)/(sumacoordindgrupo[z,1]+sumacoordindgrupo[z,2]+sumacoordindgrupo[z,3])

					}
					
										
					Nlambda<<-diag(lambda)
					invN<<-solve(Nlambda)
					Xestim<<-invN%*%Xaproxestim
				
				}else{
		
					for(z in 1:NameVal)
					{
						sumacoordvargrupo[z,1]<<-sum(covarcuad[rangos[z,1]:rangos[z,2],1])
						sumacoordvargrupo[z,2]<<-sum(covarcuad[rangos[z,1]:rangos[z,2],2])
						sumacoordvargrupo[z,3]<<-sum(covarcuad[rangos[z,1]:rangos[z,2],3])
						CRTt[z]<<-((sumacoordvargrupo[z,1]+sumacoordvargrupo[z,2]+sumacoordvargrupo[z,3])*1000)/suma2valprop
						
						CRGtFq[z,1]<<-(sumacoordvargrupo[z,1]*1000)/((descom$d[1])^2)
						CRGtFq[z,2]<<-(sumacoordvargrupo[z,2]*1000)/((descom$d[2])^2)
            CRGtFq[z,3]<<-(sumacoordvargrupo[z,3]*1000)/((descom$d[3])^2)


						CRFqGt[z,1]<<-(sumacoordvargrupo[z,1]*1000)/(sumacoordvargrupo[z,1]+sumacoordvargrupo[z,2]+sumacoordvargrupo[z,3])
						CRFqGt[z,2]<<-(sumacoordvargrupo[z,2]*1000)/(sumacoordvargrupo[z,1]+sumacoordvargrupo[z,2]+sumacoordvargrupo[z,3])
            CRFqGt[z,3]<<-(sumacoordvargrupo[z,3]*1000)/(sumacoordvargrupo[z,1]+sumacoordvargrupo[z,2]+sumacoordvargrupo[z,3])

					}
					
										
					Mlambda<<-diag(lambda)
					invM<<-solve(Mlambda)					
					Xestim<<-Xaproxestim%*%invM


				}

				colnames(CRGtFq)<<-ejes
				colnames(CRFqGt)<<-ejes

				
				SCT<<-sum(diag(t(X)%*%X))
				SCEx<<-sum(diag(t(Xestim)%*%Xestim))

				bonajusXstd<<-(SCEx/SCT)*100








			cat(gettext("File saved in:    ",domain="R-multibiplotGUI"),file=gettext("Results.xls",domain="R-multibiplotGUI"))
			cat(getwd(),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")		
			cat(gettext("CONTRIBUTIONS:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution to total variability of the row element i:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRTi,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
				
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution to total variability of the column element j:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRTj,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")

			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the row element i to the factor q-th:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CREiFq,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")

			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the column element j to the factor q-th:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CREjFq,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")

					
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the factor q-th to row element i:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRFqEi,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")


			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the factor q-th to column element j:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRFqEj,file="temp.xls", sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")



			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the set t to total variability: \n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRTt, file="temp.xls", sep="\t", dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")


			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the set t to the factor q: \n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRGtFq, file="temp.xls", sep="\t", dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")


			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Relative contribution of the factor q to the set t: \n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(CRFqGt, file="temp.xls", sep="\t", dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")



			if (tipo == "RMP"){
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Goodnes of Fit to X approximated:  ",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			cat(bonajuste,file="temp.xls")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			cat(" %",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	

			}


			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Goodness of Fit to X:  ",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			cat(bonajusXstd,file="temp.xls")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			cat(" %",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	



			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Quality of approximation for rows:  ",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			cat(calfilas,file="temp.xls")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			cat(" %",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	



			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Quality of approximation for columns:  ",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			cat(calcol,file="temp.xls")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			cat(" %",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	



			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Individual coordinates:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(coindividuosnam,file="temp.xls",sep="\t",dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			



			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Variables coordinates:\n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(covariablesnam, file="temp.xls", sep="\t", dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")
			


			cat("\n",file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")	
			cat(gettext("Eigen values: \n",domain="R-multibiplotGUI"),file="temp.xls")					
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")					
			write.table(descom$d, file="temp.xls", sep="\t", dec=",")
			file.append(gettext("Results.xls",domain="R-multibiplotGUI"),"temp.xls")


			


			file.show(gettext("Results.xls",domain="R-multibiplotGUI"))
			file.remove("temp.xls")

			datos<<-rbind(coindividuos,covariables)
			centro<<-c(0,0)
			

################ Show axes or not
	cbVal <<- as.character(tclvalue(cbValue))

		xCoords<<-datos[,1]
		yCoords<<-datos[,2]
		zCoords<<-datos[,3]
		labelsVec <<- c(textindividuos, textvariables)
		colores<<-c(colindividuos,colvariables)
		simbolos<<-c(simindividuos, simvariables)


	indexLabeled<-c(1:length(xCoords))
	indexLabeledaux<-c()
	labeledPoints <- list()

	
	wgr <- tktoplevel()
	tkwm.title(wgr,gettext("Graph",domain="R-multibiplotGUI"))
	parPlotSize <- c()
	usrCoords <- c()

	plotFunctiond <- function()
	{
  		params <- par(bg="white")
 			plot(datos,main= gettext("Graph",domain="R-multibiplotGUI"),type="n",xlab=" ",ylab=" ")
			

				if ((filas == 1) && (columnas == 0)){
					cexin<<-0.75
					cexvar<<-1
				}else{
					cexin<<-1
					cexvar<<-0.75
				}
				
				points(coindividuos[,1],coindividuos[,2],pch=simindividuos,col=colindividuos)

				arrows(centro[1],centro[2],covariables[,1],covariables[,2],col=colvariables,lty="dotted",length=0.05)

				points(centro[1],centro[2],pch=18,col="black")
				cex<<-c(rep(cexin, times=length(colindividuos)),rep(cexvar,times=length(colvariables)))

################ Show axes or not
				cbVal <<- as.character(tclvalue(cbValue))
				if (cbVal=="1"){
	
					abline(h=centro[2],v=centro[1],lty="dotted")
				}

		



  		if (length(indexLabeled)>0)
    		for (i in (1:length(indexLabeled)))
    		{
      		indexClosest <- indexLabeled[i]
      		text(xCoords[indexClosest],yCoords[indexClosest],
           		labels=labelsVec[indexClosest], col= colores[indexClosest],cex= cex[indexClosest])
    		}
  		parPlotSize <<- par("plt")
  		usrCoords   <<- par("usr")
  		par(params)
	}

	img <- tkrplot(wgr,fun=plotFunctiond,hscale=1.5,vscale=1.5)
	tkpack(img, expand="TRUE", fill="both")




	labelClosestPoint <- function(xClick,yClick,imgXcoords,imgYcoords)
	{
  		squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
  		indexClosest <<- which.min(squared.Distance)
 		if (tclvalue(mbval)=="yes"){  
			indexLabeled <<- c(indexLabeled,indexClosest)


###############################
			xCoords[indexClosest]<<-xPlotCoord
			yCoords[indexClosest]<<-yPlotCoord
###############################

 		}

		if(tclvalue(mbval)=="no"){

 		indexLabeledaux<<-c()
 		for (i in (1:length(indexLabeled)))
    		{
     			if (indexLabeled[i]!=indexClosest)
				indexLabeledaux <<- c(indexLabeledaux,indexLabeled[i])
    		}
		 indexLabeled<<-indexLabeledaux 
		}
  
	tkrreplot(img)
	}




labelClosestPointd <- function(xClick,yClick,imgXcoords,imgYcoords)
	{
  		squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
  		indexClosest <<- which.min(squared.Distance)
		mm<-tktoplevel() 	
		tkwm.title(mm, labelsVec[indexClosest])	


		framemm1<-tkframe(mm, relief = "groove", borderwidth = 2, 
        	background = "white")

		framemm2<-tkframe(mm, relief = "groove", borderwidth = 2, 
        	background = "white")

		framemm3<-tkframe(mm, relief = "groove", borderwidth = 2, 
        	background = "white")



		colori <- colores[indexClosest]
		canvasi <- tkcanvas(framemm1,width="120",height="20",bg=colori)

		ChangeColori <- function()
		{
 
			colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colores[indexClosest],title=gettext("Choose a color",domain="R-multibiplotGUI")))
 
			 if (nchar(colori)>0)
    				{
					tkconfigure(canvasi,bg=colori)
		 			colores[indexClosest]<<-colori
					colindividuos<<-colores[1:length(colindividuos)]
					colvariables<<-colores[(length(colindividuos)+1):(length(colindividuos)+length(colvariables))]
					

				}
			tkrreplot(img)
			tkdestroy(mm)
		}

		ChangeColor.buttoni<- tkbutton(framemm1,text=gettext("Change Color",domain="R-multibiplotGUI"),command=ChangeColori)
		tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")
		tclvalue(Namei) <<- labelsVec[indexClosest]
		entry.Namei <<-tkentry(framemm2,width="10",textvariable=Namei)
		NameVali <<- Namei 

		OnOKli <- function()
		{
			NameVali <<- tclvalue(Namei)
			labelsVec[indexClosest]<<-NameVali

			textindividuos<<-labelsVec[1:length(textindividuos)]
			textvariables<<-labelsVec[(length(textindividuos)+1):(length(textindividuos)+length(textvariables))]
		


#####Values of listbox###############################
		
			for (i in 1:dim(Xpon)[1])
			{
				tkdelete(tli,0)
			}

			for (i in 1:(dim(Xpon)[1]))
			{
   				tkinsert(tli,"end",textindividuos[i])
			}
	
			

			for (i in 1:(dim(Xpon)[2]))
			{
				tkdelete(tlv,0)
			}

			for (i in 1:(dim(Xpon)[2]))
			{
   				tkinsert(tlv,"end",textvariables[i])
			}



			tkrreplot(img)
			tkdestroy(mm)

	}

	OK.butli <-tkbutton(framemm2,text=gettext(" Change label",domain="R-multibiplotGUI"),command=OnOKli,width=2)
	tkbind(entry.Namei, "<Return>",OnOKli)
	tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
	


	comboBox <- tkwidget(framemm3,"ComboBox",editable=FALSE,values=symbols,width=10, text= simbolos[indexClosest])

	chang.sym <- function()
	{
   	 	simChoice <<-symbols[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
		simbolos[indexClosest]<<-simChoice
	
		simindividuos<<-simbolos[1:length(simindividuos)]
		simvariables<<-simbolos[(length(simindividuos)+1):(length(simindividuos)+length(simvariables))]
		

		tkrreplot(img)
		tkdestroy(mm)


 	}
	Change.symbol <-tkbutton(framemm3,text=gettext("   Change symbol   ",domain="R-multibiplotGUI"),command=chang.sym,width=6)
	tkpack(comboBox,Change.symbol,side="left",expand="TRUE", fill="both")


	tkpack(framemm1,framemm2,framemm3,expand = "TRUE", side="top", fill = "both")


	}



	CopyToClip <- function()
	{
  		tkrreplot(img)
	}

	copy.but <- tkbutton(wgr,text=gettext("Export",domain="R-multibiplotGUI"),command=CopyToClip)
	tkpack(img, expand="TRUE", fill="both")
#	tkpack(copy.but,expand="TRUE", fill="both")

  g3d<-function()
  {
      bg3d("white")
      if (cbVal=="1"){
          axes3d()
          }
      points3d(xCoords,yCoords,zCoords, color=colores)
      texts3d(xCoords, yCoords, zCoords,labelsVec,color=colores)
  	
      
	     for (i in 1:(dim(covariables)[1]))
	     {
	       linea<-rbind(covariables[i,],c(0,0,0))	
	       segments3d(linea[,1],linea[,2], linea[,3],color=colvariables[i])

	     }
  }
  
  but.3d <- tkbutton(wgr,text=gettext("3D",domain="R-multibiplotGUI"),command=g3d)
	tkpack(but.3d, copy.but,side="left",expand="TRUE", fill="both")



	OnLeftClick <- function(x,y)
	{
  	xClick <<- x
  	yClick <<- y
  	require(tcltk)
  	width  <<- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
  	height <<- as.numeric(tclvalue(tkwinfo("reqheight",img)))

  	xMin <<- parPlotSize[1] * width
  	xMax <<- parPlotSize[2] * width
    yMin <<- parPlotSize[3] * height
 	  yMax <<- parPlotSize[4] * height

 	 rangeX <<- usrCoords[2] - usrCoords[1]
	 rangeY <<- usrCoords[4] - usrCoords[3]

 	 imgXcoords <<- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
 	 imgYcoords <<- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin

	 xClick <<- as.numeric(xClick)+0.5
 	 yClick <<- as.numeric(yClick)+0.5
 	 yClick <<- height - yClick

	  xPlotCoord <<- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
 	  yPlotCoord <<- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)

	  msg <- (gettext("-To change the label press Yes.\n-To remove it press No.\n-If you do not want to do anything press Cancel.",domain="R-multibiplotGUI"))
	  mbval<<- tkmessageBox(title=gettext("Change of label",domain="R-multibiplotGUI"),
                       message=msg,type="yesnocancel",icon="question")

 	 labelClosestPoint(xClick,yClick,imgXcoords,imgYcoords)

	}




OnRightClick <- function(x,y)
	{
  	xClick <<- x
  	yClick <<- y
  	require(tcltk)
  	width  <<- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
  	height <<- as.numeric(tclvalue(tkwinfo("reqheight",img)))

  	xMin <<- parPlotSize[1] * width
  	xMax <<- parPlotSize[2] * width
	  yMin <<- parPlotSize[3] * height
 	  yMax <<- parPlotSize[4] * height

 	 rangeX <<- usrCoords[2] - usrCoords[1]
	 rangeY <<- usrCoords[4] - usrCoords[3]

 	 imgXcoords <<- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
 	 imgYcoords <<- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin

	 xClick <<- as.numeric(xClick)+0.5
 	 yClick <<- as.numeric(yClick)+0.5
 	 yClick <<- height - yClick

	  xPlotCoord <<- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
 	  yPlotCoord <<- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)

	
 	 labelClosestPointd(xClick,yClick,imgXcoords,imgYcoords)

	}


	tkbind(img, "<ButtonRelease-1>",OnLeftClick)
	tkconfigure(img,cursor="pencil")
	
	tkbind(img, "<Button-3>",OnRightClick)
	tkconfigure(img,cursor="pencil")

	



		}
		graphic.button <-tkbutton(framegraphic,text=gettext("    Graph    ",domain="R-multibiplotGUI"),command=Graphics)



		tkpack(graphic.button, expand="TRUE", side= "left", fill ="both")





	tkpack(tklabel(frametext,text=gettext("Transformations",domain="R-multibiplotGUI")),side="right",expand = "TRUE",fill="both")
	tkpack(tltrans,scrtrans,expand = "TRUE", side="left", fill = "both")
	tkpack.configure(scrtrans,side="left")
	tkpack(tklabel(framename,text="      ",width=25),expand = "TRUE", side="left", fill = "both")
	tkpack(tklabel(frames,text="      ",width=28),expand = "TRUE", side="left", fill = "both")




	tkpack(frametext,framet,frameok,framecol,framename,frames,framegraphic,expand = "TRUE",fill="both")



	}

	OK.but <-tkbutton(wnummat,text="   OK   ",command=OnOK)


#####Dropdown menu#############################

	topMenu <- tkmenu(wnummat)
	tkconfigure(wnummat,menu=topMenu)
	fileMenu <- tkmenu(topMenu,tearoff=FALSE)
	tkadd(fileMenu,"command",label=gettext("Some sets of individuals which have been observed in a single set of variables",domain="R-multibiplotGUI"),
  command=function() {filas<<-1; columnas<<-0})
	tkadd(fileMenu,"command",label=gettext("Some sets of variables observed in a single set of individuals",domain="R-multibiplotGUI"),
  command=function() {columnas<<-1; filas<<-0})
	tkadd(topMenu,"cascade",label=gettext("Analysis",domain="R-multibiplotGUI"),menu=fileMenu)




	tkwm.title(wnummat,gettext("Number of matrices",domain="R-multibiplotGUI"))
	nummat <- tclVar( 1 )
	entry.Name <-tkentry(wnummat,width="50",textvariable=nummat)
	tkbind(entry.Name, "<Return>",OnOK)

	tkgrid(tklabel(wnummat,text=gettext("Number of matrices to analyze:",domain="R-multibiplotGUI")))
	tkgrid(entry.Name)
	tkgrid(OK.but)
	tkfocus(wnummat)


}


	
OK.butinf <-tkbutton(winfor,text="   OK   ",command=OnOKinf)
fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
fontFixedWidth <- tkfont.create(family="courier",size=12)
tkgrid(tklabel(winfor,text="               MULTIBIPLOT               ",font=fontHeading))
tkgrid(tklabel(winfor,text="    "))
tkgrid(OK.butinf)
tkfocus(winfor)

}

