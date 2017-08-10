##########################################################################################
#
#
#  Model Kv2.1 w/Gxtx
#  D.C.Tilley - 07/2017 - drew.tilley@gmail.com
#
#  IN:         T (vector) = time duration for voltage protocol in ms to simulate
#              V (vector) = voltage step protocol
#                      dt = time step
#                       x = initial state occupancy
#                       A = transfer matrix
#               allstates = if TRUE all the states are returned
#               returnNDX = returns a single output vector, if allstates=T is overwritten
#                  gating = T/F if gating currents should be calculated
#
#  This file should be source("./TilleyK2_07092017.R") during the R-session.
##########################################################################################

TilleyK2GxTx.072017 <- function(V,x,T,dt,allstates=FALSE,returnNDX=6,gating=FALSE){

    # Define the voltage protocol Kv2-ABKK model
    if (length(T)!=length(V)){
        print("Error: Voltage and time steps do not match in length.")
        return()
    } else {
        N <- sum(T)
        TVstep <- T[1]
        for (i in seq(2,length(T))){
            TVstep[i] <- T[i]+TVstep[i-1]
        }
    }

    # Define voltage dependent rates
    # Kopen/Kclose is constant determined from fitting macroscopic time constants [citation :Tilley,et al (2017)]
    kopen <- 0.1552064
    kclose <- 0.0388016

    # Determine alpha rate vector
    gxtxAlpha <- function(v){
        alphas <- exp(4.107555-0.028739*v)
        alphas <- alphas^-1
        return(alphas)
    }
    
    # Determine beta rate vector
    gxtxBeta <- function(v){
        betas <- exp(4.584161-0.020015*v)
        betas <- betas^-1
        return(betas)
    }

    # Initialize alpha & beta
    alphas <- gxtxAlpha(V)
    betas <- gxtxBeta(V)
    # Initialize time, output, gating currents
    t <- (0:(N/dt))*dt
    O <- t(x)
    if (gating){
        ig <- rep(0,length(t))
        Glist <- list()
    }

    # Define transfer matrix and ig.matrix
    Alist <- list()

    for (i in seq(1,length(V))){
        Alist[[i]] <- initMatrix(alphas[i],betas[i],kopen,kclose)
        if (gating) {
            Glist[[i]] <- initGatrix2(alphas[i],betas[i])
        }
    }

    # loop through time
    # Voltage Step Counter
    vcounter <- 1
    for (i in seq(2,length(t))){
        if (t[i]<=TVstep[vcounter]) {
            x <- runge4kutta(A=Alist[[vcounter]],x=x,dt=dt)
            if (gating==TRUE){
                ig[i] <- sum(runge4kuttaIG(A=Glist[[vcounter]],x=x,dt=dt))
            }
        } else {
            vcounter <- vcounter+1
            x <- runge4kutta(A=Alist[[vcounter]],x=x,dt=dt)
            if (gating==TRUE){
                ig[i] <- sum(runge4kuttaIG(A=Glist[[vcounter]],x=x,dt=dt))
            }
        }

        sumofStates <- sum(x)
        x <- x*(1/sumofStates)
        error <- abs(1-sumofStates)
        if (error >= 0.01){
            print("Warning: Error exceeded 1%")
            return()
        }
        O <- rbind(O,t(x))
    }

    # format output
    if (allstates==FALSE){
        dims <- dim(O)
        O <- O[1:dims[1],returnNDX]
        if (gating==FALSE){
            output <- as.data.frame(O)
            output <- as.data.frame(cbind(t,output))
            return(output)
        } else {
            output <- as.data.frame(O)
            output <- as.data.frame(cbind(t,ig,output))
            return(output)
        }
    } else {
        if (gating==FALSE){
            output <- as.data.frame(O)
            output <- as.data.frame(cbind(t,O))
            return(output)
        } else {
            output <- as.data.frame(O)
            output <- as.data.frame(cbind(t,ig,O))
            return(output)
        }
    }

}

# 4th Order Rungekutta
runge4kutta <- function(A,x,dt){
    K1 <- dt*(A %*% x)
    K2 <- dt*(A %*% (x+0.5*K1))
    K3 <- dt*(A %*% (x+0.5*K2))
    K4 <- dt*(A %*% (x+K3))
    x1 <- x+(K1+2*K2+2*K3+K4)*(1/6)
    return(x1)
}

runge4kuttaIG <- function(A,x,dt){
    K1 <- dt*(A %*% x)
    K2 <- dt*(A %*% (x+0.5*K1))
    K3 <- dt*(A %*% (x+0.5*K2))
    K4 <- dt*(A %*% (x+K3))
    x1 <- (K1+2*K2+2*K3+K4)*(1/6)
    return(x1)
}

# Initialize State Vector for Triangle Model Gxtx Interaction
initX <- function(ndx){
    x <- matrix(ncol=1,nrow=6,data=0)
    x[ndx,1] <- 1
    return(x)
}

# Initialize Transfer Matrix for ABKK Model
initMatrix <- function(a,b,c,d){
    m <- matrix(ncol=6,nrow=6,data=0)
    # Populate transitions R4
    m[1,1] <- -(4*a)
    m[1,2] <- b
    # Populate transitions R3
    m[2,2] <- -(3*a +  b)
    m[2,1] <- 4*a
    m[2,3] <- 2*b
    # Populate transitions R2
    m[3,3] <- -(2*a +  2*b)
    m[3,2] <- 3*a
    m[3,4] <- 3*b
    # Populate transitions R1
    m[4,4] <- -(a + 3*b)
    m[4,3] <- 2*a
    m[4,5] <- 4*b
    # Populate transitions A
    m[5,5] <- -(c + 4*b)
    m[5,4] <- a
    m[5,6] <- d
    # Populate transitions 0
    m[6,6] <- -(d)
    m[6,5] <- c
    return(m)
}

# Initialize Gating Transfer Matrix for 6 State Model 
# Assumes some gating charge is apparent in Pore Opening
initGatrix <- function(a,b,c,d){
    m <- matrix(ncol=6,nrow=6,data=0)
    # Populate transitions R4
    m[1,1] <- (4*a)

    # Populate transitions R3
    m[2,2] <- (3*a - b)

    # Populate transitions R2
    m[3,3] <- (2*a - 2*b)

    # Populate transitions R1
    m[4,4] <- (a - 3*b)

    # Populate transitions A
    m[5,5] <- (c - 4*b)

    # Populate transitions 0
    m[6,6] <- -(d)

    return(m)
}

# Assumes No charge movement for pore opening
initGatrix2 <- function(a,b){
    m <- matrix(ncol=6,nrow=6,data=0)
    # Populate transitions R4
    m[1,1] <- (4*a)

    # Populate transitions R3
    m[2,2] <- (3*a - b)

    # Populate transitions R2
    m[3,3] <- (2*a - 2*b)

    # Populate transitions R1
    m[4,4] <- (a - 3*b)

    # Populate transitions A
    m[5,5] <- (-4*b)

    return(m)
}

library(minpack.lm)
exppower <- function(ds,n,tau,trace=FALSE,plot=FALSE,coefonly=TRUE,...){
    # Determine if a good fit can be found
    # Color code 1 for R goodfit=TRUE (black)
    # Color code 2 for R goodfit=FALSE (red)
    goodfit <- TRUE
    color <- 1
    N <- length(ds[[2]])
    if(max(which(ds[[2]]<0.9))==N){
        goodfit <- FALSE
        color <- 2
    }
    # Fit Between 0.1 and 0.9
    if (goodfit){
        fitrange <- c(max(which(ds[[2]]<0.1))+1,max(which(ds[[2]]<0.9))+1)
        dstrunc <- as.data.frame(cbind(ds[[1]][fitrange[1]:fitrange[2]],ds[[2]][fitrange[1]:fitrange[2]]))
        mynls <- nlsLM(V2~(1-exp(-V1/tau))^n,data=dstrunc,start=list(tau=tau,n=n),trace=trace)
        ntau <- coef(mynls)
        ntau <- as.vector(ntau)
        ntau <- c(ntau,color)
        if (plot){
            t2 <- c(rep(NA,fitrange[1]-1),dstrunc[[1]],rep(NA,N-fitrange[2]))
            y2 <- c(rep(NA,fitrange[1]-1),predict(mynls),rep(NA,N-fitrange[2]))
            par(mar=c(5.1,5.1,5.1,2.1))
            plot(ds,...)
            lines(t2,y2,col="red",type="l",lwd=2)
            eq <- bquote(y==(1-exp(-t/.(round(ntau[1],4))))^.(round(ntau[2],4)))
            text(locator(1),as.expression(eq))
        }
        if(coefonly){
            return(ntau)
        } else{
            return(mynls)
        }
    } else {
        O <- ds[[3]]
        t <- ds[[1]]
        ds <- as.data.frame(cbind(t,O))
        mynls <- nlsLM(O~(1-exp(-t/tau))^n,data=ds,start=list(tau=tau,n=n),trace=trace)
        ntau <- coef(mynls)
        ntau <- as.vector(ntau)
        ntau <- c(ntau,color)
        if (plot){
            t2 <- ds[[1]]
            y2 <- predict(mynls)
            par(mar=c(5.1,5.1,5.1,2.1))
            plot(ds,...)
            lines(t2,y2,col="red",lwd=2)
            eq <- bquote(y==(1-exp(-t/.(round(ntau[1],4))))^.(round(ntau[2],4)))
            text(locator(1),as.expression(eq))
        }
        if(coefonly){
            return(ntau)
        } else{
            return(mynls)
        }
    }
}

expIG <- function(ds,tau,trace=FALSE,plot=FALSE,coefonly=TRUE,...){
    # Fit Between 0.9 and >0
    fitrange <- c(min(which(ds[[2]]>0.9)),max(which(ds[[2]]>0.1)))
    N <- length(ds[[2]])
    dstrunc <- as.data.frame(cbind(ds[[1]][fitrange[1]:fitrange[2]],ds[[2]][fitrange[1]:fitrange[2]]))
    mynls <- nlsLM(V2~(exp(-V1/tau)),data=dstrunc,start=list(tau=tau),trace=trace)
    tau <- coef(mynls)
    if (plot){
        t2 <- c(rep(NA,fitrange[1]-1),dstrunc[[1]],rep(NA,N-fitrange[2]))
        y2 <- c(rep(NA,fitrange[1]-1),predict(mynls),rep(NA,N-fitrange[2]))
        par(mar=c(5.1,5.1,5.1,2.1))
        plot(ds,...)
        lines(t2,y2,col="red",type="l",lwd=2)
        eq <- bquote(y==(exp(-t/.(round(tau,4)))))
        text(locator(1),as.expression(eq))
    }
    if(coefonly){
        return(tau)
    } else{
        return(mynls)
    }
}


plotLogLinearLabel <- function(x,y,ylabels,xlabels,...){
    # ylabels should be a linear sequence vector
    yticks <- log(ylabels)
    xticks <- xlabels
    ymin <- min(yticks)
    ymax <- max(yticks)
    logY <- log(y)
    xmin <- min(xlabels)
    xmax <- max(xlabels)
    plot(x,logY,axes=F,ylim=c(ymin,ymax),...)
    axis(1,at=xticks,labels=xlabels)
    axis(2,at=yticks,labels=ylabels,las=2)
}

runIV <- function(Vstep,Vhold=-100,tstep=c(50,100,200),prefix="gxtxk2IV-"){
    x <- initX(1)
    for (i in seq(1,length(Vstep))){
        V <- c(Vhold,Vstep[i],Vhold)
        O <- TilleyK2GxTx.072017(V,x,T=tstep,dt=0.1,allstates=FALSE,gating=TRUE)
        name <- paste(prefix,i,".txt",sep="")
        write.table(x=O,file=name,row.names=F)
    }
}

runAct <- function(Vstep,Vhold=-100,tstep=c(50,100),prefix="gxtxk2Act-"){
    x <- initX(1)
    for (i in seq(1,length(Vstep))){
        V <- c(Vhold,Vstep[i])
        O <- TilleyK2GxTx.072017(V,x,T=tstep,dt=0.1,allstates=FALSE,gating=TRUE)
        name <- paste(prefix,i,".txt",sep="")
        write.table(x=O,file=name,row.names=F)
    }
}
