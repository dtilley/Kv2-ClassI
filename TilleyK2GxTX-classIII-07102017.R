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

TilleyK2GxTX.KKABKK.072017 <- function(V,x,T,dt,allstates=FALSE,returnNDX=6,gating=FALSE){

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
    # The KKABKK model utilizes bsa rates
    bsaAlpha <- function(v){
        alphas <- exp(1.439272-0.025609*v)
        alphas <- alphas^-1
        return(alphas)
    }
    
    # Determine beta rate vector
    # The KKABKK model utilizes bsa rates
    bsaBeta <- function(v){
        betas <- exp(2.295625-0.006848*v)
        betas <- betas^-1
        return(betas)
    }

    # Determine ktrap rate vector
    gxtxKtrap <- function(v){
        ktrap <- exp(3.292799-0.036706*v)
        ktrap <- ktrap^-1
        return(ktrap)
    }

    # Determine ktrap rate vector
    gxtxKrelease <- function(v){
        krelease <- exp(4.746323-0.040223*v)
        krelease <- krelease^-1
        return(krelease)
    }

    # Initialize alpha & beta
    alphas <- bsaAlpha(V)
    betas <- bsaBeta(V)
    ktrap <- gxtxKtrap(V)
    krelease <- gxtxKrelease(V)
    
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
        Alist[[i]] <- initMatrix(alphas[i],betas[i],kopen,kclose,ktrap[i],krelease[i])
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
    x <- matrix(ncol=1,nrow=16,data=0)
    x[ndx,1] <- 1
    return(x)
}

# Initialize Transfer Matrix for Triangle Model Gxtx Interaction
initMatrix <- function(a,b,c,d,kon,koff){
    m <- matrix(ncol=16,nrow=16,data=0)
    # Populate transitions R4
    m[1,1] <- -(4*a+4*kon)
    m[1,2] <- b
    m[1,7] <- koff
    # Populate transitions R3
    m[2,2] <- -(3*a + 3*kon + b)
    m[2,1] <- 4*a
    m[2,3] <- 2*b
    m[2,8] <- koff
    # Populate transitions R2
    m[3,3] <- -(2*a + 2*kon + 2*b)
    m[3,2] <- 3*a
    m[3,4] <- 3*b
    m[3,9] <- koff
    # Populate transitions R1
    m[4,4] <- -(a + kon + 3*b)
    m[4,3] <- 2*a
    m[4,5] <- 4*b
    m[4,10] <- koff
    # Populate transitions A
    m[5,5] <- -(c + 4*b)
    m[5,4] <- a
    m[5,6] <- d
    # Populate transitions 0
    m[6,6] <- -(d)
    m[6,5] <- c
    # Populate transitions R4tx
    m[7,7] <- -(koff + 3*a + 3*kon)
    m[7,8] <- b
    m[7,11] <- 2*koff
    m[7,1] <- 4*kon
    # Populate transitions R3tx
    m[8,8] <- -(koff + 2*a + b + 2*kon)
    m[8,7] <- 3*a
    m[8,9] <- 2*b
    m[8,12] <- 2*koff
    m[8,2] <- 3*kon
    # Populate transitions R2tx
    m[9,9] <- -(koff + a + kon + 2*b)
    m[9,8] <- 2*a
    m[9,10] <- 3*b
    m[9,13] <- 2*koff
    m[9,3] <- 2*kon
    # Populate transitions R1tx
    m[10,10] <- -(koff + 3*b)
    m[10,9] <- a
    m[10,4] <- kon
    # Populate transitions R4tx2
    m[11,11] <- -(2*koff + 2*a + 2*kon)
    m[11,12] <- b
    m[11,14] <- 3*koff
    m[11,7] <- 3*kon
    # Populate transitions R3tx2
    m[12,12] <- -(2*koff + a + kon + b)
    m[12,13] <- 2*b
    m[12,11] <- 2*a
    m[12,15] <- 3*koff
    m[12,8] <- 2*kon
    # Populate transitions R2tx2
    m[13,13] <- -(2*koff + 2*b)
    m[13,12] <- a
    m[13,9] <- kon
    # Populate transitions R4tx3
    m[14,14] <- -(3*koff + a + kon)
    m[14,15] <- b
    m[14,11] <- 2*kon
    m[14,16] <- 4*koff
    # Populate transitions R3tx3
    m[15,15] <- -(3*koff + b)
    m[15,14] <- a
    m[15,12] <- kon
    # Populate transitions R4tx4
    m[16,16] <- -(4*koff)
    m[16,14] <- kon

    return(m)
}

# Initialize Gating Transfer Matrix for Triangle Model Gxtx Interaction
initGatrix2 <- function(a,b){
    m <- matrix(ncol=16,nrow=16,data=0)
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

    # Populate transitions 0

    # Populate transitions R4tx
    m[7,7] <- (3*a)

    # Populate transitions R3tx
    m[8,8] <- (2*a - b)

    # Populate transitions R2tx
    m[9,9] <- (a - 2*b)

    # Populate transitions R1tx
    m[10,10] <- -(3*b)

    # Populate transitions R4tx2
    m[11,11] <- (2*a)

    # Populate transitions R3tx2
    m[12,12] <- (a - b)

    # Populate transitions R2tx2
    m[13,13] <- -(2*b)

    # Populate transitions R4tx3
    m[14,14] <- a

    # Populate transitions R3tx3
    m[15,15] <- -(b)

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

runIV <- function(Vstep,Vhold=-100,tstep=c(50,100,200),prefix="gxtxk2IV.KKABKK-"){
    x <- initX(16)
    for (i in seq(1,length(Vstep))){
        V <- c(Vhold,Vstep[i],Vhold)
        O <- TilleyK2GxTX.KKABKK.072017(V,x,T=tstep,dt=0.05,allstates=FALSE,gating=TRUE)
        name <- paste(prefix,i,".txt",sep="")
        write.table(x=O,file=name,row.names=F)
    }
}

runAct <- function(Vstep,Vhold=-100,tstep=c(50,100),prefix="gxtxk2Act-"){
    x <- initX(1)
    for (i in seq(1,length(Vstep))){
        V <- c(Vhold,Vstep[i])
        O <- TilleyK2GxTX.KKABKK.072017(V,x,T=tstep,dt=0.05,allstates=FALSE,gating=TRUE)
        name <- paste(prefix,i,".txt",sep="")
        write.table(x=O,file=name,row.names=F)
    }
}
