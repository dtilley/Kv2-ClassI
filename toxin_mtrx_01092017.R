#####################################################################################
#
#
#  Model Kv2.1 w/Gxtx
#  D.C.Tilley - 08/2016 - drew.tilley@gmail.com
#
#  IN:         N = number of ms to simulate
#             dt = time step
#              x = initial state occupancy
#              A = transfer matrix
#      allstates = if TRUE all the states are returned
#      returnNDX = returns a single output vector, if allstates=T is overwritten
#         gating = gating transfer matrix
#
#
#####################################################################################

# Required library for renaming dataframes #
# library(plyr)

TilleyK2Gxtx <- function(A,x,N,dt,allstates=FALSE,returnNDX=6,gating=NULL,labels=NULL){
    
    # Initialize time, output, gating currents
    t <- (0:(N/dt))*dt
    O <- t(x)
    gatrix <- gating
    if (length(gatrix)==0){
        gating <- FALSE
    } else {
        gating <- TRUE
        ig <- rep(0,length(t))
    }

    
    # loop through time
    for (i in seq(2,length(t))){
        if (gating==TRUE){
            ig[i] <- sum(runge4kuttaIG(A=gatrix,x=x,dt=dt))
        }
        x <- runge4kutta(A=A,x=x,dt=dt)
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
            if (length(labels)!=0){
                #output <- rename(output,labels)
            }
            output <- as.data.frame(cbind(t,output))
            return(output)
        } else {
            output <- as.data.frame(O)
            if (length(labels)!=0){
                #output <- rename(output,labels)
            }
            output <- as.data.frame(cbind(t,ig,output))
            return(output)
        }
    } else {
        if (gating==FALSE){
            output <- as.data.frame(O)
            if (length(labels)!=0){
                #output <- rename(output,labels)
            }
            output <- as.data.frame(cbind(t,O))
            return(output)
        } else {
            output <- as.data.frame(O)
            if (length(labels)!=0){
                #output <- rename(output,labels)
            }
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
initGatrix <- function(a,b,c,d){
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
    m[5,5] <- (c - 4*b)

    # Populate transitions 0
    m[6,6] <- -(d)

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
    # Fit Between 0.1 and 0.9
    fitrange <- c(max(which(ds[[2]]<0.1))+1,max(which(ds[[2]]<0.9))+1)
    N <- length(ds[[2]])
    dstrunc <- as.data.frame(cbind(ds[[1]][fitrange[1]:fitrange[2]],ds[[2]][fitrange[1]:fitrange[2]]))
    mynls <- nlsLM(V2~(1-exp(-V1/tau))^n,data=dstrunc,start=list(tau=tau,n=n),trace=trace)
    ntau <- coef(mynls)
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


