#####################################################################################
#
#
#  Model Kv2.1 Class I Model 
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

SigmaTau2rateMatrix <- function(sigma,tau){
    # Tau should be in passed in ms
    # model function to transform sigma -> alpha/kopen 
    model.AoC<- c(1.28,2.56,5.12,10.24,20.48,40.96,81.92,163.84)
    model.sigmas <- c(5.02802,3.58649,2.28934,1.60055,1.28164,1.13578,1.06646,1.03299)
    sigma2AoC <- splinefun(x=log(model.sigmas),y=log(model.AoC))

    # model function to transform AoC -> rates
    model.AlphaTau <- c(1.47633,2.31074,4.37302,9.04926,18.9965,39.299,80.1666,169.245)
    AoC2AT <- splinefun(x=log(model.AoC),y=log(model.AlphaTau))

    # Calculate rates
    if (sigma < 4 && sigma > 1.03299){
        AoC <- exp(sigma2AoC(log(sigma)))
        AT <- exp(AoC2AT(log(AoC)))
        a <- AT/tau
        c <- (AoC^-1)*a
        rates <- c(a,c)
        return(rates)
    } else {
        print("Sigma Value is out of range.")
        return()
    }

}


TilleyK2Gxtx.classI <- function(A,x,N,dt,allstates=FALSE,returnNDX=6,gating=NULL){
    # Take this out
    labels <- NULL

    
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

# Initialize State Vector for 6 State Gating Model
initX <- function(ndx=1){
    x <- matrix(ncol=1,nrow=6,data=0)
    x[ndx,1] <- 1
    return(x)
}

# Initialize Transfer Matrix for 6 State Gating Model
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
    m[3,3] <- -(2*a + 2*b)
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

# Initialize Gating Transfer Matrix for 6 State Gating Model
initGatrix <- function(a,b){
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
    m[5,5] <- (- 4*b)

    return(m)
}

# Analysis Functions 
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


