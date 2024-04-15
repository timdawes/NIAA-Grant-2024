# Post-operative LCOSpo fromin neonates
# Tim Dawes
# February 2024

# Code to simulate sample sizes:

library(colorRamps)
library(readxl)
library(Hmisc)
library(RColorBrewer)
library(ggsci)
library(viridis)

cols<- viridis(10, option = "turbo")
# Choose some colours
      #cols<- c(brewer.pal(3,"Reds")[-1], "black", rev(brewer.pal(3, "Blues")[-1]))
      #cols<- brewer.pal(n=8, name='Dark2')
      cols<- pal_npg("nrc")(5)
      

# Functions
      logistic <- function(x) {1 / (1 + exp(-x))}

      log_func <- function(N, organ, weight.log.mean, weight.log.sd, intersect.mat, grad.mat) {
        
        log.weight <- rnorm(N, weight.log.mean, weight.log.sd)
        intersect<- rnorm(N, intersect.mat[organ,1], intersect.mat[organ,2])
        grad<- rnorm(N, grad.mat[organ,1], grad.mat[organ,2])
        
        complication<- round(logistic(intersect + log.weight * grad))
        
        df<- data.frame(log.weight=log.weight, complication=complication)
        fit<- glm(complication ~ log.weight, family=binomial(link="logit"), data=df)
        
        fit2<- summary(fit)
        p<- fit2$coefficients[2,4]
        stat <- fit2$coefficients[2,3]
        
        return(c(t=stat, p=p, sig=(p < .05)))
      }
      
      
      log2_func <- function(N, organ, weight.log.mean, weight.log.sd, intersect.mat, grad.mat) {
        
        log.weight <- rnorm(N, weight.log.mean, weight.log.sd)
        intersect<- rnorm(N, intersect.mat[organ,1], intersect.mat[organ,2])
        grad<- rnorm(N, grad.mat[organ,1], grad.mat[organ,2])
        
        complication<- round(logistic(intersect + log.weight * grad))
        
        df<- data.frame(log.weight=log.weight, complication=complication)
        fit<- glm(complication ~ log.weight, family=binomial(link="logit"), data=df)
        
        fit2<- summary(fit)
        p<- fit2$coefficients[2,4]
        stat <- fit2$coefficients[2,3]
        
        return(c(t=stat, p=p, sig=(p < .05)))
      }


# Read in morbidity data
    d<- data.matrix(read.csv("morbidities.csv"))

        
fit<- list()
par(mfrow=c(2,3))
morbidities<- c("Vent","AKI","NEC","Lactate")
intersect.mat<- grad.mat<- matrix(0, nrow=length(morbidities), 2, dimnames=list(morbidities, c("Mean","SD")))
col.match<- match(morbidities, colnames(d))

# Extract the data which you need for the power calculation
      for (organ in 1:length(morbidities))
      {
        value.rows<- which(is.na(d[,col.match[organ]])==FALSE)
        y<- d[value.rows,col.match[organ]]
        if (identical(sort(unique(y)), c(1,2))) {y<- y - 1}
        
        x<- d[value.rows,1]
        boxplot(log(x) ~ y, xlab="", ylab="", main=colnames(d)[col.match[organ]])
        mtext(side=1, "Complication", line=2, cex=1)
        mtext(side=2, "Log Weight (kg)", line=2, cex=1)
        
        fit[[organ]]<- glm(y ~ log(x), family=binomial(link="logit"))
        
        coefficients.GOSH<- summary(fit[[organ]])$coefficients
        intersect.mat[organ,]<- c(coefficients.GOSH[1,1], coefficients.GOSH[1,2]*1.96)
        grad.mat[organ,]<- c(coefficients.GOSH[2,1], coefficients.GOSH[2,2]*1.96)
            
      }

        weight.log.mean<- mean(na.omit(log(d[,1])))
        weight.log.sd<- sd(na.omit(log(d[,1])))

# Loops to check whether effect is detected in successive trials
      df<- list()
      par(mfrow=c(1,1), mar=c(6,6,2,3), mgp=c(1,1,1))
      calc<- TRUE
      trialSize<- c(10,20,30,40,50,75,100,150,200,300,400,500,600,700,800,900,1000)
      spans<- c(0.2,0.3,0.6,0.8)
      maxTrials<- 1000
      sample.size<- rep(0, 4)
      options(warn=-1)
      
      for (k in 1:length(morbidities))
      {
        cat(k)
        sig<- rep(0, length(trialSize))
        
      if (calc==TRUE){
            for (i in trialSize)
            {
              cat(".")
              for (j in 1:maxTrials)
              {
              corr.results<- log_func(i, k, weight.log.mean, weight.log.sd, intersect.mat, grad.mat)
              sig[match(i, trialSize)]<- sig[match(i, trialSize)] + corr.results[3] 
              }
              df[[k]]<- data.frame(X=trialSize, Y=100 * sig / maxTrials)
            }
      }
        
        if (k==1) {plot(trialSize, df[[k]]$Y, type='p', lwd=1, col=cols[k], pch=19, cex=1, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', xlim=c(0,max(trialSize)), ylim=c(0,100))}
        if (k>1) {for (l in 2:k) {points(df[[l]]$X, df[[l]]$Y, type='p', lwd=1, col=cols[l], pch=19, cex=1, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', xlim=c(0,max(trialSize)), ylim=c(0,100))}}
        
        axis(side=1, lwd=5, at=round(seq(0,max(trialSize),length.out=5)), line=0.5, cex.axis=2, font=2)
        mtext("Number of patients in study cohort", side=1, line=3.5, cex=2, font=2)
        axis(side=2, lwd=5, at=seq(0,100,25), line=-1, cex.axis=2, font=2, las=2)
        mtext("Power (%)", side=2, line=3.5, cex=2, font=2)
        lo<- loess(Y~X, span=spans[k], data=df[[k]])
        sample.size[k]<- df[[k]]$X[which.min(abs(predict(lo,df[[k]]$X)-90))]
        
        new.y<- with(df[[k]], predict(lo, X))
        new.y[new.y>100]<- 100
        new.x<- df[[k]]$X
        points(new.x, new.y, type='l', lwd=7, col=cols[k])
      } 
      
# Sample size needed in each group for 90% power
      
      abline(h=90, lwd=5, col="black", lty=2)
      abline(v=sample.size[2], lwd=5, col="black", lty=2)
      
      text(850, 85, "90% Power", cex=1.8, font=2)
      step<- 7
      text(800, (5*step)-2, "Outcome", cex=1.8, font=2)
      
      for (k in 1:length(morbidities))
          {
          segments(600,(k*step)-5,1000,(k*step)-5, lwd=20, col=cols[k])
          t<- morbidities[k]
          text(800,(k*step)-5, labels=t, col="white", cex=1.5)
      }
      
      