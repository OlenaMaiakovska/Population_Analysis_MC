#!/bin/R
#
# AUTHOR : Carine Legrand
#
# Date   : 4.12.2019
#
# Object : Bootstrap and sampling procedure for estimate of population size variability
#

D <- c(11,6,4,1,0,1,1,0,1)
T <- c(11,6,4,1,0,1,1,0,1)
L <- c(2,0,1,0,0,0,0,0,0)
a.m <- 593.0743
a.se <- 21.5518
I.m <- 40.467
I.se <- 4.467

#### Bootstrap, cath efficiency
      T.sur.L.collect <- c()
      set.seed(56178961)
      for (b in 1:10000) {
         ix.b <- sample(1:9, replace=TRUE, size=9)
         T.b <- T[ix.b]
         L.b <- L[ix.b]
         if (sum(L.b) > 0) {
           T.sur.L.b <- sum(T.b) / sum(L.b)
         } else {
           T.sur.L.b <- NA
         }
         T.sur.L.collect <- c(T.sur.L.collect, T.sur.L.b)
      }

      summary(T.sur.L.collect)
      hist(T.sur.L.collect, col="red")

      T.sur.L.median <- median(T.sur.L.collect, na.rm=TRUE)
      T.sur.L.mad    <- mad(T.sur.L.collect, na.rm=TRUE)
      c(T.sur.L.median, T.sur.L.mad)
      # 8.000           2.471
      # center          mad, with constant for equivalence with sd


#### Bootstrap, population size, >4cm or all
      N.ge.4cm.collect <- c()
      N.all.collect <- c()
      set.seed(56178962)
      for (b in 1:10000) {
         # Catch efficiency
         ix.b <- sample(1:9, replace=TRUE, size=9)
         T.b <- T[ix.b]
         L.b <- L[ix.b]
         # Density
         iD.b <- sample(1:9, replace=TRUE, size=9)
         D.m  <- mean(D[iD.b])
         D.se <- sd(D[iD.b]) / sqrt(9)
         # Asymptote
         a.se2 <- sqrt(a.se**2 + (D.se/D.m*a.m)^2)
         A.b <- rnorm(n=1,mean=a.m, sd=a.se2)
         # Influence zone
         I.b <- rnorm(n=1,mean=I.m, sd=I.se)
         # Nb animaux >4cm
         N.ge.4cm.b <- A.b * I.b
         N.ge.4cm.collect <- c(N.ge.4cm.collect, N.ge.4cm.b)
         # Nb
         if (sum(L.b) > 0) {
           N.all.b <- N.ge.4cm.b * (sum(T.b) / sum(L.b))
         } else {
           N.all.b <- NA
         }
         N.all.collect <- c(N.all.collect, N.all.b)
      }
      ####
      summary(N.ge.4cm.collect)
      hist(N.ge.4cm.collect,col='red')
      median(N.ge.4cm.collect)
      mad(N.ge.4cm.collect)
      #
      summary(N.all.collect)
      hist(N.all.collect,col='red')
      median(N.all.collect, na.rm=TRUE)
      mad(N.all.collect, na.rm=TRUE)
      #
      list(N.ge4.median=median(N.ge.4cm.collect),
           N.ge4.mad=mad(N.ge.4cm.collect),
           N.ge4.mad.pc=100*mad(N.ge.4cm.collect)/median(N.ge.4cm.collect),
           N.all.median= median(N.all.collect, na.rm=TRUE),
           N.all.mad=mad(N.all.collect, na.rm=TRUE),
           N.all.mad.pc=100*mad(N.all.collect, na.rm=TRUE)/median(N.all.collect, na.rm=TRUE))

#$N.ge4.median
#[1] 23608.47

#$N.ge4.mad
#[1] 9876.645

#$N.all.median
#[1] 191977

#$N.all.mad
#[1] 103165.8


