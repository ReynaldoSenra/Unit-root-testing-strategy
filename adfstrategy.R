MIT License

Copyright (c) [2017] [Reynaldo Senra]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
  
  The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


# All the ADF tests are at 0.05 percent of significance
status = c("unknown", "nogrow", "grow")
order =5 # check autocorrelation upto 5th order
q = 40 # Upto how many lags you are willing to go to correct autocorrelation
pvalu = 0.05 # significance level of the breush godfrey that you want
adfstrategy<-function(y, status = "unknown", lags = 14, selectlags = "AIC", order =5, order.by = NULL, q = 40, pvalu = 0) {
  alladf<-list()
  totadf<-list()
  for (i in 2:ncol(y)){
    nombre<-colnames(y)[i]
    x<-na.omit(y[,i])
    breush<-c()
    lagss<-lags
    selectlagss<-selectlags
    if (status == "grow") {
      adfresult<-ur.df(x, type = "trend", lags, selectlags)
      z<-diff(x)
      n <- length(z)
      b<-(nrow(adfresult@testreg$coefficients)-2)
      w <- embed(z, b)
      z.diff <- w[, 1]
      z.lag.1 <- x[b:n]
      trend <- b:n
      z.diff.lag = w[, 2:b]
      regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
      for (k in 1:order){
        godfrey<-bgtest(regre, order = k, order.by, type = "Chisq")
        breush[k]<-godfrey$p.value
      }
      if (any(breush <= pvalu)) {
        repeat {
          adfresult<-ur.df(x, type = "trend", lags = b, selectlags = "Fixed")
          b<-b+1
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          trend <- b:n
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
          for (a in 1:order){
            godfrey<-bgtest(regre, order = a, type = "Chisq")
            breush[a]<-godfrey$p.value
          }
          if ((b == (q+1)) | (all(breush > pvalu))){
            if ((any(breush <= pvalu))) {
              warning(paste("Significant autocorrelation problems in the initial ADF test for the series", nombre, "that is growing"))
            }
            break
          }
        }
      }
      purinto<-breush
      purinte<-adfresult
      if (adfresult@teststat[3] <= adfresult@cval[6]) {
        adfresult<-ur.df(x, type = "drift", lags = lagss, selectlags = selectlagss)
        z<-diff(x)
        n <- length(z)
        b<-(nrow(adfresult@testreg$coefficients)-1) 
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
        for (j in 1:order){
          godfrey<-bgtest(regre, order = j, order.by, type = "Chisq")
          breush[j]<-godfrey$p.value
        }
        if (any(breush <= pvalu)) {
          repeat {
            adfresult<-ur.df(x, type = "drift", lags = b, selectlags = "Fixed")
            b<-b+1
            w <- embed(z, b)
            z.diff <- w[, 1]
            z.lag.1 <- x[b:n]
            z.diff.lag = w[, 2:b]
            regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
            for (l in 1:order){
              godfrey<-bgtest(regre, order = l, type = "Chisq")
              breush[l]<-godfrey$p.value
            }
            if ((b == (q+1)) | (all(breush > pvalu))){
              if ((any(breush <= pvalu))) {
                warning(paste("Significant autocorrelation problems the ADF test with drift for the series", nombre, "that is growing"))
              }
              break
            }
          }
        }
        printi<-breush
        summa<-adfresult
        alladf<-list(series=nombre, first="Initial ADF test result for a growing series", firstBreush=purinto, firstADF=summary(purinte), final="The series has a unit root with drift", finalBreush=printi, finalADF=summary(summa))
        print(alladf)
      } else {
        adfresult<-ur.df(x, type = "trend", lags = lagss, selectlags = selectlagss)
        z<-diff(x)
        n <- length(z)
        b<-(nrow(adfresult@testreg$coefficients)-2)
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        trend <- b:n
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
        for (a in 1:order){
          godfrey<-bgtest(regre, order = a, order.by, type = "Chisq")
          breush[a]<-godfrey$p.value
        }
        if (any(breush <= pvalu)) {
          repeat {
            adfresult<-ur.df(x, type = "trend", lags = b, selectlags = "Fixed")
            b<-b+1
            w <- embed(z, b)
            z.diff <- w[, 1]
            z.lag.1 <- x[b:n]
            trend <- b:n
            z.diff.lag = w[, 2:b]
            regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
            for (m in 1:order){
              godfrey<-bgtest(regre, order = m, type = "Chisq")
              breush[m]<-godfrey$p.value
            }
            if ((b == (q+1)) | (all(breush > pvalu))){
              if ((any(breush <= pvalu))) {
                warning(paste("Significant autocorrelation problems the ADF test with trend for the series", nombre, "that is growing"))
              }
              break
            }
          }
        }
        alladf<-list(series=nombre, first="Initial ADF test result for a growing series", firstBreush=purinto, firstADF=summary(purinte), final="The series is stationary around a trend", finalBreush=breush, finalADF=summary(adfresult))
        print(alladf)
      }
    } else if (status == "nogrow") {
      adfresult<-ur.df(x, type = "drift", lags, selectlags)
      z<-diff(x)
      n <- length(z)
      b<-(nrow(adfresult@testreg$coefficients)-1)
      w <- embed(z, b)
      z.diff <- w[, 1]
      z.lag.1 <- x[b:n]
      z.diff.lag = w[, 2:b]
      regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
      for (k in 1:order){
        godfrey<-bgtest(regre, order = k, order.by, type = "Chisq")
        breush[k]<-godfrey$p.value
      }
      if (any(breush <= pvalu)) {
        repeat {
          adfresult<-ur.df(x, type = "drift", lags = b, selectlags = "Fixed")
          b<-b+1
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
          for (a in 1:order){
            godfrey<-bgtest(regre, order = a, type = "Chisq")
            breush[a]<-godfrey$p.value
          }
          if ((b == (q+1)) | (all(breush > pvalu))){
            if ((any(breush <= pvalu))) {
              warning(paste("Significant autocorrelation problems in the initial ADF test for the series", nombre, "that is not growing"))
            }
            break
          }
        }
      }
      purinto<-breush
      purinte<-adfresult
      if (adfresult@teststat[2] <= adfresult@cval[4]) {
        adfresult<-ur.df(x, type = "none", lags = lagss, selectlags = selectlagss)
        z<-diff(x)
        n <- length(z)
        b<-(nrow(adfresult@testreg$coefficients)) 
        w <- embed(z, b)
        z.diff <- w[, 1]
        z.lag.1 <- x[b:n]
        z.diff.lag = w[, 2:b]
        regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
        for (j in 1:order){
          godfrey<-bgtest(regre, order = j, order.by, type = "Chisq")
          breush[j]<-godfrey$p.value
        }
        if (any(breush <= pvalu)) {
          repeat {
            adfresult<-ur.df(x, type = "none", lags = b, selectlags = "Fixed")
            b<-b+1
            w <- embed(z, b)
            z.diff <- w[, 1]
            z.lag.1 <- x[b:n]
            z.diff.lag = w[, 2:b]
            regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
            for (l in 1:order){
              godfrey<-bgtest(regre, order = l, type = "Chisq")
              breush[l]<-godfrey$p.value
            }
            if ((b == (q+1)) | (all(breush > pvalu))){
              if ((any(breush <= pvalu))) {
                warning(paste("Significant autocorrelation problems in the ADF test without drift for the series", nombre, "that is not growing"))
              }
              break
            }
          }
        }
        printuu<-breush
        summaa<-adfresult
        alladf<-list(series=nombre, first="Initial ADF test result for a not growing series", firstBreush=purinto, firstADF=summary(purinte), final="The series has a unit root without drift", finalBreush=printuu, finalADF=summary(summaa))
        print(alladf)
      } else {
        alladf<-list(series=nombre, first="Initial ADF test result for a not growing series", firstBreush=purinto, firstADF=summary(purinte), final="The series is stationary around a drift", finalBreush=purinto, finalADF=summary(purinte))
        print(alladf)
      }
    } else {
      adfresult<-ur.df(x, type = "trend", lags, selectlags)
      z<-diff(x)
      n <- length(z)
      b<-(nrow(adfresult@testreg$coefficients)-2) 
      w <- embed(z, b)
      z.diff <- w[, 1]
      z.lag.1 <- x[b:n]
      trend <- b:n
      z.diff.lag = w[, 2:b]
      regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
      for (d in 1:order){
        godfrey<-bgtest(regre, order = d, order.by, type = "Chisq")
        breush[d]<-godfrey$p.value
      }
      if (any(breush <= pvalu)) {
        repeat {
          adfresult<-ur.df(x, type = "trend", lags = b, selectlags = "Fixed")
          b<-b+1
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          trend <- b:n
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 + 1 + trend + z.diff.lag)
          for (e in 1:order){
            godfrey<-bgtest(regre, order = e, type = "Chisq")
            breush[e]<-godfrey$p.value
          }
          if ((b == (q+1)) | (all(breush > pvalu))){
            if ((any(breush <= pvalu))) {
              warning(paste("Significant autocorrelation problems in the initial ADF test for the series", nombre, "which growth status is unknoun"))
              war<-warning(paste("Significant autocorrelation problems in the ADF test for the series", nombre, "which growth status is unknown and it may be stationary around a trend"))
            }
            break
          }
        }
      }
      purinto<-breush
      purinte<-adfresult
      pri<-breush
      adfre<-adfresult
      if (adfresult@testreg[4]$coefficients[2,3] > adfresult@cval[4]) {
        z<-diff(x)
        outpute<-lm(z ~ 1)
        if (summary(outpute)$coefficients[4] <= 0.05) {
          adfresult<-ur.df(x, type = "drift", lags = lagss, selectlags = selectlagss)
          z<-diff(x)
          n <- length(z)
          b<-(nrow(adfresult@testreg$coefficients)-1) 
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          trend <- b:n
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
          for (f in 1:order){
            godfrey<-bgtest(regre, order = f, order.by, type = "Chisq")
            breush[f]<-godfrey$p.value
          }
          if (any(breush <= pvalu)) {
            repeat {
              adfresult<-ur.df(x, type = "drift", lags = b, selectlags = "Fixed")
              b<-b+1
              w <- embed(z, b)
              z.diff <- w[, 1]
              z.lag.1 <- x[b:n]
              z.diff.lag = w[, 2:b]
              regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
              for (g in 1:order){
                godfrey<-bgtest(regre, order = g, type = "Chisq")
                breush[g]<-godfrey$p.value
              }
              if ((b == (q+1)) | (all(breush > pvalu))){
                if ((any(breush <= pvalu))) {
                  warning(paste("Significant autocorrelation problems the ADF test of the series", nombre, "with unit root with drift which growing was unknoun"))
                }
                break
              }
            }
          }
          alladf<-list(series=nombre, first="The initial ADF test result for a series which growth status is unknoun", firstBreush=purinto, firstADF=summary(purinte), final="The series wich growth status was unknoun has a unit root with drift", finalBreush=breush, finalADF=summary(adfresult))
          print(alladf)
        } else {
          adfresult<-ur.df(x, type = "none", lags = lagss, selectlags = selectlagss)
          z<-diff(x)
          n <- length(z)
          b<-nrow(adfresult@testreg$coefficients) 
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          trend <- b:n
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
          for (h in 1:order){
            godfrey<-bgtest(regre, order = h, order.by, type = "Chisq")
            breush[h]<-godfrey$p.value
          }
          if (any(breush <= pvalu)) {
            repeat {
              adfresult<-ur.df(x, type = "none", lags = b, selectlags = "Fixed")
              b<-b+1
              w <- embed(z, b)
              z.diff <- w[, 1]
              z.lag.1 <- x[b:n]
              z.diff.lag = w[, 2:b]
              regre<- lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
              for (w in 1:order){
                godfrey<-bgtest(regre, order = w, type = "Chisq")
                breush[w]<-godfrey$p.value
              }
              if ((b == (q+1)) | (all(breush > pvalu))){
                if ((any(breush <= pvalu))) {
                  warning(paste("Significant autocorrelation problems the ADF test of the series", nombre, "with a unit root without drift which growth was unknoun"))
                }
                break
              }
            }
          }
          alladf<-list(series=nombre, first="The initial ADF test result for a series which growth status is unknoun", firstBreush=purinto, firstADF=summary(purinte), final="The series which growth status was unknoun has a unit root without drift", finalBreush=breush, finalADF=summary(adfresult))
          print(alladf)
        }
      } else {
        if (adfresult@testreg[4]$coefficients[3,4] <= 0.05) {
          alladf<-list(series=nombre, first="The initial ADF test result for a series which growth status is unknoun", firstBreush=purinto, firstADF=summary(purinte), final="The series which growth status was unknoun does not have a unit root but has a trend", finalBreush=pri, finalADF=summary(adfre))
          print(alladf)
        } else {
          adfresult<-ur.df(x, type = "drift", lags = lagss, selectlags = selectlagss)
          z<-diff(x)
          n <- length(z)
          b<-(nrow(adfresult@testreg$coefficients)-1)
          w <- embed(z, b)
          z.diff <- w[, 1]
          z.lag.1 <- x[b:n]
          trend <- b:n
          z.diff.lag = w[, 2:b]
          regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
          for (u in 1:order){
            godfrey<-bgtest(regre, order = u, order.by, type = "Chisq")
            breush[u]<-godfrey$p.value
          }
          if (any(breush <= pvalu)) {
            repeat {
              adfresult<-ur.df(x, type = "drift", lags = b, selectlags = "Fixed")
              b<-b+1
              w <- embed(z, b)
              z.diff <- w[, 1]
              z.lag.1 <- x[b:n]
              z.diff.lag = w[, 2:b]
              regre<- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
              for (o in 1:order){
                godfrey<-bgtest(regre, order = o, type = "Chisq")
                breush[o]<-godfrey$p.value
              }
              if ((b == (q+1)) | (all(breush > pvalu))){
                if ((any(breush <= pvalu))) {
                  warning(paste("Significant autocorrelation problems the ADF test of the series", nombre, "without a unit root and without drift which growth was unknoun"))
                }
                break
              }
            }
          }
          alladf<-list(series=nombre, first="The initial ADF test result for a series which growth status is unknoun", firstBreush=purinto, firstADF=summary(purinte), final="The series which growth status was unknoun does not have a unit root but has a drift", finalBreush=breush, finalADF=summary(adfresult))
          print(alladf)
        }
      } 
    }
  }
}