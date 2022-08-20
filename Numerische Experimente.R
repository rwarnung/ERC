#f?r Inverse
library(matlib)
#f?r Rechenaufwand
library(microbenchmark)
#f?r Kovarianzmatrix
library(corpcor)
library(tseries)
library(quantmod)
#Excel einlesen
library(xlsx)

library(tidyquant)
#Zusammensetzung eines Index direkt von yahoo, zum Beispiel df_sp500 <- yf_index_composition("SP500")
library(yfR)

x = read.xlsx2("C:\\Users\\bodoj\\OneDrive\\Dokumente\\Mathematik\\Bachelorarbeit\\Stocks in the Dow Jones Industrial Average.xlsx","Sheet1",colClasses=c("character"),startRow = 2,colIndex = 1,stringsAsFactors=FALSE)
#x = read.xlsx("C:\\Users\\bodoj\\OneDrive\\Dokumente\\Mathematik\\Bachelorarbeit\\Stocks in the NASDAQ-100 Index.xlsx","Sheet1",colClasses=c("character"),startRow = 2,colIndex = 1,stringsAsFactors=FALSE)
#x = read.xlsx("C:\\Users\\bodoj\\OneDrive\\Dokumente\\Mathematik\\Bachelorarbeit\\Stocks in the SP 500 Index.xlsx","Sheet1",colClasses=c("character"),startRow = 2,colIndex = 1,stringsAsFactors=FALSE)
#x = read.xlsx("C:\\Users\\bodoj\\OneDrive\\Dokumente\\Mathematik\\Bachelorarbeit\\Stocks in the Russell 1000 Index.xlsx","Sheet1",colClasses=c("character"),startRow = 2,colIndex = 1,stringsAsFactors=FALSE)
#x = read.xlsx("C:\\Users\\bodoj\\OneDrive\\Dokumente\\Mathematik\\Bachelorarbeit\\Stocks in the Russell 2000 Index.xlsx","Sheet1",colClasses=c("character"),startRow = 2,colIndex = 1,stringsAsFactors=FALSE)

covmatrix = function(x){
  c = tq_get(x[1,1])#Daten von erster Aktie
  c = tail(c,780)#52 Wochen im Jahr, 3 Jahre, das sind 156, wollen wöchentliche Daten (5 Handelstage die Woche), damit erhalten wir 156*5 = 780, nehmen daher die letzten 780 Daten
  c = data.frame(c[seq(1,nrow(c),5),"close"])#alle 5 Handelstage (also wöchentlich, zirka)
  c = diff(log(c[,1]),1)#berechne Differenz
  c = data.frame(c)#erstelle Liste
  portfolio = c
  j = 0#werden dann Daten aus x löschen, j zählt mit, wie viele gelöscht wurden
  for(i in 2:length(x[,1])) {#mit for-Schleife gehts weiter f?r restliche Aktien
    c = tq_get(x[i-j,1])
    if(is.null(ncol(c))){#falls Daten nicht vorhanden
       x = x[-(i-j),1]
       x = data.frame(x)
       j = j+1
      next
    }
    c = tail(c,780)
    c = c[seq(1,nrow(c),5),"close"]#wieder nehmen wir jedes 5. Element
    if(nrow(c) < 156){#nicht bei jeder Aktie sind die letzten 3 Jahre vorhanden, in diesem Fall wird sie ausgelassen
      x = x[-(i-j),1]
      x = data.frame(x)
      j = j+1
      next
    }
    c = c$close#wir fügen die entsprechende Zeile der data frame hinzu
    c = diff(log(c),1)
    portfolio[,ncol(portfolio) + 1] <- c#f?ge immer eine weitere Spalte hinzu
  }
  names(portfolio) = unlist(x)#geben Namen der Spalten, ist zwar unn?tig, ich lasse es aber jetzt noch so, bei Tests besser lesbar
  covariance = matrix(c(cov(portfolio)), nrow = nrow(x), ncol = nrow(x))#erstelle Kovarianzmatrix, Dimension!
  dimnames(covariance) = list(unlist(x), unlist(x))#benenne Zeilen und Spalten der Matrix
  # if(ncol(covariance)>=100){#Shrinkage Schätzer
  #   n = nrow(covariance)
  #   covariance = cov.shrink(covariance)
  #   covariance = covariance[1:n,1:n]
  # }
  return(covariance)
}

#erster Versuch mit CCD:
covariance = covmatrix(x)
print("covariance matrix done")
#Dimension
n = nrow(covariance)
#daraus Korrelationsmatrix
correlation = cov2cor(covariance)
#berechne mir den Anfangswert, der ist x_0 = (1/\sigma_i) / (\sum 1/\sigma_k)
sigma = diag(covariance)
sigma = sqrt(sigma)

x = (1 / sigma) / (sum(1 / sigma))
Sigmax = function(x,y,z,i){#berechnet neues Sigma*x, verbessert Rechenleistung
  return(z - covariance[,i] * y[i] + covariance[,i] * x[i])
}

sigmax = sqrt(t(x) %*% covariance %*% x)#Risikoma?

#Iteration von CCD
ccd = function(x){
  z = covariance %*% x
  RC = c(1:n) #Risikoanteil vom i-ten Asset als Vektor
  for(h in 1:n){ 
    RC[h] = x[h] * z[h]}
  while(max(abs(RC / as.vector(sigmax) - as.vector((1/n) * sigmax)))>0.0001){#Abbruchskriterium ist Definition von Risikoparit?t
    for(i in 1:n){
      y = x
      x[i] = (-z[i] + x[i] * covariance[i,i] + sqrt((z[i] - x[i] * covariance[i,i])^2 + 4 * covariance[i,i] * (1 / n) * sigmax) ) / (2 * covariance[i,i])
      z = Sigmax(x,y,z,i)
      sigmax = sqrt(sigmax^2 - 2 * y[i] * (covariance[i,] %*% y) + y[i]^2 * covariance[i,i] + 2 * x[i] * (covariance[i,] %*% x) - x[i]^2 * covariance[i,i])
    }
    for(h in 1:n){ #berechne Risikoanteil neu
      RC[h] = x[h] * z[h]}
  }
  x = x / (sum(x))
  return(x)
}

ccd = ccd(x)
print("CCD done")
#verbessertes CCD

#Anfangswert
sumcorr = function(correlation){
  c = 0
  for(i in 1:n){
    for(j in 1:n){
      c = c +correlation[i,j]
    }
  }
  return(c)
}

x = c(rep(1,n))
x = x * 1 / (sqrt(sumcorr(correlation)))

a = function(x,i){
  return(((correlation %*% x)[i] - x[i]) / 2)
}

RC = function(x){
  z = c(rep(1,n))
  for(i in 1:n){
    z[i] = x[i] * (correlation %*% x)[i]
  }
  return(z)
}

impccd = function(x){
  while(max(abs(RC(x)  - as.vector((1 / n)))) > 0.0001){
    for(i in 1:n){
      a = a(x,i)
      x[i] = sqrt(a^2 + 1 / n) - a
    }
    x = x / (as.vector(sqrt(t(x) %*% correlation %*% x)))
  }
  return((x / sigma) / sum(x / sigma))
}

impccd = impccd(x)
print("improved CCD done")

# #weiter mit Newton
#Anfangswert
x = rep(1 / (sum(correlation)),n)


lambdastar = 0.95 * (3-sqrt(5))/2

b = rep(1 / n, n)

#erste Ableitung
derivative = function(x){
  return(correlation %*% x - b * x^(-1))
}

#zweite Ableitung
twoderivative = function(x){
  return(correlation + (b %*% x^(-2))[1,1] * diag(n))
}

#Newtonverfahren
newton = function(x){
  first = derivative(x)
  second = twoderivative((x))
  bdelta = solve(second) %*% first
  ldelta = max(abs(bdelta / x))
  lambda = sqrt(t(first) %*% bdelta)
  while(lambda > lambdastar) {
    x = x - (1 / (1 + ldelta)) * bdelta
    first = derivative(x)
    second = twoderivative((x))
    bdelta = solve(second) %*% first
    ldelta = max(abs(bdelta / x))
    lambda = sqrt(t(first) %*% bdelta)
  }
  while(lambda > 0.01){ #da das lambda nur langsam fällt ist für das Verfahren entscheidend, wie die untere Schranke gewählt wird
    x = x - bdelta
    first = derivative(x)
    second = twoderivative((x))
    bdelta = solve(second) %*% first
    lambda = sqrt(t(first) %*% bdelta)
  }
  y = x
  for(i in 1:n){
    c = 0
    for(j in 1:n){
      c = c +sqrt(covariance[j,j])^(-1)*y[j]
    }
    x[i] = sqrt(covariance[i,i])^(-1) * y[i] / c
  }
  return(x)
}

newton = newton(x)
print("newton done")

#verbessertes Newton
#Anfangswert
x = rep(1 / (sum(correlation)),n)

impnewton = function(x){
  a = (correlation %*% c(rep(1,n)) - c(rep(1,n))) / as.vector(2 * sqrt(t(c(rep(1,n))) %*% correlation %*% c(rep(1,n))))
  x = sqrt(a * a + b) - a#eine ccd Iteration
  first = derivative(x)
  second = twoderivative((x))
  bdelta = solve(second) %*% first
  lambda = sqrt(t(first) %*% bdelta)
  #normales newton
  while(lambda > 0.01){
    x = x - bdelta
    first = derivative(x)
    second = twoderivative((x))
    bdelta = solve(second) %*% first
    lambda = sqrt(t(first) %*% bdelta)
  }
  y = x
  for(i in 1:n){
    c = 0
    for(j in 1:n){
      c = c +sqrt(covariance[j,j])^(-1)*y[j]
    }
    x[i] = sqrt(covariance[i,i])^(-1) * y[i] / c
  }
  return(x)
}

impnewton = impnewton(x)
print("improved newton done")