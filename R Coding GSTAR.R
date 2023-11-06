GSTAR<-function(Data,W,p){
  Series<-ts(Data)
  k<-ncol(Data) # Number Location
  n<-nrow(Data) # Number Time
  zt<-(stack(as.data.frame(t(Data)))[,1]) # Zt with Lag
  
  # Function Create Lag Data
  ###########################################
  MD<-function(zt){
    M<-matrix(0,(n*k),k)
    z=0
    for( i in 1:(n)){
      for(j in 1:k){
        z<-z+1
        M[z,j]<-zt[z]
      }
    }
    M
  }
  ############################################
  M1<-MD(zt)
  MA<-matrix(0,(n*k-k*p),k*p)
  MAW<-matrix(0,(n*k-k*p),k*p)
  W1<-kronecker(diag(n), W)
  ztw<-W1%*%zt
  M2<-MD(ztw)
  zt<-as.matrix(zt)
  ZT<-zt[-(1:(k*p)),]
  for (i in 1:p){
    MA[(1:(n*k-k*p)),(k*(i-1)+1):(k*((i-1)+1))] <- M1[(k*(p-i)+1):(n*k-(k*i)),]
    MAW[(1:(n*k-k*p)),(k*(i-1)+1):(k*((i-1)+1))] <- M2[(k*(p-i)+1):(n*k-(k*i)),]
  }
  XT<-MA
  WXT<-MAW
  GSTAR<-data.frame(ZT,XT,WXT)
  GSTARfit<-lm(ZT~.-1,data=GSTAR)
  fit <-summary(GSTARfit)
  # Function Create Lag Data_diff
  ###########################################
  MD<-function(zt){
    M<-matrix(0,(n*k),k)
    z=0
    for( i in 1:(n)){
      for(j in 1:k){
        z<-z+1
        M[z,j]<-zt[z]
      }
    }
    M
  }
  ############################################
  M1<-MD(zt)
  MA<-matrix(0,(n*k-k*p),k*p)
  MAW<-matrix(0,(n*k-k*p),k*p)
  W1<-kronecker(diag(n), W)
  ztw<-W1%*%zt
  M2<-MD(ztw)
  zt<-as.matrix(zt)
  ZT<-zt[-(1:(k*p)),]
  for (i in 1:p){
    MA[(1:(n*k-k*p)),(k*(i-1)+1):(k*((i-1)+1))] <- M1[(k*(p-i)+1):(n*k-(k*i)),]
    MAW[(1:(n*k-k*p)),(k*(i-1)+1):(k*((i-1)+1))] <- M2[(k*(p-i)+1):(n*k-(k*i)),]
  }
  XT<-MA
  WXT<-MAW
  GSTAR<-data.frame(ZT,XT,WXT)
  GSTARfit<<-lm(ZT~.-1,data=GSTAR)
  MSE<-(fit$sigma)^2
  MAE<-sum(abs(GSTARfit$residuals))/(n*k-k*p)
  MAPE<-sum(abs(GSTARfit$residuals)/ZT)/(n*k-k*p)
  MADP<-sum(abs(GSTARfit$residuals))/sum(ZT)
  AIC<-AIC(GSTARfit)
  R2<-fit$r.squared
  R2a<-fit$adj.r.squared
  Coef<-GSTARfit$coefficients
  Error<-data.frame(MSE=MSE, MAE=MAE, MAPE=MAPE,
                    MADP=MADP, AIC=AIC, Rsquare=R2, AdjRSquare=R2a)
  Coefficient<<-data.frame(Coef)
  Zt<-as.data.frame(matrix(ZT,n-p,k,byrow=T))
  Zhat1<<-as.data.frame(matrix(GSTARfit$fitted.values,(n-p),k,byrow = T))
  Residual<<-as.data.frame(matrix(GSTARfit$residuals,(n-p),k,byrow = T))
  Result<<-data.frame(Zt,Zhat1, Residual)
  Result<-ts(Result)
  #windows()
  #plot(Series,plot.type="single", lty=1:3, col = 4:2,xlab="Time",main="Data Series Plot")
  #windows()
  #plot(Result,plot.type="multiple", lty=1:3, col = 4:2,main="Predictive vs Observed")
  #windows()
  #plot(Result,plot.type="single", lty=1:3, col = 4:2,main="Predictive vs Observed")
  Residual<-ts(Residual)
  #windows()
  #plot(Residual,plot.type="single", lty=1:3, col =4:2,main="Residual")
  Summary<-list(GOF=Error,Fit=fit)
  return(Summary)
}

#' Select the Best p
#' @param p
#' @return p
#' @export

Model<-function(p){
  Error<-matrix(0,p,7)
  for (i in 1:p){
    Error[i,]<-as.matrix(GSTAR(Data, W,i)[[1]])
  }
  colnames(Error)<-c("MSE", "MAE", "MAPE", "MADP","AIC", "Rsquare", "AdjRSquare")
  Result<-list(Error=Error)
  return(Result)
}

## fungsi rootMSE
RootMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

#### FUNGSI GSTAR UNTUK KASUS 3 LOKASI STASIUN CUACA
Data_Tot <- read.csv("dataset/origin/data.csv", header=TRUE)
Data_Asli <- as.matrix(Data_Tot)
Data_Tot <- as.matrix(Data_Tot)
dif <- 12
Data_Tot <- as.matrix(diff(Data_Tot, lag=dif, differences=1))
n_train <- floor(nrow(Data_Asli)*0.75)
Data <- as.matrix(Data_Tot[1:(n_train-dif), ])

### Stasioneritas dan Diagnosa model
library(tseries)
for (lok in 1:3)
{
  ADF_data_asli <- adf.test(Data_Asli[ ,lok], k=24)
  Test_ADF <- adf.test(Data[ ,lok], k=24)
  print(ADF_data_asli)
  print(Test_ADF)
  AutoCorrelation <- acf(Data[ ,lok], plot = FALSE)
  plot(AutoCorrelation, main = "Curah Hujan di 3 lokasi Stasiun Cuaca", lwd=3, ci.type = "ma")
  Partial_acf <- pacf(Data[ , lok], plot=FALSE)
  plot(Partial_acf, main = "CUran Hujan di 3 lokasi Stasiun Cuaca", lwd=3)
}
W<-matrix(c(0,0.5,0.5,0.5,0,0.5,0.5,0.5,0),3,3)


## Perhitungan weight berdasar normalisasi korelasi
r12 <- 0.6896245
r13 <- 0.7749152
r23 <- 0.8109261

W <- matrix(c(0,r12/(r12+r13), r13/(r12+r13), r12/(r12+r13), 0, r23/(r12+r23), r13/(r13+r12), r23/(r12+r23),0),3,3)

p <- 2
gstar <- GSTAR(Data, W, p)

Train_Result_temp <- array(data=NA, dim=c(nrow(Result), 3, 3), dimnames=NULL)
Train_Result_GSTAR <- array(data=NA, dim=c(n_train, 3, 3), dimnames=NULL)
for (lok in 1:3){                          
  for (col in 1:3){
     Train_Result_temp[ ,col,lok] <- Result[ ,(lok+(col-1)*3)]
  }
}
## result adalah hasil training gstar berdimensi jumlah_train x 9
## kolom 1-3 = data, kolom 4-6=prediksi, kol 7-9= error

### Tambah Data unt Initial
Tambah <- array(data=NA, dim=c(dif+p, 3, 3), dimnames=NULL)
for (lok in 1:3){                          
  Tambah[ , 1,lok] <- Data_Asli[1:(dif+p),lok] 
  Tambah[ , 2,lok] <- Data_Asli[1:(dif+p),lok]
}

## Invers Difference dari Train_Result_GSTAR
for (lok in 1:3) {
  Train_Result_GSTAR[(p+1):n_train,1,lok]<-diffinv(Train_Result_temp[ ,1,lok], lag=dif, diffrerences = 1, xi=Tambah[(p+1):(dif+p),1,lok])
  Train_Result_GSTAR[1:p,1,lok] <- Tambah[1:p,1,lok]
  Train_Result_GSTAR[(p+1):n_train,2,lok]<-diffinv(Train_Result_temp[ ,2,lok], lag=dif, diffrerences = 1, xi=Tambah[(p+1):(dif+p),2,lok])
  Train_Result_GSTAR[1:p,2,lok] <- Tambah[1:p,2,lok]
  Train_Result_GSTAR[ ,3,lok] <- Train_Result_GSTAR[ ,1,lok] - Train_Result_GSTAR[ ,2,lok]
}

## Train_Result_GSTAR adalah set training yang lengkap

#Menyiapkan Data Testing
n_validation <- nrow(Data_Asli)-n_train
Data_Test <- as.matrix(Data_Tot[(n_train+1):nrow(Data_Tot), ])

## transformasi koef model ke matrix
Coefficient <- as.data.frame(Coefficient)
Coef_par <- array(data=NA, dim=c(3, 2, p), dimnames=NULL)
Matrix_par1 <- array(data=NA, dim=c(3, 3, p), dimnames=NULL)
Matrix_par <- array(data=NA, dim=c(3, 3, p, 2), dimnames=NULL)
for (kol in 1:2)
{
  for (i in 1:p)
  {
    Coef_par[ ,1, i] <- Coefficient[((i-1)*3+1):(i*3),1]
    Coef_par[ ,2, i] <- Coefficient[((i-1)*3+(3*p+1)):(i*3+3*p),1]
    for (j in 1:3)
    {
      Matrix_par1[j, j, i] <- Coef_par[j,kol, i]
    }
  }
  Matrix_par[, , ,kol] <- Matrix_par1
}
Matrix_par[is.na(Matrix_par)]=0


##Proses prediksi 1 step ke depan
Data_prediksi <- matrix(data=NA, (n_validation-dif), 3, dimnames=NULL)
for (i in (p+1):(n_validation-dif))
{
  temp <- matrix(data=0, 3, 1)  
    for (j in 1:p)
    {
      temp <- temp + Matrix_par[ , , j, 1]%*% t(t(Data_Test[i-j, ])) + Matrix_par[ , , j, 2]%*% W %*% t(t(Data_Test[i-j, ]))
    }
  Data_prediksi[i, ]<-t(temp)
}

## Data prediksi 1 step, periode 1 s/d p
for (i in 1:p)
{
  Data_prediksi[i, ]<- Data_Tot[n_train+i, ]
}

##Data_prediksi[is.na(Data_prediksi)]=0

## MELAKUKAN INVERS DIFF dATuntuk perhitungan kinerja
Valid_Result_GSTAR <- array(data=NA, dim=c(n_validation, 3, 3), dimnames=NULL)
for (lok in 1:3)
{ 
  Valid_Result_GSTAR[,1,lok]<-diffinv(Data_Test[ ,lok], lag=dif, diffrerences = 1, xi=Data_Asli[(n_train+1):(n_train+dif),lok])
  Valid_Result_GSTAR[,2,lok]<-diffinv(Data_prediksi[,lok], lag=dif, diffrerences = 1, xi=Data_Asli[(n_train+1):(n_train+dif),lok])
  Valid_Result_GSTAR[,3,lok] <- Valid_Result_GSTAR[,1,lok]-Valid_Result_GSTAR[,2,lok]
  }


## Perhitungan Kinerja Training dan Validasi per lokasi KESELURUHAN DATA
korelasi <- matrix(data=NA, 3, 2)
R_sq <- matrix(data=NA, 3, 2)
RMSE <- matrix(data=NA, 3, 2)

for (lok in 1:3)
{korelasi[lok, 1] <- cor(Train_Result_GSTAR[ ,1,lok], Train_Result_GSTAR[ ,2,lok])
R_sq[lok, 1] <- (korelasi[lok,1])^2
RMSE[lok, 1] <- RootMSE(Train_Result_GSTAR[ ,1,lok], Train_Result_GSTAR[ ,2, lok])
}
for (lok in 1:3)
{
  korelasi[lok, 2] <- cor(Valid_Result_GSTAR[(p+dif+1):n_validation,1, lok], Valid_Result_GSTAR[(p+dif+1):n_validation,2,lok])
  R_sq[lok, 2] <- (korelasi[lok,2])^2
  RMSE[lok, 2] <- RootMSE(Valid_Result_GSTAR[(p+dif+1):n_validation,1, lok], Valid_Result_GSTAR[(p+dif+1):n_validation ,2, lok]) 
}

## Perhitungan Kinerja Gabungan
correlation <- matrix(data=NA, 1, 2)
RMSE_total<- matrix(data=NA, 1, 2)
R_sq_total <- matrix(data=NA, 1, 2)

for (kol in 1:2)
{
  correlation [1, kol]<- mean(korelasi[ ,kol])
  R_sq_total[1, kol]  <- mean(R_sq[ , kol])
  RMSE_total[1, kol]  <- sum(RMSE[ , kol])
}

## Perhitungan AIC correction
## fungsi AIC correction
AICC = function(MSE, n, p)
{
  n*log(MSE)+2*(p+1) + 2*(p+1)*(p+2)/(n-p)
}

AIC_correction <-matrix(data=NA, 3,1)
N <- n_train
for (lok in 1:3)
{
  MSE_value = (RootMSE(Train_Result_GSTAR[,1,lok], Train_Result_GSTAR[,2,lok]))^2
  AIC_correction[lok,1] <- AICC(MSE_value,N, p)
}
AIC_correction <-(AIC_correction[1,1]+AIC_correction[2,1]+AIC_correction[3,1])/3
print(correlation)
print(R_sq_total)
print(RMSE_total)
print(AIC_correction)


#### Plot Aktual vs Training Dari Data Differencing
for (lok in 1:3)
{
  plot(Train_Result_GSTAR[(p+1):n_train,1, lok], col="blue", type="l",main =c("Training GSTAR: Data Aktual vs Prediksi", lok), xlab="bulan ke-", ylab="curah hujan (mm)", lwd = 2) 
  lines(Train_Result_GSTAR[(p+1):n_train,2,lok], type = "l", col = "red", lwd = 1)
  legend("topleft", legend=c("aktual", "prediksi"),
         col=c("blue","red"),lty=1:0.5, cex=0.6, box.lty=1, box.lwd=0.5)
}

#### Plot Aktual vs Testing
for (lok in 1:3)
{
  plot(Valid_Result_GSTAR[(p+dif+1):n_validation,1, lok], col="blue", type="l",main =c("Testing GSTAR: Data Aktual vs Prediksi", lok), xlab="bulan ke-", ylab="curah hujan (mm)", lwd = 2) 
  lines(Valid_Result_GSTAR[(p+dif+1):n_validation,2,lok], type = "l", col = "red", lwd = 1)
  legend("topleft", legend=c("aktual", "prediksi"),
         col=c("blue","red"),lty=1:0.5, cex=0.6, box.lty=1, box.lwd=0.5)
}


## Analisis Residual
for (lok in 1:3)
{
AutoCorrelation <- acf(Train_Result_temp[ ,3,lok], plot = FALSE, lag.max = 24)
plot(AutoCorrelation, main = c("Error dari Model GSTAR", lok), lwd=3,  ci.type = "ma")
}
