setwd("/Users/rayneyang/Documents/work/test")

library(rmgarch)
library(rugarch)
library(parallel)
library(openxlsx)
library(tseries)
library(FinTS)
library(forecast)
library(vars)
library(timeDate)
library(timeSeries)
library(fBasics)
library(fGarch)
dat1 = read.csv("coin.csv",header=T)
head(dat1)
x = dat1[,2:length(dat1)]
name = c("BTC","ETH","XRP","LTC","XLM")
colnames(x) = name
head(x)
p1cor=cor(x)             ##correlation coefficient matrix
p1cor


#################Descriptive Statistics################################

png(file ="norm.jpg",width=1000,height=750)
par(mfrow=c(3,3))
for( i in 1:length(x[1,])){
	hist(x[,i],breaks =seq(from=-80,to=80),ann = F)
	par(new=T)
	curve( dnorm(x,0),ann = F, xlim=c(-10,10),xaxt="n",yaxt="n",col='red', lwd=2 )
	title(main=name[i],cex.main=2)	
}
dev.off()

#Jarque-bera test
X_squared = 0
df=0
JB_p=0
for( i in 1:length(x)){
  JB = jarque.bera.test(x[,i])
  X_squared[i] = JB$statistic[1]
  JB_p[i] = JB$p.value[1]
}
JB_result = cbind(name,X_squared,JB_p)
write.csv(JB_result,"Jarque-bera.csv")

#ADF test
Dickey_Fuller = 0
ADF_p = 0
for( i in 1:length(x)){
	adf = adf.test(x[,i])
	Dickey_Fuller[i] = adf$statistic[1]
	ADF_p[i] = adf$p.value[1]
}
adf_result = cbind(name,Dickey_Fuller,ADF_p)
write.csv(adf_result,"adf.csv")


# PP test
Dickey_Fuller = 0
PP_p = 0
for( i in 1:length(x)){
	pp = PP.test(x[,i])
	Dickey_Fuller[i] = pp$statistic[1]
	PP_p[i] = pp$p.value[1]
}
PP_result = cbind(name,Dickey_Fuller,PP_p)
write.csv(PP_result,"PP.csv")

# LM test
Lmtest = function(series,nm,max_lag){
	chi = 0
	lm_p = 0
	chi[1] = nm
	lm_p[1] = nm
	chi[2] = "chi"
	lm_p[2] = "p"
	for( lags in 1:max_lag)
	  {
		lm = ArchTest(x=series,lag=lags)
		chi[lags+2] = lm$statistic[1]
		lm_p[lags+2] = lm$p.value[1]}
	once_mat = rbind(chi,lm_p)
}

max_lag = 5
lm_mat = matrix(nrow=(2*length(x)),ncol=(max_lag+2))  
for (j in 1:length(x)){	
	start = 1+2*(j-1)
	end = 2*j
	lm_mat[start:end,] = Lmtest(x[,j],name[j],max_lag)
}
colnames(lm_mat) = append(c("cryptocurrency","parameter"),seq(from = 1, to = max_lag, by = 1))
write.csv(lm_mat,"lm_mat.csv")
#Box.test (x[1], lag = 10, type = "Ljung")
# Q-stat Ljung-Box
Boxtest = function(series,nm,max_lag){
	X_squared = 0
	box_p = 0
	X_squared[1] = nm
	box_p[1] = nm
	X_squared[2] = "X_2"
	box_p[2] = "p"
	for( lags in 1:max_lag){
		x.fit = garch(x = series,order=c(0,lags))
		res_sqr = x.fit$residuals**2
		box = Box.test (res_sqr, lag = 1, type = "Ljung")
		X_squared[lags+2] = box$statistic[1]
		box_p[lags+2] = box$p.value[1]
	}
	once_mat = rbind(X_squared,box_p)
}

max_lag = 5
box_mat = matrix(nrow=(2*length(x)),ncol=(max_lag+2))  
for (j in 1:length(x)){	
	start = 1+2*(j-1)
	end = 2*j
	box_mat[start:end,] = Boxtest(x[,j],name[j],max_lag)
}
colnames(box_mat) = append(c("coin","parameter"),seq(from = 1, to = max_lag, by = 1))
write.csv(box_mat,"box_mat.csv")
#################################################################

############################var##########################

dat1 <- read.csv("coin.csv",header=T)
date=as.Date(as.character(dat1[,1]),"%d/%m/%Y")
x <- abs(dat1[,-1])
#View(x)

colnames(x) = c("BTC","ETH","XRP","LTC","XLM")
N=nrow(x)
dim(x)

VARselect(x, lag.max = 5, type="const")
var1<-VAR(x, p = 1, type = "const")
summary(var1)

varfit<-varxfit(x,1, constant = TRUE,postpad = c("constant"))
show(varfit)
################################################################################


######################var-dcc-garch#############################################

spec <- ugarchspec( variance.model = list(model = "sGARCH", garchOrder = c(1, 1)) ,
                    mean.model = list(armaOrder = c(1, 1), include.mean = FALSE,archm = TRUE))


muspec = multispec(replicate(spec,n=5))
fitlist = multifit(multispec = muspec, data = x , out.sample = 0, solver = list("solnp"))
spec1 = dccspec(uspec = muspec,VAR=TRUE,lag=1, 
                dccOrder = c(1, 1), 
                model=c("DCC"),
                distribution = "mvnorm")

fit1 = dccfit(spec1, x, fit.control = list(eval.se = TRUE), 
              VAR.fit = varfit,solver = "solnp", out.sample = 0)
show(fit1)

#################################################################################


##############dcc2#################
dcc_once=function(xx){
  spec1 = ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)), mean.model = list(armaOrder=c(1,1)))
  spec2 = ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)), mean.model = list(armaOrder=c(1,1)))
  mspec = multispec( c( spec1, spec2 ) )
  spec = dccspec(mspec, VAR = FALSE, robust = FALSE,external.regressors = NULL, 
                 dccOrder = c(1,1), model = "DCC", distribution = "mvnorm") 
  m1=dccfit(spec, xx, out.sample = 0, solver = "solnp", solver.control = list(), fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL)
  show(m1)
  matcoef=m1@mfit$matcoef 
  rho=m1@mfit$R
  n=length(rho)
  rho.t=0
  for(i in 1:n){
    r=rho[[i]]
    rho.t[i]=r[1,2]
  }
  return(rho.t)
}
cormat = matrix(nrow=length(x[,1]),ncol=((length(x)-1)))  
for( i in 2:length(x)){
  xx = cbind(x[,1],x[,i])
  cormat[,(i-1)] = dcc_once(xx)
}

png(file ="BTC.jpg",width=750,height=750)
par(mfrow=c(2,2))
for( i in 1:length(cormat[1,])){
  plot(ts(cormat[,i],start=c(2015,9),frequency = 250), ylim = c(-1, 1),type = "l",lty=1,ylab= 'correlation',xlab= 'Time')
  title(main=paste(main=paste("BTC & ",name[i+1],sep="")),cex.main=2)	
}
dev.off()

write.csv(cormat,"cormat_ne.csv")

data <- read.csv("cormat_ne.csv",header = T)
head(data)
data$Date = as.Date(data$Date)
ggplot(data,aes(Date,ne_it))+geom_line()+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

ggplot(data,aes(Date,oil_it))+geom_line()+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())


BTC_arma <- auto.arima(x[1], stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##btc ARMA(1,2)
ETH_arma <- auto.arima(x[2], stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##eth ARMA(1,1)
XRP_arma <- auto.arima(x[3], stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##xrp  ARMA(2,2)
LTC_arma <- auto.arima(x[4], stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##ltc ARMA(4,3)
XLM_arma <- auto.arima(x[5], stationary = T, seasonal = F, ic = "aic", allowdrift = FALSE, trace = TRUE)
##xlm ARMA(5,0)
checkresiduals(BTC_arma)
shapiro.test(BTC_arma$residuals)
##The residuals are not normally distributed.



r_BTC <- BTC_arma$residuals
r_ETH <- ETH_arma$residuals
r_XRP <- XRP_arma$residuals
r_LTC <- LTC_arma$residuals
r_XLM <- XLM_arma$residuals

ugarch.spec1 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(1,2)),
                           distribution.model = "sstd")

BTC_garch <- ugarchfit(spec = ugarch.spec1, data = x[1])
BTC_garch
ugarch.spec2 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(1,1)),
                           distribution.model = "sstd")

ETH_garch <- ugarchfit(spec = ugarch.spec2, data = x[2])
ETH_garch

ugarch.spec3 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(2,2)),
                           distribution.model = "sstd")

XRP_garch <- ugarchfit(spec = ugarch.spec3, data = x[3])
XRP_garch

ugarch.spec4 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(4,3)),
                           distribution.model = "sstd")

LTC_garch <- ugarchfit(spec = ugarch.spec4, data = x[4])
LTC_garch

ugarch.spec5 <- ugarchspec(variance.model=list(garchOrder=c(1,1)), 
                           mean.model=list(armaOrder=c(5,0)),
                           distribution.model = "sstd")

XLM_garch <- ugarchfit(spec = ugarch.spec2, data = x[5])
XLM_garch

dcc.garch.spec = dccspec(uspec = multispec(c(ugarch.spec1, 
                                             ugarch.spec2, 
                                             ugarch.spec3, 
                                             ugarch.spec4,
                                             ugarch.spec5)),
                         dccOrder = c(1,1),
                         distribution = "mvt")

dcc.garch.spec
dcc.fit = dccfit(dcc.garch.spec, data = x)
show(dcc.fit)
