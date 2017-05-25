#'Beveridge & Nelson Decomposition
#'
#'Implements Beveridge and Nelson decomposition of series into a permanent
#'component and a cyclical one.
#'
#' @param y:   the series to decompose
#' @param pm:  maximum lag for AR (default=2)
#' @param qm:  maximum lag for MA (default=2)
#'
#' @return bnd returns a 2-column matrix of, respectively, the permanent and cyclical components in y
#' @examples #' bndco(gl)   # Quarterly US log-GDP(1947:2008)
#' y <- rnorm(1:100)   # Level-Stationary --> Stop
#' y <- cumsum(rnorm(1:100)) + (1:100)**3   # Non-Stationary 1st difference --> Stop
#'
#' @format
#' \describe{
#'   \item{}{log(RGDP)}
#' }
#' @source Example data (rlgdp):log series of U.S Real GDP (1947:1 to 2008:2)
#'  dowloaded from \url{http://time-series.net/data_sets}
#'
#' @author Deepankar Basu, \email{dbasu@@econs.umass.edu}
#' @author Fouad Ferdi, \email{fouad.ferdi@@univ-paris13.fr}


#as.person(c(
 # "Fouad Ferdi <fouad.ferdi@univ-paris13.fr> [aut, cre]",
  #"Deepankar Basu <dbasu@econs.umass.edu> [aut]"))

#' @export
bndco<- function(mod_data,pm=2,qm=2){

  #Differentiate the series
  mod_dataD<- diff(mod_data)

  # Plot the series for visual examination
  ts.plot(mod_data,xlab="time",ylab="Data",col="blue")
  abline(h=0,col="red")

  # Ask if a trend should be included in ??
  answer<-readline("Any visible time trend ? Y/N : ")
  if (substr(answer, 1, 1) == "y"){
    # Trend : Test for stationarity (ADF+KPSS)
              #tadf<-ur.df(mod_data, type="trend",selectlags="AIC")
              #tkpss<-ur.kpss(mod_data, type="tau",lags="long")

    walou<-readline("answer is yes : Trend in ADF ...")

    tadf<-tseries::adf.test(mod_data)
    tkpss<-tseries::kpss.test(mod_data,"Trend")
  }
  else {
    # No-Trend : Test for stationarity (ADF+KPSS)
                #tadf<-ur.df(mod_data, type="drift",selectlags="AIC")
                #tkpss<-ur.kpss(mod_data, type="mu",lags="long")

    walou<-readline("answer is no : notrend in ADF ...")

    tadf<-tseries::adf.test(mod_data)
    tkpss<-tseries::kpss.test(mod_data)
  }



  walou<-readline("Check stationarity test...")

  # non-stationary levels : then 1st-level differentiation
  if ((tadf$p.value)>.05 & (tkpss$p.value)<.05){
     # Test for difference-stationnarity
                  #tadf<-ur.df(mod_dataD, type="drift",selectlags="AIC")
                  #tkpss<-ur.kpss(mod_dataD, type="mu",lags="long")
    tadfD<-tseries::adf.test(mod_dataD)
    tkpssD<-tseries::kpss.test(mod_dataD,"Trend")
     # if 1st-level differences are non-stationary then STOP
    if ((tadfD$p.value)>.05 & (tkpssD$p.value)<.05) return("STOP! the series are non-stationnary in first difference")
  }
  # stationnary Level- thus use LEVELS (ANY SENS ??? ########### ?????)
    if ((tadf$p.value)<.05 & (tkpss$p.value)>.05) return("STOP! Series Level-stationary: decomposition still relevant ???" )




  # Order Selection ---------------------------------------------------------
  pvmax<-0
  for(i in 0:pm){
    for(j in 0:qm){
      fit <- arima(mod_dataD, order=c(i,0,j),method = "ML")
      # p-value of Lung-Box t-stat Ho: independance of residuals
      pvalij <- Box.test(fit$residuals,lag = 8, type="Ljung-Box",fitdf=i+j+1)$p.value
      # Computed AIC
      aicij <--2*fit$loglik+(2*(length(fit$coef)+1))
      # Computed BIC
      bicij <--2*fit$loglik+log(fit$nobs)*(length(fit$coef)+1)
      text1<-paste0("ARMA(",i,",",j,")-> ","pv:",round(pvalij,digits = 3),
                    " - aic:",round(aicij,digits = 3),
                    " - AIC:",round(fit$aic,digits = 3),
                    " - bic:",round(bicij,digits = 3),
                    " - BIC:",round(BIC(fit),digits = 3),
                    "\n")
      cat(text1)

      if (pvalij > pvmax){
        p<-i
        q<-j
        pvmax <- pvalij
      }
    }
  }
  text2<-paste0("\n","Best model fit is ARMA(",p,",",q,")","\n")
  cat(text2)

  #Arma regression
  m_arma<-arma(mod_dataD,order = c(p,q))
  summary(m_arma)
  # Plot Diff. Series
  plot.default(mod_dataD,type = "l")
  walou<-readline("Press any key to resume ...")

  # get theta & phi
  m_coef<-coef(m_arma)
  m_coef

  # compute psi
  if (p==0) PHI1<-1 else PHI1<- 1-sum(m_coef[1:p])
  MU<- m_coef["intercept"]/(PHI1)
  if (q==0) THETA1<-1 else THETA1<- 1+sum(m_coef[(p+1):(length(m_coef)-1)])              #[(paste0("ma",1)):(paste0("ma",q))])
  PSI1 = THETA1/PHI1

  # get residuals
  m_res<-residuals(m_arma)

  # FAO: First Available Observation
  fao<-max(p,q)+1

  # Set initial value
  mod_data0 <- mod_data[fao]

  # compute YP(t)
  PermComp<- MU*(fao:(length(mod_dataD)))+
    (PSI1*cumsum(m_res[(fao:(length(mod_dataD)))]))+
    mod_data0

  # compute YS(t)
  StochComp <- mod_data[-(1:fao)]-PermComp

  # plot Permanent vs. Stochastic components
  par(mfrow=c(2,1))
  plot.default(PermComp,type = "l")
  plot.default(StochComp,type = "l")
  par(mfrow=c(1,1))
  return(invisible(cbind(PermComp,StochComp)))
}

