#### PORTFOLIO ####
rm(list = ls())
library(xts)
library(grid)
library(stringr)
library(highcharter)
library(lubridate)
library(tidyverse)
library(jrvFinance)
library(RColorBrewer)
cols = brewer.pal(n = 3,name = "Set1")
# We load a dataframe of 26 US Gov Bonds retrived from MarketInsider.com
# using a Python script.
# Then, we calculate the cash flows of each bond and yield to maturity
# and dtore them in a list separately.
# We spot some differences in the yield posted and the one calculated.
# The objective is to create one unique vector of cash flows which
# should be mapped to the KEY RATES (or vertices) over which PCA has been performed.

bonds <- read.csv("/Users/Andrea/Desktop/thesis/data/bonds_dataframe.csv")
colnames(bonds)
#### Data transformation of df containing bonds ####
bonds <- bonds %>% 
  #as_tibble(bonds) %>% 
  mutate(PURCHASE_DAY = "2020-03-10") %>% 
  mutate_at(c("ISSUE_DATE","COUPON_PYMT_DATE","MATURITY","STARTCOUPON_DATE","FINALCOUPON_DATE"), as.Date, format = "%m/%d/%Y") %>% 
  mutate_at(c("ISSUE_DATE","COUPON_PYMT_DATE","MATURITY","STARTCOUPON_DATE","FINALCOUPON_DATE","PURCHASE_DAY"), as.Date, format = "%Y-%m-%d") %>%
  #mutate(FACE_VAL = FACE_VAL*10) %>% 
  #mutate(ISSUE_PRICE = ISSUE_PRICE*10) %>% 
  #mutate(MRKT_PRICE = MRKT_PRICE*10) %>% 
  mutate(DAYS_TO_COUPON = COUPON_PYMT_DATE - PURCHASE_DAY) %>% 
  #mutate(YTM = YTM/100) %>% 
  mutate(YTM_2 = sqrt(1+YTM)-1) %>% 
  mutate(SEMI_COUPON_AMOUNT = ((FACE_VAL*(CP_RATE/100))/2)) %>% # semi annual coupon payment amount
  mutate(NUMB_PAYMENTS = 2) %>% 
  rename(POSTED_YTM = YTM)

# Calculate the number of coupons remaining
for (i in 1:dim(bonds)[1]){
  bonds[i,"N_CFs"] <- jrvFinance::coupons.n(settle = bonds[i,"PURCHASE_DAY"],
                                            mature = bonds[i,"MATURITY"],
                                            freq = 2)
}

#### Create empty list for semi annual cash flows of each bond ####
cash_flows <- vector(mode = "list", length = dim(bonds)[1])
# Populate list with cash flows for each bonds. Firts is negative, last = Face val + Coupon
for (cp in 1:dim(bonds)[1]) {
  cash_flows[[cp]] <- c(-bonds[cp,c("MRKT_PRICE")], 
                        rep(bonds[cp,c("SEMI_COUPON_AMOUNT")],bonds[cp,c("N_CFs")]-1),
                        bonds[cp,c("SEMI_COUPON_AMOUNT")] + bonds[cp,c("FACE_VAL")])
}
cash_flows

#### Calculate semi annual ytm for each set of cash flows ####
yields_to_mat_semi <- c()
for (i in 1:length(cash_flows)){
  yields_to_mat_semi[i] <- jrvFinance::irr(cash_flows[[i]],
                                           interval = c(-1,1),
                                           cf.freq = 2,
                                           comp.freq = 2)
}
#((1007.5/993)^2)-1 0.02941766
# irr(c(-903.10,rep(100,3),1000),interval=c(0,1),cf.freq = 1,comp.freq = 1)
yields_to_mat_semi
bonds[,"SEMI_ANN_YTM"] <- yields_to_mat_semi
bonds <- bonds %>% 
  mutate(ANNUAL_YTM=((1+SEMI_ANN_YTM)^2)-1)

# bonds[,"DELTA_YTM"] <- round((bonds[,"ANNUAL_YTM"]*100)-(bonds[,"POSTED_YTM"]),5)
# #### Bar plot of differences in the calculated yield and posted one ####
# ggplot(bonds,aes(x = NAME, y = DELTA_YTM)) + 
#   geom_bar(stat="identity") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle=90),
#         plot.title = element_text(hjust = 0.5),
#         panel.grid.major.x = element_blank(),
#         axis.ticks.y = ) +
#   ylab("Obtained YTM - Posted YTM") +
#   xlab("Bond") +
#   ggtitle("Difference between 'obtained' yearly YTM and 'posted' YTM on Market Insider. March 10 2020")
#ggsave("/Users/Andrea/Desktop/thesis/images/delta_btw_YTM.png", device = "png", width = 297, height = 210, units = "mm")

#### Recover dates associated to KEYRATES of US TREASURY ####
# Starting from purchase day 2020-03-10
s <- ymd("2020-03-10") #start date
KEYRATES <- c(s+ddays(31),s+3*ddays(31),s+6*ddays(31), s+dyears(1), s+2*dyears(1),s+3*dyears(1),
              s+5*dyears(1),s+7*dyears(1),s+10*dyears(1),s+20*dyears(1),s+30*dyears(1))

# Create Factor from COUPON_PYMT_DATE 
# we see there are 7 common dates of first payments.
# we want 7 different groups of bonds depending on
# the their coupon payment dates
coupon_payment_date <- sort(ymd(bonds[,"COUPON_PYMT_DATE"]))
fct <- factor(coupon_payment_date, labels = c("I","II","III","IV","V","VI","VII"))
bonds <- bonds %>% 
  mutate(GROUP = fct)

# Order the seven first date of payments 
c_pmyt_dates_7 <- sort(unique(ymd(bonds[,"COUPON_PYMT_DATE"])))
names(c_pmyt_dates_7) <- levels(fct)

#### Create seven series of dates to be mapped to each group of bonds ####
# We associate a date to each cash flow of each bond.
s1 <- c_pmyt_dates_7[1]
dates_GI <- c(ymd("2020-03-10"),s1)
for(i in 1:60){
  dates_GI[i+2] <- s1+(i*(6*ddays(31)))
}

s2 <- c_pmyt_dates_7[2]
dates_GII <- c(ymd("2020-03-10"),s2) 
for(i in 1:60){
  dates_GII[i+2] <- s2+(i*(6*ddays(31)))
}

s3 <- c_pmyt_dates_7[3]
dates_GIII <- c(ymd("2020-03-10"),s3) 
for(i in 1:60){
  dates_GIII[i+2] <- s3+(i*(6*ddays(31)))
}

s4 <- c_pmyt_dates_7[4]
dates_GIV <- c(ymd("2020-03-10"),s4) 
for(i in 1:60){
  dates_GIV[i+2] <- s4+(i*(6*ddays(31)))
}

s5 <- c_pmyt_dates_7[5]
dates_GV <- c(ymd("2020-03-10"),s5) 
for(i in 1:60){
  dates_GV[i+2] <- s5+(i*(6*ddays(31)))
}

s6 <- c_pmyt_dates_7[6]
dates_GVI <- c(ymd("2020-03-10"),s6)
for(i in 1:60){
  dates_GVI[i+2] <- s6+(i*(6*ddays(31)))
}

s7 <- c_pmyt_dates_7[7]
dates_GVII <- c(ymd("2020-03-10"),s7) 
for(i in 1:60){
  dates_GVII[i+2] <- s7+(i*(6*ddays(31)))
}
#### Set the names of the elements of the list cash_flows ####  
# in order to match them with the series of dates
cash_flows_df <- lapply(cash_flows,as.data.frame)
names(cash_flows_df) <- bonds[,"GROUP"]

#### Create a unique vector which aggregates all the 26 cash flows ####
for (i in 1:length(cash_flows)) {
  if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[1]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GI[1:dim(cash_flows_df[[i]])[1]])
  } 
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[2]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GII[1:dim(cash_flows_df[[i]])[1]])
  }
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[3]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GIII[1:dim(cash_flows_df[[i]])[1]])
  }
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[4]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GIV[1:dim(cash_flows_df[[i]])[1]])
  }
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[5]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GV[1:dim(cash_flows_df[[i]])[1]])
  }
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[6]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GVI[1:dim(cash_flows_df[[i]])[1]])
  }
  else if (names(cash_flows_df)[i] == names(c_pmyt_dates_7)[7]) {
    cash_flows_df[[i]] <- xts(cash_flows_df[[i]],dates_GVII[1:dim(cash_flows_df[[i]])[1]])
  }
}
cash_flows_df
for(i in 1:length(cash_flows_df)) {
  colnames(cash_flows_df[[i]]) <- c("CF")
}

CFs_all <- merge(cash_flows_df[[1]],cash_flows_df[[2]],cash_flows_df[[3]],cash_flows_df[[4]],
                 cash_flows_df[[5]],cash_flows_df[[6]],cash_flows_df[[7]],
                 cash_flows_df[[8]],cash_flows_df[[9]],cash_flows_df[[10]],cash_flows_df[[11]],
                 cash_flows_df[[12]],
                 cash_flows_df[[13]],cash_flows_df[[14]],cash_flows_df[[15]],cash_flows_df[[16]],
                 cash_flows_df[[17]],cash_flows_df[[18]],cash_flows_df[[19]],cash_flows_df[[20]],
                 cash_flows_df[[21]],cash_flows_df[[22]],cash_flows_df[[23]],cash_flows_df[[24]],
                 cash_flows_df[[25]],cash_flows_df[[26]], all = TRUE)
dim(CFs_all)
CFS <- apply(CFs_all,1,sum, na.rm = TRUE)
#### Cash Flows distribution along the years - Graph ####
# Dates of the cash flows
dates_CFS <- as.Date(names(CFS))
CFS_df <- as.data.frame(CFS)
CFS_ts <- xts(CFS, dates_CFS)

# We add the dates corresponding to the key rates
dates_cfs_key <- sort(c(ymd(rownames(CFS_df)),KEYRATES))
length(dates_cfs_key)

CF_plot <- xts(rep(0,length(dates_cfs_key)),dates_cfs_key)
CFS_df_plot <- as.data.frame(merge(CFS_ts,CF_plot))
CFS_df_plot <- apply(CFS_df_plot,1,sum)
CFS_df_plot <- as.data.frame(CFS_df_plot)
CFS_df_plot[is.na(CFS_df_plot)] <- 0
CFS_df_plot[,"Date"] <- rownames(CFS_df_plot)
colnames(CFS_df_plot)[1] <- "CF"
head(CFS_df_plot)
str(CFS_df_plot)
rownames(CFS_df_plot) <- seq(1,dim(CFS_df_plot)[1],1)
key_rates_index <- which(CFS_df_plot[,"CF"] == 0.000)

KEY_df <- as.data.frame(KEYRATES)
ggplot(CFS_df_plot, aes(x=Date,y=CF)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=90,size=rel(0.5)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1"),
        strip.background = element_rect("#EDF0F1")) +
  ggtitle("Distribution of the portfolio cash flows along the years.\n In blue, the vertices of the cash flows map rates:\n1m, 3m, 6m, 1y, 2y, 3y, 5y, 7y, 10y, 20y, 30y") +
  ylab("($)") +
  xlab("Day") +
  geom_vline(xintercept = c(key_rates_index),linetype = 2, colour = 'blue', size = 0.25)
ggsave("/Users/Andrea/CF_series_COLOUR_CORRECT.png",units = "mm",device = "png", width = 297,height = 210, dpi = 320)
#ggsave("/Users/Andrea/Desktop/thesis/images/CF_series.png",units = "mm",device = "png", width = 297,height = 210, dpi = 320)

#### Calculate the interest rates at intermediate maturities ####
# given bonds market prices and key rates or vertices
yield_curve_20200310 <- c(0.57,0.44,0.43,0.43,0.5,0.58,0.63,0.73,0.76,1.16,1.28)
maturities.Tres <- c(1/12,3/12,6/12,1,2,3,5,7,10,20,30)
yield_curve_20200310_df <- as.data.frame(cbind(maturities.Tres,yield_curve_20200310),
                                         stringsAsFactors = FALSE)
colnames(yield_curve_20200310_df) <- c("Maturity","Yield")
yield_curve_20200310_df
ggplot(yield_curve_20200310_df,aes(x=Maturity,y=Yield)) +
  geom_line() +
  geom_point() +
  theme_light() +
  labs(y = "Yield (%)", x = "Maturity (Years)") +
  theme(plot.title = element_text(hjust= 0.5),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1")) +
  ggtitle("Spot Rate US Treasury Curve on March 10, 2020")
ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct6.png",device="png",units = "mm",height = 210,width = 297,dpi=320)

#### Svensson Parametric model for spot curve ####
# loss of accuracy a short and long terms
# less flexible than splines
# In lenghts_from_t0 we store the number of days that separates 
# each cash flow to 2020-03-10

# lenghts_from_t0 <- list()
# for(i in 1:length(dates_cfs_key)){
#   lenghts_from_t0[[i]] <- dates_cfs_key[i] - ymd("2020-03-10")
# }

lenghts_from_t0 <- list()
for(i in 1:length(dates_CFS)){
  lenghts_from_t0[[i]] <- dates_CFS[i] - ymd("2020-03-10")
}
#### INTEREST RATE DATA ####
maturities.CFS <- unlist(lenghts_from_t0)[-1]/360 # years
library(YieldCurve)
X <- read.table("/Users/Andrea/Desktop/thesis/data/USTREASURY_YIELD_UPDATED.csv")
X <- xts(X,ymd(rownames(X)))
X <- window(X, start="2006-02-09",end="2020-03-10")
dim(X)
#### Estimate of the Svensson paramters for March 10, 2020 ####
maturities.Tres <- c(1/12,3/12,6/12,1,2,3,5,7,10,20,30)
SvensonParameters <- Svensson(rate = X[dim(X)[1],], maturity = maturities.Tres)
Svensson.rate <- Srates(SvensonParameters, maturities.CFS, "Spot")
Svensson_rates_df <- as.data.frame(cbind(maturities.CFS,as.vector(Svensson.rate)))
colnames(Svensson_rates_df) <- c("Maturity","Yield")
ggplot(data=Svensson_rates_df,aes(x=Maturity,y=Yield)) +
  #geom_point(size= 0.5) +
  geom_line() +
  theme_light() +
  labs(x = "Maturity (Years)", y = "Yield (%)") +
  theme(plot.title = element_text(hjust= 0.5),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1")) +
  ggtitle("Estimated 'Svensson' Spot Curve on March 10, 2020")
ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct5.png",device="png",units = "mm",height = 210,width = 297,dpi=320)


#### Svensson: Calculate the absolute interest rates change in basis points of the interest rates estimated with Svensson's model ####
#SvensonParametersAll <- Svensson(rate = X, maturity = maturities.Tres)
#SvensonParametersAll <- as.data.frame(SvensonParametersAll)
#write.table(SvensonParametersAll,"/Users/Andrea/Desktop/thesis/data/svensson_pars_all.csv",sep=",")
SvensonParametersAll <- read.csv("/Users/Andrea/Desktop/thesis/data/svensson_pars_all.csv")
SvensonParametersAll <- xts(SvensonParametersAll,index(X))

#Svensson.ratesAll <- Srates(SvensonParametersAll, maturities.CFS, "Spot")
labelsSvenssonAll <- str_c(round(maturities.CFS,2),"YR")
#colnames(Svensson.ratesAll) <- labelsSvenssonAll
#Svensson.rateAll <- as.data.frame(Svensson.ratesAll)
#write.table(Svensson.ratesAll,"/Users/Andrea/Desktop/thesis/data/svensson_rates_all.csv",sep=",")
Svensson.ratesAll <- read.csv("/Users/Andrea/Desktop/thesis/data/svensson_rates_all.csv")
colnames(Svensson.ratesAll) <- labelsSvenssonAll
Svensson.ratesAll <- xts(Svensson.ratesAll,index(X))

#### Svensson: Absolute Interest Rates changes in Bps ####
XSvensson <- apply(Svensson.ratesAll,2,diff)*100
library(RColorBrewer)
cols = brewer.pal(n = dim(XSvensson)[2],name = "Greys")
#png("/Users/Andrea/Desktop/thesis/images/svensson_rates.png",height=210,width=297,units="mm",res=320)
plot(Svensson.ratesAll, 
     main = "Daily Yield Curve estimated with Svensson model\nacross each cash flow maturitiy of the portfolio",
     col = cols)
#dev.off()
# compare with:
#png("/Users/Andrea/Desktop/thesis/images/US_rates_XTS.png",height=210,width=297,units="mm",res=320)
plot(X, main= "Spot Interest Rates, US Treasury",
     col = cols)
#dev.off()

#### Svensson: Volatilities of Rates at different maturities - Absolute rates changes in Bps ####
XSvenssonSD <- as.data.frame(apply(XSvensson,2,sd))
colnames(XSvenssonSD) <- "Bps"
XSvenssonSD[,"Maturity"] <- labelsSvenssonAll

ggplot(XSvenssonSD, aes(x=Maturity,y=Bps,group=1)) + 
  geom_line() +
  scale_x_discrete(limits=labelsSvenssonAll) +
  theme_light() +
  ggtitle("Standard deviations in Basis Points of the Svensson Interest rates Absolute changes\nUS spot curve 2006-2020") +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=90,size=rel(0.6))) +
  ylab("Volatility (Bps)") +
  xlab("Maturity (Years)")
#ggsave("/Users/Andrea/Desktop/thesis/images/SVENSSON_volatilities_abs_changes_x_100.png",height = 210,width = 297,units = "mm",dpi = 320)

#### Svensson: Volatilities of Rates at different maturities - Svensson's Interest rates ####
Svensson.ratesAll_df_sd <- as.data.frame(apply(Svensson.ratesAll,2,sd))
Svensson.ratesAll_df_sd[,"Maturity"] <- labelsSvenssonAll
colnames(Svensson.ratesAll_df_sd)[1] <- "Sd"

ggplot(Svensson.ratesAll_df_sd, aes(x=Maturity,y=Sd,group=1)) + 
  geom_line() +
  scale_x_discrete(limits=labelsSvenssonAll) +
  theme_light() +
  ggtitle("Standard deviations in Basis Points of the Svensson Interest rates\nUS spot curve 2006-2020") +
  theme(plot.title = element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=90,size=rel(0.6))) +
  ylab("Volatility (Sd)") +
  xlab("Maturity (Years)")
#ggsave("/Users/Andrea/Desktop/thesis/images/SVENSSON_volatilities_plain_interest_rates.png",height = 210,width = 297,units = "mm",dpi = 320)

#### PCA on SVENSSON Yield Curves ####
#### PCA on the COVARIANCE matrix - Absolute rates changes in Bps ####
CovXSvensson <- cov(XSvensson)
W <- eigen(CovXSvensson,symmetric=TRUE)$vectors
Y <- XSvensson%*%W
Y <- xts(Y,index(X)[-1]) # Because of diff()
colnames(Y) <- str_c("MAT",round(maturities.CFS,2))

#### Graph First PC Cov. Mat. - Absolute rates changes in Bps ####
highchart(type = "stock") %>% 
  hc_title(text = "First three PCs - Cov. Matrix - Svensson Absolute rates changes in Bps") %>% 
  hc_add_series(Y[,1], name = "PC1") %>% 
  hc_add_series(Y[,2], name = "PC2") %>%   
  hc_add_series(Y[,3], name = "PC3") %>%
  hc_add_theme(hc_theme_flat()) %>% 
  hc_navigator(enabled = FALSE) %>% 
  hc_scrollbar(enabled = FALSE) %>% 
  hc_exporting(enabled = TRUE) %>% 
  hc_legend(enabled = TRUE) 

#### Graph of the Eigenvectors estiamted on COVARIANCE ####
# Absolute changes in basis points #
maturities_labs <- list(rep(round(maturities.CFS,2),3)) # *3 for faceting
colnames(W) <- str_c("PC",seq(1,dim(W)[2],1))
W <- as.data.frame(W)
W3 <- W[,c(1:3)]
W3 <- stack(list(PC1 = W3[,1], PC2 = W3[,2], PC3 = W3[,3]))
W3[,"Maturity"] <- rep(str_c(round(maturities.CFS,2),"Y"),3)
colnames(W3) <- c("LoadingValue", "Eigenvectors", "Maturity")

ggplot(W3, aes(x = Maturity, y = LoadingValue, colour = Eigenvectors, group = Eigenvectors)) + 
  geom_line() +
  scale_x_discrete(limits = str_c(round(maturities.CFS,2),"Y")) +
  #geom_point() +
  theme_light() +
  theme(legend.position="top") +
  ggtitle("Eigenvectors of the US daily spot rate Covariance Matrix") +
  ylim(c(-0.25,0.25)) +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle=90,size=rel(0.4)),
        plot.title = element_text(hjust = 0.5)) 
#ggsave("3_eigenvec_COV_int_perc_changes_SVENSSON.png",device="png",dpi=320,height=210,width=297,units="mm")

#### Variance Explained SVENSSON BAR PLOT Covariance ####
SvenLamdasCov <- eigen(CovXSvensson,symmetric=TRUE)$values
SUMSvenLamdasCov <- sum(SvenLamdasCov)
VarExplainedSvenCov <- as.data.frame((SvenLamdasCov/SUMSvenLamdasCov)*100)
VarExplainedSvenCov[,"PC"] <- str_c("PC",seq(1,dim(XSvensson)[2],1))
colnames(VarExplainedSvenCov)[1] <-  "VarExp"
VarExplainedSvenCov <- VarExplainedSvenCov %>% 
  mutate(VarExp = round(VarExp,2))

ggplot(VarExplainedSvenCov, aes(x=PC,y=VarExp)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=VarExp),vjust=-0.3) +
  theme_light() +
  scale_x_discrete(limits = str_c("PC",seq(1,dim(XSvensson)[2],1))[1:10],
                  labels = c(expression(lambda[1]),
                             expression(lambda[2]),
                             expression(lambda[3]),
                             expression(lambda[4]),
                             expression(lambda[5]),
                             expression(lambda[6]),
                             expression(lambda[7]),
                             expression(lambda[8]),
                             expression(lambda[9]),
                             expression(lambda[10]))) +
  scale_y_continuous(limits = c(0,100)) +
  theme(axis.text.x = element_text(size=rel(1.5)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  ggtitle("Variance explained (%) by first 10 PCs out of 226 - Svensson - Covariance Matrix") +
  ylab("Variance Exp. (%)") +
  xlab("PCs")
#ggsave("/Users/Andrea/Desktop/thesis/images/VAREXPLAINED_Svensson_Cov.png",device="png",dpi=320,height=210,width=297,units="mm")

#### PCA REPRESENTATION first three PCs COVARIANCE ####
Y3 <- Y[,1:3]
# Xhat contains in each column the vectors of the Absolute changes in basis points
Xhat <- Y3%*%t(W[,1:3])
colnames(Xhat) <- str_c("MAT",round(maturities.CFS,2))
Xhat <- xts(Xhat,index(X)[-1])
Mat_labs <- str_c("MAT",round(maturities.CFS,2))
plot(Xhat,
     main = "PCA Representation by means of first three PCs - Cov. Matrix - Absolute changes in basis points")
# compare with the original svensson:
plot(xts(XSvensson,index(X)[-1]),
     main = "Svensson absolute interest rates changes in Bps")

#### PCA on the CORRELATION matrix - Absolute rates changes in Bps ####
CorrXSvensson <- cov2cor(CovXSvensson)
W2 <- eigen(CorrXSvensson,symmetric = "TRUE")$vectors
J <- XSvensson%*%W2
J <- xts(J,index(X)[-1]) 
colnames(J) <- str_c("MAT",round(maturities.CFS,2))

#### Graph First three PCs Corr. Mat. - Absolute rates changes in Bps ####
highchart(type = "stock") %>% 
  hc_title(text = "First three PCs - Corr. Matrix - Svensson Absolute rates changes in Bps") %>% 
  hc_add_series(J[,1], name = "PC1") %>% 
  hc_add_series(J[,2], name = "PC2") %>%   
  hc_add_series(J[,3], name = "PC3") %>%
  hc_add_theme(hc_theme_flat()) %>% 
  hc_navigator(enabled = FALSE) %>% 
  hc_scrollbar(enabled = FALSE) %>% 
  hc_exporting(enabled = TRUE) %>% 
  hc_legend(enabled = TRUE) 

#### Graph of the Eigenvectors estiamted on CORRELLATION ####
# Absolute changes in basis points #
Wcr <- eigen(cor(XSvensson))$vectors
colnames(Wcr) <- str_c("PC",seq(1,dim(Wcr)[2],1)) 
Wcr <- as.data.frame(Wcr)
W3cr <- Wcr[,c(1:3)]
W3cr <- stack(list(PC1 = W3cr[,1], PC2 = W3cr[,2], PC3 = W3cr[,3]))
W3cr[,"Maturity"] <- rep(str_c(round(maturities.CFS,2),"Y"),3) # *3 labs for faceting
colnames(W3cr) <- c("LoadingValue", "Eigenvectors", "Maturity")

ggplot(W3cr, aes(x = Maturity, y = LoadingValue, colour = Eigenvectors, group = Eigenvectors)) + 
  geom_line() +
  scale_x_discrete(limits = str_c(round(maturities.CFS,2),"Y")) +
  #geom_point() +
  theme_light() +
  theme(legend.position="top") +
  ggtitle("Eigenvectors of the US daily spot rate Correlation Matrix") +
  ylim(c(-0.25,0.25)) +
  theme(panel.grid.major = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=90,size=rel(0.4))) 
#ggsave("/Users/Andrea/Desktop/thesis/images/3_eigenvec_CORR_int_perc_changes_SVENSSON.png",device="png",dpi=320,height=210,width=297,units="mm")

#### Variance Explained SVENSSON BAR PLOT Correlation ####
SvenLamdasCoR <- eigen(cov2cor(CovXSvensson),symmetric=TRUE)$values
SUMSvenLamdasCoR <- sum(SvenLamdasCoR)
VarExplainedSvenCoR <- as.data.frame((SvenLamdasCoR/SUMSvenLamdasCoR)*100)
VarExplainedSvenCoR[,"PC"] <- str_c("PC",seq(1,dim(XSvensson)[2],1))
colnames(VarExplainedSvenCoR)[1] <-  "VarExp"
VarExplainedSvenCoR <- VarExplainedSvenCoR %>% 
  mutate(VarExp = round(VarExp,2))

ggplot(VarExplainedSvenCoR, aes(x=PC,y=VarExp)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=VarExp),vjust=-0.3) +
  theme_light() +
  scale_x_discrete(limits = str_c("PC",seq(1,dim(XSvensson)[2],1))[1:10],
                   labels = c(expression(lambda[1]),
                              expression(lambda[2]),
                              expression(lambda[3]),
                              expression(lambda[4]),
                              expression(lambda[5]),
                              expression(lambda[6]),
                              expression(lambda[7]),
                              expression(lambda[8]),
                              expression(lambda[9]),
                              expression(lambda[10]))) +
  scale_y_continuous(limits = c(0,100)) +
  theme(axis.text.x = element_text(size=rel(1.5)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  ggtitle("Variance explained (%) by first 10 PCs out of 226 - Svensson - Correlation Matrix") +
  ylab("Variance Exp. (%)") +
  xlab("PCs")
#ggsave("/Users/Andrea/Desktop/thesis/images/VAREXPLAINED_Svensson_CoRR.png",device="png",dpi=320,height=210,width=297,units="mm")

#### PCA REPRESENTATION first three PCs CORRELATION ####
J3 <- J[,1:3]
# Xhat contains in each column the vectors of the Perc. int. rate changes at distinct mat.
Xhat2 <- J3%*%t(W2[,1:3])
colnames(Xhat2) <- Mat_labs
Xhat <- xts(Xhat2,index(X)[-1])
head(Xhat2)
Xhat2 <- xts(Xhat2,index(X)[-1])
plot(Xhat2,
     main = "PCA Representation by means of first three PCs - Corr Matrix - Absolute changes in basis points")

#########################
options(pillar.sigfig = 6)
t <- c(1,2,3,4)
cf <- c(6,6,6,106)
int <- c(4.5,4.75,4.85,5)
rm(d)
d <- as.data.frame(cbind(t,cf,int))
d <- d %>%
  as_tibble(d) %>% 
  mutate(pval=cf/(1+(int/100))^(t)) %>% 
  mutate(int01=int-0.01) %>% 
  mutate(pval01=cf/(1+int01/100)^(t)) %>% 
  mutate(bp01 = (((1+(int/100)-(0.01/100))^(-t))-(1+int/100)^(-t))) %>% 
  mutate(pv01=bp01*cf)
PV01 <- sum(d[,"pval01"])-sum(d[,"pval"])
##########################

#### Calculate Present Value of each Cash flow ####
rm(Int_CF)
options(scipen=999)
options(pillar.sigfig = 6)
t = dim(Svensson.ratesAll)[1]
Svensson.ratesAll[t]
Int_CF <- cbind(CFS[-1],maturities.CFS,as.vector(Svensson.ratesAll[t]),as.vector(Svensson.ratesAll[t-1]))
dim(Int_CF)
colnames(Int_CF) <- c("C","Maturity","R10032020","R09032020")
head(Int_CF)
Int_CF <- Int_CF %>% 
  as_tibble() %>% 
  mutate(Yield = R10032020/100) %>% 
  mutate(PVal = C/(1+Yield)^Maturity) %>% 
  mutate(Yield_PV01 = Yield-0.01/100) %>% 
  mutate(PVal001 = C/(1+Yield_PV01)^Maturity) %>% 
  mutate(BP01 = (((1+Yield-(0.01/100))^(-Maturity))-((1+Yield)^(-Maturity)))) %>% 
  mutate(PV01 = C*BP01) %>% 
  mutate(DeltaR = R10032020 - R09032020) %>% 
  select(Maturity, C, R10032020, R09032020, DeltaR, everything()) %>% 
  mutate(TdotPVal = Maturity*PVal) %>% 
  mutate(PV01dotDeltaR = PV01*DeltaR)

sum(Int_CF[,"PVal"])
sum(Int_CF[,"PVal001"])
MaculayDuration <- sum(Int_CF[,"TdotPVal"])/(100*26) #years
y <- jrvFinance::irr(cf = CFS, interval = c(0,1),cf.t=c(0,maturities.CFS)) #1.07893%
ModDuration <- MaculayDuration/(1+y)
# If the yield increase from 1.07% to 2.07% (100bps=0.010)
# the approximate percentage change in the value of the portfolio is:
-ModDuration*(0.010)*100 # -16.55%
# Dollar duration:
p <- abs(CFS[1])
-(ModDuration)*p
# Dollar price change
-(ModDuration)*p*(0.010)
# if the required yield would increase by 1.00% (100bps)
# then the estimated price change of the portfolio per $2,600.00 face value
# is -$574.80 

head(Int_CF)
tail(Int_CF)
apply(Int_CF,2,sum)
PV01 <- sum(Int_CF[,"PVal001"])-sum(Int_CF[,"PVal"])
PV01x <- sum(Int_CF[,"PV01"])
#IT TELLS US THAT IF THE ZERO RATES WERE TO SHIFT DOWN BY ONE BASIS
#POINT, THE Portfolio value WOULD RISE BY $4.31

x <- Int_CF %>% 
  filter(C > 15) %>% 
  mutate(R001=R10032020-0.01) %>% 
  filter(between(Maturity,1.45,1.48) | between(Maturity,3,3.05) | between(Maturity,5.5,5.7) | between(Maturity,24.10,24.30)) %>%
  mutate(Yield_PV01 = Yield_PV01*100) %>% 
  select(Maturity,C,R10032020,PVal,Yield_PV01,PVal001,PV01) 
  
apply(x,2,sum)[6]-apply(x,2,sum)[4]

# Principal component approximation for each interest rate changes
# Interest Rate changes
rm(SvenssonX)
Svensson.ratesAll <- xts(Svensson.ratesAll,index(X))
dim(Svensson.ratesAll)
Svensson.ratesAll[1,]
Svensson.ratesAll[3521,]
SvenssonX <- apply(Svensson.ratesAll,2,diff)
SvenssonX <- xts(SvenssonX,index(X)[-1])
SvenssonX[1,]
SvenssonX[3520,]
V <- cov(SvenssonX)
W <- eigen(V)$vectors
P <- SvenssonX%*%W
colnames(P) <- str_c("PC",seq(1,dim(P)[2]))
colnames(W) <- str_c("w",seq(1,dim(P)[2]))
rownames(W) <- Mat_labs
P[dim(P)[1],1:3]
W[,1:3]
# Rchanges contains the approximated Yield curves for each day from 2006 to 2020
# using only the first three principal components
Rchanges <- P[,1:3]%*%t(W[,1:3])
Rchanges <- xts(Rchanges,index(X)[-1])
Rchanges["2020-03-10"] # changes in the yield curve between 2020-03-10 and 2020-03-09
ApproxYieldChange10032020 <- t(Rchanges["2020-03-10"])
dim(ApproxYieldChange10032020)
ApproxYieldChange10032020 <- cbind(ApproxYieldChange10032020,maturities.CFS)
colnames(ApproxYieldChange10032020) <- c("YieldChange","Maturity")
ApproxYieldChange10032020 <- as.data.frame(ApproxYieldChange10032020)
class(ApproxYieldChange10032020)

ggplot(ApproxYieldChange10032020, aes(x=Maturity,y=YieldChange)) +
  geom_line() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Abs. Yield change", x="Maturity (Years)") +
  ggtitle("Principal component approximation of the the Yield Curve change\nbetween 2020-03-10 and 2020-03-09\nusing the first three PCs")
#ggsave("/Users/Andrea/Desktop/thesis/images/PCA_Yield_change_20200310.png",device="png",dpi=320,height=210,width=297,units="mm")

# Compare with:
SvenssonYieldChange <- cbind(t(SvenssonX[3520,]),maturities.CFS)
colnames(SvenssonYieldChange) <- c("YieldChange","Maturity")
SvenssonYieldChange <- as.data.frame(SvenssonYieldChange)

ggplot(SvenssonYieldChange,aes(x=Maturity,y=YieldChange)) +
  geom_line() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Abs. Yield Change", x="Maturity (Years)") +
  ggtitle("Svensson Yield Curve Change\nbetween 2020-03-10 and 2020-03-09")
#ggsave("/Users/Andrea/Desktop/thesis/images/SVENSSON_Yield_change_20200310.png",device="png",dpi=320,height=210,width=297,units="mm")

# Compare with:
rm(Xdiff)
Xdiff <- apply(X,2,diff)
Xdiff[3250,]
TrueYieldChange <- cbind(Xdiff[3250,],maturities.Tres)
colnames(TrueYieldChange) <- c("YieldChange","Maturity")
TrueYieldChange <- as.data.frame(TrueYieldChange)
ggplot(TrueYieldChange,aes(x=Maturity,y=YieldChange)) +
  geom_point() +
  geom_line() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Abs. Yield Change", x="Maturity (Years)") +
  ggtitle("True Yield Curve Change between 2020-03-10 and 2020-03-09")

#### Principal component Risk factor sensitivity ####
Int_CF <- Int_CF %>% 
  mutate(w1 = W[,1]) %>% 
  mutate(w2 = W[,2]) %>% 
  mutate(w3 = W[,3]) %>% 
  mutate(rf1 = PV01*w1) %>% 
  mutate(rf2 = PV01*w2) %>% 
  mutate(rf3 = PV01*w3)

# Risk Factor sensitivity
RF1 = -sum(Int_CF[,"rf1"])
RF2 = -sum(Int_CF[,"rf2"])
RF3 = -sum(Int_CF[,"rf3"])

# In alternative the risk factor RF sensitivity are given by:
# notice: that the risk factor sensitivity are constant!
# i.e. one vector w with three constant component
W3 <- t(W[,1:3])
PV01t <- as.matrix(as.data.frame(Int_CF[,"PV01"]))
RF <- W3%*%PV01t

#### Principal factor model representation of the Profit and loss dist of the portfolio ####
# for each day
PL <- -P[,1:3]%*%RF
colnames(PL) <- "PL"
PL <- xts(PL,index(X)[-1])
tail(PL)

PL_plot <- as.data.frame(PL)
head(PL_plot)
ggplot(PL_plot,aes(x=PL)) +
  geom_line(stat="density") +
  theme_light()

hist(PL)

#### EXCEL YTM ####
YTM_CALCULATED_EXCEL <- c(0.01511,0.00705,0.01716,0.02343,0.01860,0.02198,0.04491,0.04562,0.00493,0.04069,0.03955,0.03503,0.02830,
                          0.02949,0.02769,0.02798,0.02540,0.02693,0.02071,0.02167,0.02338,0.01987,0.01975,0.02002,0.02233,0.01838)
bonds[,"YTM_CALCULATED_EXCEL"] <- YTM_CALCULATED_EXCEL
























