library(tidyverse)
library(lubridate)
library(broom)
library(jrvFinance)
library(YieldCurve)

# We load a data set of 26 US Gov Bonds retrieved from MarketInsider.com
# using a Python script.
# Then, we calculate the cash flows of each bond and yield to maturity
# and store them in a list separately.
# We spot some differences in the yield posted and the one calculated.
# The objective is to create one unique vector of cash flows which
# should be mapped to the KEY RATES (or vertices) over which PCA has been performed.

bonds <- read_csv("data/bonds-data.csv")

#### Bonds ####
bonds <- bonds %>% 
  # rename columns
  rename(MATURITY_DATE = MATURITY) %>% 
  rename(POSTED_YTM = YTM) %>% 
  # Parsing dates
  mutate(PURCHASE_DATE = "03-10-2020", across(contains("DATE"), mdy),
         DAYS_TO_COUPON = COUPON_PYMT_DATE - PURCHASE_DATE,
         # Semi-annual yield to maturity
         YTM2_POSTED = sqrt(1+POSTED_YTM)-1,
         # Semi-annual coupon payment amount
         SEMI_COUPON_AMOUNT = ((FACE_VAL*(CP_RATE/100))/2),
         NUMB_PAYMENTS = 2,
         # Calculate the number of coupons remaining for each bond
         N_CFs = coupons.n(settle = PURCHASE_DATE, mature = MATURITY_DATE, 
                           freq = 2))

#### Semi-annual cash flows ####
# First cash flow is negative, last = Face val + Coupon
cash_flows <- map(1:nrow(bonds), ~ c(-bonds[.x, c("MRKT_PRICE")], # t=t0
                       rep(bonds[.x, "SEMI_COUPON_AMOUNT"], bonds[.x, "N_CFs"]-1),
                       # t=tn
                       bonds[.x, "SEMI_COUPON_AMOUNT"] + bonds[.x, "FACE_VAL"])) %>% 
  set_names(nm = bonds$NAME)

yields <- enframe(cash_flows) %>% 
  unnest(cols = value) %>% 
  mutate(value = unlist(value)) %>% 
  rename(NAME = name) %>% 
  rename(CF = value) %>% 
  group_by(NAME) %>% 
  summarise(SEMI_ANN_YTM = irr(CF, interval = c(-1,1), cf.freq = 2, comp.freq = 2))

bonds <- left_join(bonds, yields, by = "NAME") %>% 
  mutate(ANNUAL_YTM = ((1 + SEMI_ANN_YTM)^2) - 1,
         DELTA_YTM = ANNUAL_YTM*100 - POSTED_YTM,
         # Create 7 groups of bonds that have common first payment date
         GROUP = factor(COUPON_PYMT_DATE, 
                        labels = c("I","II","III","IV","V","VI","VII")))

# Examples:
((1007.5/993)^2) - 1 # 0.02941766
jrvFinance::irr(c(-903.10, rep(100,3), 1000),
                interval = c(0,1), cf.freq = 1, comp.freq = 1)

#### Differences in YTM plot ####
ggplot(bonds,aes(x = NAME, y = DELTA_YTM)) +
  geom_col(colour = "grey60", fill = "lightblue") +
  theme_light() +
  theme(axis.text.x = element_text(angle=90),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank()) +
  ylab("Obtained YTM - Posted YTM") +
  xlab("Bond") +
  ggtitle("Differences between 'calculated' yearly YTM and 'posted' YTM on Market Insider.\n Quotations date: March 10th, 2020")
# ggsave("delta-annual-ytm.png", path = "plots", device = "png",
#        width = 297, height = 210, units = "mm")

#### Keyrates dates ####
# Starting from the purchase day 2020-03-10 we recover the actual dates 
# associated to the 11 different interest rates maturities 
s <- ymd("2020-03-10")
keyrates <- c(s + ddays(31), s + 3*ddays(31), s + 6*ddays(31), 
              s + dyears(1), s + 2*dyears(1), s + 3*dyears(1), 
              s + 5*dyears(1), s + 7*dyears(1), s + 10*dyears(1), 
              s + 20*dyears(1), s + 30*dyears(1)) %>% 
  as_tibble_col(column_name = "PAYMENT_DATE")

#### Bond portfolio cash flows ####
cf <- enframe(cash_flows) %>% 
  unnest(cols = value) %>% 
  mutate(value = unlist(value)) %>% 
  rename(NAME = name) %>% 
  rename(CF = value) %>% 
  left_join(bonds[, c("NAME", "MATURITY_DATE", "PURCHASE_DATE")], 
            by = "NAME") %>% 
  group_by(NAME) %>% 
  nest(CF = CF) %>% 
  mutate(across(.cols = contains("DATE"), as.character),
         # Find the dates of the cash flows
         PAYMENT_DATE = map2(MATURITY_DATE, PURCHASE_DATE, 
                              ~ jrvFinance::coupons.dates(settle = .y,
                                                          mature = .x,
                                                          freq = 2))) %>% 
  mutate(CF = map(CF, ~ slice(.x, -1))) %>% 
  unnest(c(CF, PAYMENT_DATE)) %>% 
  ungroup() %>% 
  group_by(PAYMENT_DATE) %>% 
  arrange(PAYMENT_DATE) %>% 
  summarise(CF = sum(CF)) %>% 
  ungroup() %>%
  full_join(keyrates, by = "PAYMENT_DATE") %>% 
  arrange(PAYMENT_DATE)

# Graph portfolio cash flows
cf %>% 
  mutate(PAYMENT_DATE = factor(PAYMENT_DATE)) %>% 
  ggplot(aes(x = PAYMENT_DATE, y = CF)) +
  geom_col(colour = "grey60", fill = "lightblue") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, size = rel(0.6)),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(title = "Portfolio cash flows with keyrates: 1m, 3m, 6m, 1y, 2y, 3y, 5y, 7y, 10y, 20y, 30y",
  y = "($)",
  x = "Day") +
  geom_vline(xintercept = which(is.na(cf$CF)), 
             linetype = 1, 
             colour = 'grey', size = 0.4)
# ggsave("portfolio-cash-flows.png", path = "plots", units = "mm",device = "png", 
#        width = 297, height = 210, dpi = 320)

#### Sample spot yield curve ####
rates <- read_csv("data/USTREASURY.csv") %>% 
  filter(DATE <= "2020-03-10")

rates %>% 
  filter(DATE == "2020-03-10") %>% 
  pivot_longer(!DATE, names_to = "MATURITY", values_to = "RATE") %>% 
  mutate(X = c(1/12,3/12,6/12,1,2,3,5,7,10,20,30)) %>% 
  ggplot(aes(x = X, y = RATE)) +
  geom_line() +
  geom_point() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "Yield (%)", x = "Maturity (Years)",
       title = "Spot Rate US Treasury Curve on March 10, 2020")
# ggsave("sample-yield-curve.png", path = "plots", 
#        device = "png", units = "mm", height = 210, width = 297, dpi = 320)

#### Svensson Parametric model ####
# loss of accuracy a short and long terms
# less flexible than splines
# In lenghts_from_t0 we store the number of days that separates 
# each cash flow to 2020-03-10

#### Svensson ####
maturities <- c(1/12, 3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30)
t0 <- ymd("2020-03-10")

maturitiesCF <- cf %>% 
  # Remove the keyrates 
  filter(!is.na(CF)) %>% 
  pull(PAYMENT_DATE)

# Payments maturities as integers
matSven <- as.numeric((maturitiesCF - t0)/360)
#ratesforSvensson <- rates %>% select(-1) %>% as.matrix() 
#svenssonPars <- Svensson(rate = ratesforSvensson, maturity = maturities)
svenssonPars <- read.csv("data/svensson-pars.csv")
rownames(svenssonPars) <- rates$DATE
#write_csv(as_tibble(svenssonPars), path = "data/svensson-pars.csv")
svenssonRates <- Srates(Coeff = xts(svenssonPars, rates$DATE), 
                        maturity = matSven, whichRate = "Spot")
colnames(svenssonRates) <-  map2_chr(rep("MAT", length(matSven)), round(matSven, 2), 
                                     ~ str_c(.x, .y))
#### Sample Svensson curve ####
# Estimate of the Svensson parameters for March 10, 2020 
svenssonRates %>% 
  as_tibble() %>% 
  slice_tail() %>% 
  pivot_longer(cols = everything(), names_to = "MATURITY", 
               values_to = "SVN_RATE") %>% 
  bind_cols(matSven) %>% 
  ggplot(aes(x = ...3, y = SVN_RATE)) +
  geom_line() +
  theme_light() +
  theme(plot.title = element_text(hjust= 0.5)) +
  labs(x = "Maturity (Years)", y = "Yield (%)",
       title = "Estimated 'Svensson' Spot Curve on March 10, 2020")
# ggsave("sample-svensson-curve.png", path = "plots", device = "png", units = "mm",
#        height = 210, width = 297, dpi = 320)

# Svensson Term structure
sample_svensson_yield <- svenssonRates %>% 
  as_tibble() %>% 
  mutate(DATE = rates$DATE) %>% 
  slice_tail() %>% 
  select(DATE) %>% 
  pull()

svenssonRates %>% 
  as_tibble() %>% 
  mutate(DATE = rates$DATE) %>% 
  relocate(DATE, .before = everything()) %>% 
  pivot_longer(cols = !DATE, names_to = "MATURITY", values_to = "RATE") %>% 
  ggplot(aes(x = DATE, y = RATE, group = MATURITY, colour = MATURITY)) +
  geom_line(alpha = 0.3) +
  geom_vline(xintercept = sample_svensson_yield, colour = "blue") +
  theme_light() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Recovered Term Structure by means of Svensson's model",
    x = "Year",
    y = "Rate")
# ggsave("svensson-term-structure.png", path = "plots", device = "png", units = "mm",
#        height = 210, width = 297, dpi = 320)

# Svensson: Absolute Interest Rates changes in Bps
svenAbsIntRatesChange <- as_tibble(svenssonRates) %>% 
  mutate(across(where(is.numeric), ~ (.x - lag(.x, n = 1))*100)) %>% 
  mutate(DATE = rates$DATE) %>% 
  slice(-1) %>% 
  relocate(DATE, .before = everything())

#### PCA on Svensson Yield Curves ####
# PCA - Absolute rates changes in Bps
# prcomp uses SVD

pca <- svenAbsIntRatesChange %>% 
  mutate(across(.cols = where(is.numeric), .fns = ~ (. - mean(.))/sd(.))) %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = FALSE)

# Eigenvalues
pca %>% 
  tidy(matrix = "eigenvalues") %>% 
  slice(1:10) %>% 
  ggplot(aes(x = PC, y = percent)) +
  geom_col(fill = "lightblue", colour = "grey", alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PCs",
       y = "(%)",
       title = "Variance explained (%) by first 10 PCs out of 226 - Svensson - Covariance Matrix")
# ggsave("variance-explained-svensson.png", path = "plots", device = "png", dpi = 320,
#        height = 210, width = 297, units = "mm")

# Eigenvectors
pca %>% 
  tidy(matrix = "loadings") %>% 
  arrange(PC) %>%
  filter(PC %in% c(1,2,3)) %>% 
  mutate(column = rep(colnames(svenssonRates), 3),
         PC = factor(PC)) %>% 
  filter(PC %in% c(1, 2, 3)) %>% 
  ggplot(aes(x = column, y = value, colour = PC, group = PC)) +
  geom_line() +
  labs(x = "Maturity",
       y = "Loadings",
       title = "Eigenvectors of the US daily spot rate Covariance Matrix") +
  scale_x_discrete(limits = colnames(svenssonRates)) +
  theme_light() +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, size = rel(0.4)),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  ylim(c(-0.25,0.25))
# ggsave("three-eigenvectors-svensson.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units="mm")

# Principal components
pca %>% 
  tidy(matrix = "scores") %>% 
  filter(PC %in% c(1,2,3)) %>% 
  arrange(PC) %>% 
  mutate(row = rep(rates$DATE[-1], 3)) %>% 
  rename(DATE = row) %>% 
  mutate(PC = factor(PC)) %>% 
  ggplot(aes(x = DATE, y = value, fill = PC, colour = PC)) +
  geom_line(alpha = 0.6) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "First three PCs - Svensson absolute rates changes in Bps",
       x = "Years",
       y = "Value")
# ggsave("three-pcs-svensson.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units = "mm")

#### PCA Representation first three PCs  ####
pcsSven <- pca %>% 
  tidy(matrix = "scores") %>% 
  filter(PC %in% c(1,2,3)) %>% 
  arrange(PC) %>% 
  mutate(row = rep(rates$DATE[-1], 3)) %>% 
  rename(DATE = row) %>% 
  pivot_wider(names_from = PC, names_prefix = "PC", values_from = value) %>% 
  select(where(is.numeric)) %>% 
  as.data.frame()

eigenvalSven <- pca %>% 
  tidy(matrix = "loadings") %>% 
  arrange(PC) %>%
  filter(PC %in% c(1,2,3)) %>% 
  mutate(column = rep(colnames(svenssonRates), 3)) %>% 
  pivot_wider(names_from = PC, names_prefix = "PC", values_from = value) %>% 
  select(where(is.numeric)) %>% 
  as.data.frame()

as.matrix(pcsSven)%*%t(as.matrix(eigenvalSven)) %>% 
  as_tibble(.name_repair = ~ colnames(svenssonRates)) %>% 
  mutate(DATE = rates$DATE[-1]) %>% 
  relocate(DATE, .before = everything()) %>% 
  pivot_longer(cols = !DATE, names_to = "MATURITY", values_to = "RATE") %>% 
  ggplot(aes(x = DATE, y = RATE, group = MATURITY, colour = MATURITY)) +
  geom_line(alpha = 0.2) +
  theme_light() +
  guides(colour = F) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Principal component representation Svensson",
       x = "Years", y = "Rate")

#### Calculate Present Value of each Cash flow ####
R10032020 <- as_tibble(svenssonRates) %>% 
  slice_tail() %>% 
  pivot_longer(cols = contains("MAT")) %>% 
  pull()

t <- nrow(svenssonRates)
R09032020 <- as_tibble(svenssonRates) %>% 
  slice(t-1) %>% 
  pivot_longer(cols = contains("MAT")) %>% 
  pull()

cf <- cf %>% 
  filter(!is.na(CF)) %>% 
  bind_cols(matSven, .name_repair = "universal") %>% 
  rename(Maturity = ...3) %>% 
  rename(C = CF) %>% 
  mutate(R10032020 = R10032020,
         R09032020 = R09032020,
         Yield = R10032020/100,
         PVal = C/(1+Yield)^Maturity,
         Yield_PV01 = Yield-0.01/100,
         PVal001 = C/(1+Yield_PV01)^Maturity,
         BP01 = (((1+Yield-(0.01/100))^(-Maturity))-((1+Yield)^(-Maturity))),
         PV01 = C*BP01,
         DeltaR = R10032020 - R09032020) %>% 
  relocate(Maturity, C, R10032020, R09032020, DeltaR, .before = everything()) %>% 
  mutate(TdotPVal = Maturity*PVal,
         PV01dotDeltaR = PV01*DeltaR)

PV01 <- sum(cf[,"PVal001"]) - sum(cf[,"PVal"])
PV01x <- sum(cf[,"PV01"])
#IT TELLS US THAT IF THE ZERO RATES WERE TO SHIFT DOWN BY ONE BASIS
#POINT, THE Portfolio value WOULD RISE BY $4.257581

# Principal component approximation for each interest rate changes
pcStats <- svenssonRates %>% 
  as_tibble() %>% 
  mutate(across(everything(), ~ (.x - lag(.x, n = 1)))) %>% 
  pivot_longer(cols = everything(), names_to = "MATURITY", values_to = "RATES") %>% 
  group_by(MATURITY) %>% 
  summarise(mean = mean(RATES, na.rm = TRUE), sd = sd(RATES, na.rm = TRUE)) %>% 
  arrange(MATURITY)

pcaRiskFactorModel <- svenssonRates %>% 
  as_tibble() %>% 
  mutate(across(everything(), ~ (.x - lag(.x, n = 1))),
         DATE = ymd(index(svenssonRates))) %>% 
  relocate(DATE, .before = everything()) %>% 
  slice(-1) %>% 
  mutate(across(where(is.numeric), ~ (.x - mean(.x))/(sd(.x)))) %>% 
  select(!DATE) %>% 
  prcomp(scale = FALSE)

pcsSvenRFM <- pcaRiskFactorModel %>% 
  tidy(matrix = "scores") %>% 
  filter(PC %in% c(1,2,3)) %>% 
  arrange(PC) %>% 
  mutate(row = rep(rates$DATE[-1], 3)) %>% 
  pivot_wider(names_from = PC, names_prefix = "PC", values_from = value) %>% 
  rename(DATE = row)

eigenvecSvenRFM <- pcaRiskFactorModel %>% 
  tidy(matrix = "loadings") %>% 
  arrange(PC) %>% 
  filter(PC %in% c(1,2,3)) %>% 
  mutate(column = rep(colnames(svenssonRates), 3)) %>% 
  pivot_wider(names_from = PC, names_prefix = "w", values_from = value) %>% 
  rename(MATURITY = column)

means <- pcStats %>% 
  select(mean) %>% 
  pull() %>% 
  as.list() %>% 
  set_names(nm = pcStats$MATURITY)

sd <- pcStats %>% 
  select(sd) %>% 
  pull() %>% 
  as.list() %>% 
  set_names(nm = pcStats$MATURITY)

# Graph Svensson yield change on March 10, 2020 approx. by first three PCs
as.matrix(pcsSvenRFM[,2:4])%*%t(as.matrix(eigenvecSvenRFM[,2:4])) %>% 
  as_tibble(.name_repair = ~ colnames(svenssonRates)) %>% 
  mutate(DATE = ymd(index(svenssonRates)[-1])) %>% 
  relocate("DATE", .before = everything()) %>% 
  # Add the mean and multiply by sd
  mutate(across(all_of(names(means)),
                ~ (.x + means[[cur_column()]]))) %>% 
  mutate(across(all_of(names(sd)), 
                ~ (.x * sd[[cur_column()]]))) %>% 
  # Changes in the Svensson yield curve between 2020-03-10 and 2020-03-09
  # which corresponds to the settling day of the bond portfolio.
  slice_tail() %>% 
  pivot_longer(cols = !DATE, names_to = "Maturity", values_to = "Rate") %>% 
  select(-DATE) %>% 
  mutate(Maturity = str_remove_all(Maturity, "MAT"),
         Maturity = as.numeric(Maturity)) %>% 
  ggplot(aes(x = Maturity, y = Rate)) +
  geom_line() + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "Abs. Yield change", x = "Maturity (Years)") +
  ggtitle("Principal component approximation of the the Yield Curve change\nbetween 2020-03-10 and 2020-03-09\nusing the first three PCs")
# ggsave("pca-svensson-yield-change-20200310.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width=297, units = "mm")

#### Principal component Risk factor sensitivity ####
cf <- cf %>% 
  mutate(w1 = eigenvecSvenRFM$w1,
         w2 = eigenvecSvenRFM$w2,
         w3 = eigenvecSvenRFM$w3) %>% 
  mutate(rf1 = PV01*w1,
         rf2 = PV01*w2,
         rf3 = PV01*w3)

RF1 = -sum(cf[,"rf1"])
RF2 = -sum(cf[,"rf2"])
RF3 = -sum(cf[,"rf3"])

# In alternative, the risk factor RF sensitivity are given by:
# notice: that the risk factor sensitivity are constant!
# i.e. one vector w with three constant component
riskFactors <- t(as.matrix(eigenvecSvenRFM[,c("w1","w2","w3")]))%*%as.matrix(cf$PV01)

#### Principal factor model rep. of the Profit and Loss of the bond port. ####
# for each day
-(as.matrix(pcsSven)%*%as.matrix(riskFactors)) %>% 
  as_tibble() %>% 
  ggplot(aes(x = V1)) +
  geom_histogram(colour = "grey60", fill = "lightblue", bins = 150) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Principal component factor model representation\nProfit and Loss of the bond portfolio at time t = '2020-03-10'",
       x = "")

