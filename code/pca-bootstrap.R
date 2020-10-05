# PACKAGES
library(tidyverse)
library(lubridate)
library(broom)
library(lattice)

#### US TREASURY ####
mat <- c("MAT1MO","MAT3MO","MAT6MO","MAT1YR","MAT2YR",
         "MAT3YR","MAT5YR","MAT7YR","MAT10YR","MAT20YR","MAT30YR")

rates <- read_csv("data/USTREASURY.csv") %>% 
  filter(DATE <= "2020-03-10") %>% 
  mutate(DATE = ymd(DATE))

rates %>% 
  pivot_longer(cols = where(is.numeric), 
               names_to = "Maturity", values_to = "RATE") %>% 
  mutate(Maturity = factor(Maturity,
                           levels = mat)) %>%
  ggplot(aes(x = DATE, y = RATE, group = Maturity, colour = Maturity)) +
  geom_line(alpha = 0.7) +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "US Treasury Interest Rates",
       x = "Year",
       y = "Rate")
# ggsave("us-treasury-yields.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units = "mm")

## INTEREST RATES CHANGES IN BASIS POINTS ##
ratesChangesBps <- rates %>% 
  mutate(across(where(is.numeric), ~ (.x - lag(.x, n = 1))*100))

## CENTER and SCALE, SINCE SVD IS USED IN PRCOMP ##
pca <- ratesChangesBps %>% 
  select(where(is.numeric)) %>% 
  slice(2:nrow(rates)) %>% 
  mutate(across(.cols = everything(), .fns = ~ (. - mean(.))/sd(.))) %>% 
  prcomp()

# Check centering and scaling
ratesChangesBps %>% 
  select(where(is.numeric)) %>% 
  slice(2:nrow(rates)) %>% 
  mutate(across(.cols = everything(), .fns = ~ (. - mean(.))/sd(.)),
         DATE = rates$DATE[-1]) %>% 
  pivot_longer(cols = !DATE, names_to = "MATURITY", values_to = "RATES") %>%
  select(!DATE) %>% 
  group_by(MATURITY) %>% 
  summarise(count = n(), mean = round(mean(RATES), 3), sd = round(sd(RATES), 2))

#### EIGENVALUES ####
pca %>% 
  tidy(matrix = "eigenvalues") %>% 
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
       title = "Variance explained (%) by the 11 PCs")
# ggsave("variance-explained-eigenvalues.png", path = "plots", device = "png", 
#        dpi = 320, height = 210, width = 297, units = "mm")

# GRAPH FIRST 3 PCS
pca %>% 
  tidy(matrix = "scores") %>% 
  arrange(PC) %>% 
  mutate(row = rep(rates$DATE[-1], 11),
         PC = factor(PC)) %>% 
  rename(DATE = row) %>% 
  filter(PC %in% c(1,2,3)) %>% 
  ggplot(aes(x = DATE, y = value, fill = PC, colour = PC)) +
  geom_line(alpha = 0.5) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "First three PCs - Absolute rates changes in Bps",
       x = "Years",
       y = "Value")
# ggsave("three-pcs.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units = "mm")

#### EIGENVECTORS ####
pca %>% 
  tidy(matrix = "loadings") %>% 
  arrange(PC) %>% 
  mutate(column = rep(mat, 11),
         PC = factor(PC)) %>% 
  rename(MATURITY = column) %>% 
  filter(PC %in% c(1, 2, 3)) %>% 
  ggplot(aes(x = MATURITY, y = value, colour = PC, group = PC)) +
  geom_line() +
  labs(x = "Maturity",
       y = "Loadings",
       title = "Eigenvectors from SVD of the US daily spot rate") +
  scale_x_discrete(limits = mat) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
# ggsave("three-eigenvectors.png", path = "plots",
#        device = "png", dpi = 320, height = 210, width = 297, units="mm")

#### PCA REPRESENTATION ####

# First three eigenvectors
eigenvectors <- pca %>% 
  tidy(matrix = "loadings") %>% 
  arrange(PC) %>% 
  mutate(column = rep(colnames(rates)[-1], 11)) %>% 
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>% 
  select(PC1, PC2, PC3)

# First three principal components
prComponents <- pca %>% 
  tidy(matrix = "scores") %>% 
  arrange(PC) %>% 
  mutate(row = rep(rates$DATE[-1], 11)) %>% 
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>% 
  select(PC1, PC2, PC3)

# Principal component representation using the first three components
prRepresentation <- as.matrix(prComponents)%*%t(as.matrix(eigenvectors)) %>% 
  as_tibble(.name_repair = ~ mat) %>% 
  mutate(DATE = rates$DATE[-1]) %>% 
  relocate(DATE, .before = everything())

prRepresentation %>% 
  pivot_longer(cols = !DATE, names_to = "MATURITY", values_to = "RATE") %>% 
  mutate(MATURITY = factor(MATURITY, 
                           levels = mat)) %>% 
  ggplot(aes(x = DATE, y = RATE, colour = MATURITY, group = MATURITY)) +
  geom_line(alpha = 0.5) +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Principal component representation using the first three PCs",
       x = "Year",
       y = "Rate")
# ggsave("pca-representation.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units="mm")

ratesChangesBpsScaled <- ratesChangesBps %>% 
  select(where(is.numeric)) %>% 
  slice(2:nrow(rates)) %>% 
  mutate(across(.cols = everything(), .fns = ~ (. - mean(.))/sd(.)),
         DATE = rates$DATE[-1]) %>% 
  relocate(DATE, .before = everything()) 

ratesChangesBpsScaled %>% 
  pivot_longer(cols = !DATE, names_to = "MATURITY", values_to = "RATE") %>% 
  ggplot(aes(x = DATE, y = RATE, colour = MATURITY, group = MATURITY)) +
  geom_line(alpha = 0.5) +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "US Treasury Interest rates centered and scaled",
       x = "Year",
       y = "Rate")
# ggsave("us-treasury-centeres-scaled.png", path = "plots", 
#        device = "png", dpi = 320, height = 210, width = 297, units="mm")

#### COMPARISON CORRELATION MATRICES ####

# Correlation matrix representation
corrRepres <- prRepresentation %>% 
  select(where(is.numeric)) %>% 
  cor()

# Correlation matrix 
corrActual <- ratesChangesBpsScaled %>% 
  select(where(is.numeric)) %>% 
  cor()

round(abs(corrRepres - corrActual), 2)

#### PCA ON WINDOWS ####
# We investigate the persistence of the fundamental structure:
# parallel shift, tilt, curvature
us_rates <- read_csv("data/USTREASURY.csv")
window <- (nrow(us_rates) - 2)/6

# prcomp function on numeric variables
principalComp <- function(df) {
  df = select(df, where(is.numeric))
  pca = prcomp(df, scale = FALSE)
}

pca_windows <- us_rates %>% 
  mutate(across(where(is.numeric), ~ (.x - lag(.x, n = 1))*100)) %>% 
  slice(-c(1,2)) %>% 
  mutate(across(.cols = everything(), .fns = ~ (. - mean(.))/sd(.)),
         DATE = us_rates$DATE[-c(1,2)]) %>% 
  mutate(WINDOW = unlist(map(c("A","B","C","D","E","F"), ~ rep(., window)))) %>% 
  group_by(WINDOW) %>% 
  nest() %>% 
  mutate(pca = map(data, principalComp),
         eigenvalues = map(pca, ~ tidy(.x, matrix = "eigenvalues")),
         eigenvectors = map(pca, ~ tidy(.x, matrix = "loadings")),
         pcomponents = map(pca, ~ tidy(.x, matrix = "scores")))
  
# PLOTS WINDOWS
pca_windows <- pca_windows %>% 
  mutate(plt_eigenval = map2(eigenvalues, data, ~ ggplot(.x, aes(.x$PC, .x$percent)) +
                            geom_col(fill = "lightblue", colour = "grey", alpha = 0.8) +
                            scale_x_continuous(breaks = 1:11) +
                            scale_y_continuous(labels = scales::percent_format()) +
                            theme_light() +
                            theme(panel.grid.major.x = element_blank(),
                                  panel.grid.minor.x = element_blank(),
                                  plot.title = element_text(hjust = 0.5),
                                  plot.subtitle = element_text(hjust = 0.5)) +
                            labs(x = "PCs", y = "(%)",
                                 title = "Variance explained (%) by the 11 PCs",
                                 subtitle = str_c("from", 
                                                  as.character(.y$DATE[1]),
                                                  "to",
                                                  as.character(.y$DATE[nrow(.y)]),
                                                  sep = " "))),
         plt_eigenvec = map2(eigenvectors, data, ~ dplyr::arrange(.x, .x$PC) %>% 
                               mutate(column = rep(unique(column), 11),
                                      PC = factor(PC)) %>% 
                               filter(PC %in% c(1,2,3)) %>% 
                               ggplot(aes(x = column, y = value, 
                                          colour = PC, group = PC)) +
                               geom_line() +
                               labs(x = "Maturity",
                                    y = "Loadings",
                                    title = "Eigenvectors of the US daily spot rate (SVD)",
                                    subtitle = str_c("from",
                                                     as.character(.y$DATE[1]),
                                                     "to",
                                                     as.character(.y$DATE[nrow(.y)]),
                                                     sep = " ")) +
                               scale_x_discrete(limits = unique(.x$column)) +
                               theme_light() +
                               theme(panel.grid.major = element_blank(),
                                     plot.title = element_text(hjust = 0.5),
                                     plot.subtitle = element_text(hjust = 0.5),
                                     legend.position = "top")))

#### BOOTSTRAP ####
# Start from data TS: Absolute changes in basis points
# We create 10.000 bootstrap samples each made of 582 obs.
# Ratio 1/6 = 582/3520

boot_rates <- us_rates %>% 
  mutate(across(where(is.numeric), ~ (.x - lag(.x, n = 1))*100)) %>% 
  slice(-1) %>% 
  mutate(across(.cols = everything(), .fns = ~ (. - mean(.))/sd(.)),
         DATE = us_rates$DATE[-1]) %>% 
  filter(DATE <= "2020-03-10")

set.seed(1234)
reps <- 10000
boot_samples <- 
  tibble(sample_id = 1:reps,
         samples = map(1:reps, ~ slice_sample(boot_rates, prop = 1/6, 
                                              replace = T))) %>% 
  group_by(sample_id) 

pca_boot <- boot_samples %>% 
  mutate(pca = map(samples, principalComp),
         eigenvalues = map(pca, ~ tidy(.x, matrix = "eigenvalues")),
         eigenvectors = map(pca, ~ tidy(.x, matrix = "loadings")))

##### BOOT: DENSITIES EIGENVECTORS ####
pca_boot %>% 
  select(eigenvalues) %>% 
  unnest(cols = "eigenvalues") %>% 
  mutate(PC = as.character(PC)) %>% 
  filter(PC %in% c(1,2,3)) %>% 
  ggplot(aes(x = std.dev, y = ..density.., fill = PC)) +
  geom_histogram(bins = 100) +
  geom_density(size = 0.5, colour = "grey60") +
  facet_grid(PC ~ .) +
  theme_light() +
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.position="top",
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank()) +
  labs(title = "Eigenvalues Bootstrap Estimates",
       x = expression(lambda))
# ggsave("densities-bootstrap-eigenvalues.png", path = "plots",
#        device = "png", dpi = 320, height = 210, width = 210, units="mm")

#### BOOT: BOX PLOT EIGENVECTORS ####
pca_boot %>% 
  ungroup() %>% 
  select(eigenvalues) %>% 
  unnest(cols = "eigenvalues") %>% 
  mutate(PC = factor(PC)) %>% 
  filter(PC %in% c(1,2,3)) %>% 
  ggplot(aes(x = PC, y = std.dev)) + 
  geom_boxplot(outlier.shape = 21) + 
  theme_light() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c(expression(lambda[1]),
                              expression(lambda[2]),
                              expression(lambda[3]))) +
  labs(title = "Eigenvalues Bootstrap Estimates, Num of reps. = 10.000, Prop = 1/6")
# ggsave("box-plots-bootstrap-eigenvalues.png", path = "plots",
#        device = "png", dpi = 320, height = 210, width = 210, units="mm")

#### BOOT: DENSITY VARIANCE EXPLAINED ####
# TABLE
pca_boot %>% 
  ungroup() %>% 
  select(eigenvalues) %>% 
  unnest(cols = "eigenvalues") %>% 
  select(PC, cumulative) %>% 
  filter(PC == 3) %>% 
  summary()

# PLOT
pca_boot %>% 
  ungroup() %>% 
  select(eigenvalues) %>% 
  unnest(cols = "eigenvalues") %>% 
  select(PC, cumulative) %>% 
  filter(PC == 3) %>% 
  ggplot(aes(x = cumulative, y = ..density..)) + 
  geom_histogram(bins = 50, colour = "grey60", fill = "lightblue") +
  geom_density() +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Variance Explained by the first three PCs on 10.000 Bootstrap samples",
       x = expression((lambda[1]+lambda[2]+lambda[3])/(lambda[1]+lambda[2]+...+lambda[11])))
# ggsave("density-var-explained.png", path = "plots",device = "png", dpi = 320, 
#        height = 210, width = 210, units = "mm")


##### BOOT: LOADINGS ####
pca_boot %>% 
  ungroup() %>% 
  select(eigenvectors) %>% 
  unnest(cols = "eigenvectors") %>% 
  arrange(PC) %>% 
  mutate(column = rep(mat, 11*reps)) %>% 
  filter(PC %in% c(1, 2, 3)) %>% 
  mutate(PC = factor(PC, labels = c("w1", "w2", "w3"))) %>% 
  ggplot(aes(x = value, fill = PC, group = PC)) + 
  geom_histogram(bins = 300) + 
  facet_grid(column ~ PC) +
  theme_light() +
  guides(fill = FALSE) +
  scale_y_continuous(breaks = c(0, 500, 1000)) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eeeeee", colour = "grey"),
        strip.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Value of the Loadings", y = "Frequency",
       title = "Histogram of the estimated first three eigenvectors loadings\n obtained from 10000 boostrap samples")
# ggsave("bootstrap-loadings.png", path = "plots",
#        device = "png", dpi = 320, height=297, width=210, units = "mm")
