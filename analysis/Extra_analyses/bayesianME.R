
pacman::p_load(ggplot2, dplyr, reshape2, stringr, tidyr, tidyverse, tibble, brms, cmdstanr,txtplot)

rm(list = ls()); gc()

source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")

# Simulate data ----
n_patients <- 500
n_timepoints <- 6
EID <- rep(1:n_patients, each = n_timepoints)
age <- rep(seq(20, 80, length.out = n_timepoints), n_patients)

# True parameters 
#+ Amax The maximum asymptotic value of the function. As age increases, the output approaches this value.
#+ A0: The minimum value (near the origin or lower bound). It helps determine the steepness and the scale of the curve relative to Amax.
#+ k:The growth rate or steepness of the curve. Higher k values make the sigmoid transition faster (i.e., sharper rise or fall around the inflection point).
#+ offset:A horizontal shift of the sigmoid curve. Controls the age at which the function rapidly increases—it’s essentially the inflection point.
#+

Amax_true <- 5
A0_true <- .001
k_true <- .2


set.seed(42)
offset_true <- rnorm(n_patients, mean = 60, sd = 20)
offset <- rep(offset_true, each = n_timepoints)

txtdensity(offset_true)

# Sigmoid-like function
ez_func <- function(age, Amax, A0, k, offset) {
  Amax / (1 + ((Amax - A0) / A0) * exp(-k * (age + offset)))
}

0:100 %>% 
  {x <- ez_func(., Amax_true,A0_true,k_true,0);
  plot(.,x,type = "l",lwd = 5, xlab = "age", ylab = "ezloss")
  }


# Simulate response with noise
set.seed(42)
ezloss <- ez_func(age, Amax_true, A0_true, k_true, offset) + rnorm(length(age),1, sd = 0.5) %>% 
  pmax(0)

# Combine into a data frame
data <- data.frame(ezloss, age, eid = factor(EID))

data$offset_true <- offset
data$true_offset_age <- data$offset_true+data$age


data %>% 
  pivot_longer(matches("age")) %>% 
  ggplot(aes(value, ezloss, color =  eid))+
  # geom_line(show.legend = FALSE)+
  # stat_smooth(aes( group =  eid),method = "loess", show.legend = FALSE, se = FALSE,   color =  rgb(0,0,0,.1))+
  stat_smooth(aes( group =  eid),method = "gam",formula = y~s(x, k = 4), show.legend = FALSE, se = FALSE,   color =  rgb(0,0,0,.1))+
  facet_wrap(~name, scales = "free")


# # Parameters from prior knowledge or model
# Amax_fixed <- 5
# A0_fixed <- 0.01
# k_fixed <- 0.2
# 
# # Estimate offset for each patient
# offsets <- data %>%
#   group_by(eid) %>%
#   group_map(~{
#     patient_data <- ungroup(.)
#     
#     loss_fn <- function(offset) {
#       pred <- ez_func(patient_data$age, Amax_fixed, A0_fixed, k_fixed, offset)
#       sum((patient_data$ezloss - pred)^2)
#     }
#     
#     opt <- optimize(loss_fn, interval = c(0, 100))  # Adjust bounds as needed
#     
#     opt$minimum
#   }) %>% 
#   unlist()
# 
# data$optim_offset <- rep(offsets, each = n_timepoints)
# data$optim_offset_age <-data$optim_offset+data$age


# model ez loss using sigmoid curve  ----
#+ You are modeling EZ.loss.Size (the outcome variable) as a sigmoid-shaped function of age, where the shape of the curve is governed by three global parameters:
#+   
#+ Amax: the maximum value the curve approaches as age increases
#+ 
#+ A0: the minimum value at early ages (near origin)
#+ 
#+ k: the steepness of the curve (how quickly it rises)
#+ 
#+ However, each individual (identified by EID) has their own personalized age shift
#+  — this is the offset — representing the age at which their sigmoid curve begins to rise.
#+ 
#+ Interpretation:  The model assumes that all individuals share the same maximum value, minimum value, and growth rate of EZ.loss.Size, 
#+ but differ in when the process starts — modeled by a random offset for each person.
#+ 

formula <- bf(ezloss ~ Amax  / (1 + ((Amax - A0) / A0) * exp(-k * (age+offset) )),
              offset ~ (1 | eid),
              Amax + A0 + k ~1,
              nl = TRUE)

prior <- prior(normal(5, 5), nlpar = 'Amax') +
  prior(normal(0, 0.2), nlpar = 'A0') +
  prior(normal(60, 10), nlpar = 'offset') +
  prior(uniform(0, 10), nlpar = 'k')

# Fit model
fit <- brm(
  formula = formula,
  data = data,
  prior = prior,
  chains = 4,
  cores = parallel::detectCores()-1,
  iter = 1000,
  backend = "cmdstanr"
)

# saveRDS(fit,"/var/folders/0l/1sf1trts12jgjwwysb4qrk5h0000gn/T//Rtmp91daPk/brmfitbe636cf4a29")
fit <- readRDS("/var/folders/0l/1sf1trts12jgjwwysb4qrk5h0000gn/T//Rtmp91daPk/brmfitbe636cf4a29")

plot(fit)
summary(fit)

offset_estimate <- coef(fit)$eid[,,"offset_Intercept"][, "Estimate"]

plot(offset_true, offset_estimate); add_lines()

txtdensity(offset_true)
txtdensity(offset_estimate)

# add to data
data$offset <- rep(offset_estimate, each = n_timepoints)
data$offset_age <- data$offset+data$age


data %>% 
  pivot_longer(c(age, true_offset_age,offset_age)) %>% 
  mutate(name = factor(name, levels = c("age", "true_offset_age", "offset_age"))) %>% 
  ggplot(aes(value, ezloss))+
  stat_smooth(
    aes(group = eid),
    method = "nls",
    formula = y ~ Amax_true / (1 + ((Amax_true - A0_true) / A0_true) * exp(-k_true * x)),
    method.args = list(start = list(Amax_true = Amax_true, A0_true = A0_true, k_true = k_true)),
    # formula = y ~ 5 / (1 + ((5 - .001) / .001) * exp(-.2 * x)),
    show.legend = FALSE,
    se = FALSE,
    color = rgb(0, 0, 0, 0.4)
  )+
  facet_wrap(~name, scales = "free")+
  theme_bw()

# Simulate from prior
set.seed(1)
mu <- rnorm(1000, mean = 7, sd = 1)
sigma <- abs(rnorm(1000, mean = 0, sd = 1))  # Half-normal

txtdensity(mu)
txtdensity(sigma)

# Simulate ratings from prior
prior_samples <- rnorm(5e3, mean = mu, sd = sigma)

# Keep within 0–10 (optional truncation)
prior_samples <- prior_samples[prior_samples >= 0 & prior_samples <= 10]

# Plot
ggplot(data.frame(x = prior_samples), aes(x)) +
  geom_density(fill = "lightblue") +
  labs(title = "Prior Predictive Distribution of Customer Experience",
       x = "Customer Rating", y = "Density")


# Load package
library(kernlab)

# create regression data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(length(x),sd=0.03)

# regression with gaussian processes
foo <- gausspr(x, y)
foo

# predict and plot
ytest <- predict(foo, x)
plot(x, y, type ="l")
lines(x, ytest, col="red")


#predict and variance
x = c(-4, -3, -2, -1,  0, 0.5, 1, 2)
y = c(-2,  0,  -0.5,1,  2, 1, 0, -1)
plot(x,y)
foo2 <- gausspr(x, y, variance.model = TRUE)
xtest <- seq(-4,2,0.2)
lines(xtest, predict(foo2, xtest))
lines(xtest,
      predict(foo2, xtest)+2*predict(foo2,xtest, type="sdeviation"),
      col="red")
lines(xtest,
      predict(foo2, xtest)-2*predict(foo2,xtest, type="sdeviation"),
      col="red")
