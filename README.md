# Folded Normal RNG

A multivariate folded normal distribution is a distribution of the absolute value of a Gaussian random vector. In thie repo, we provide several instances to generate multivariate folded normal random numbers.

## Test for 2 dimensions Gibbs Method with naive R

```R
## Indicator Function
Indicator <- function(x){
  return(ifelse(FALSE %in% ifelse(x >= 0, TRUE, FALSE), 0, 1))
}

## Folded Normal Sampling Funcion
FN_Sampling <- function(number, mu, sig){
  sample <- abs(rnorm(number, mean = mu, sd = sqrt(sig)))
  density <- function(x, u, s){
    1/sqrt(2 * pi)/sqrt(s) * (exp(-(x - u)^2/(2 * s)) + exp(-(x + u)^2/(2 * s)) * Indicator(x))
  }
  return(list(sample,density(sample, mu, sig)))
}

## Gibbs Sampling Function
Gibbs <- function(N,x_init,sig1,sig2,rho){
  x_seq <- matrix(0, N, 2)
  x_seq[1, ] <- x_init
  for (i in 2:N){
    ## x1|x2
    x_seq2 <- x_seq[i - 1, 2]
    x_seq[i, 1] <- FN_Sampling(1, x_seq2 * rho * sig1/sig2, sig1 * (1 - rho^2))[[1]]
    ## x2|x1
    x_seq1 <- x_seq[i, 1]
    x_seq[i, 2] <- FN_Sampling(1, x_seq1 * rho * sig2/sig1, sig2 * (1 - rho^2))[[1]]
  }
  return(x_seq)
}

## Input parameters
N <- 10000 
Chain_n <- 4
sig1 <- 1
sig2 <- 1
rho <- 0.5
x_initial <- list(first = 0, second = 0, third = 0, fourth = 0)
x_init <- x_initial[[1]] <-matrix(c(0.5, 0.5), ncol = 2)
x_initial[[2]] <- matrix(c(0.5, 1.5), ncol = 2)
x_initial[[3]] <- matrix(c(1.5, 0.5), ncol = 2)
x_initial[[4]] <- matrix(c(1.5, 1.5), ncol = 2)

## Simulation 
set.seed(1010)
sample_p1 <- sample_p2 <- 0
Chain_all <- list(Chain1 = 0, Chain2 = 0, Chain3 = 0, Chain4 = 0)
t_start <- Sys.time()
for(j in 1:Chain_n){
  result = Gibbs(N,x_init,sig1,sig2,rho)
  Chain_all[[j]] <- result
  if(j == 1){
    sample_p1 <- result[,1]
    sample_p2 <- result[,2]
  }
  else{
    sample_p1 <- c(sample_p1, result[,1])
    sample_p2 <- c(sample_p2, result[,2])
  }
  if(j != Chain_n){
    x_init = x_initial[[j+1]]
  }
}
t_end <- Sys.time()
t_end - t_start
#> Time difference of 0.3977759 secs

## Data Clean
dim2 <- paste0("chain:", 1:4)
dim3 <- paste0("theta", 1:2)
sample_array <- array(c(sample_p1,sample_p2),c(N,4,2),dimnames = list(NULL,dim2,dim3))
names(dimnames(sample_array)) <- c("Iteration", "Chain", "Parameter")
```

We notice the runningtime is relatively long (0.3977759 secs). Therefore, the `rfnorm` package with a C++ implement is provided (under development).

| Desc | logs |
|:-----|:-----|
| Package | rfnorm |
| maintainer  | Yifan Yang <mailto: yfyang.86@hotmail.com> |
| Version | 0.1-0 |
| Date | 2023-01-12 |


## Install `rfnorm` package

- Clone `rfnorm`
- Redirect to `rfnorm/Rpack/rfnorm` directory
- Install from `devtools`

```R
libray(devtools)
install()
```

## Test for 2 dimensions Gibbs Method with the `rfnorm` package

```R
#### Load Packages
## Main package
library(rfnorm)
#### Settings ####
N <- 10000 
Chain_n <- 4
Sig = matrix(
  c(1  , 0.5,
    0.5, 1  ), byrow = T, ncol =2)

#### RNG ####
set.seed(1010)

t_start <- Sys.time()

result = mvrfnorm(N,
c(0,0),
Sig,
init=matrix(c(
    0.5, 1.5,
    0.5, 1.5,
    1.5, 0.5,
    1.5,1.5
    ), byrow = T, ncol = 2), gibbs_chain =  4)

t_end <- Sys.time()
t_end - t_start
#> Time difference of 0.02137399 secs

sample_array <- result[["value"]]
```

We find the running time is now 0.02137399 secs comparing to 0.3977759 secs.

### Performance Illustration: Ergodic Mean Plot

We could verify the RNG (random number generation) quality of the proposed Gibbs by several methods. The following code shows an example to draw the Ergodic mean plot.

One could find more evaluations in the next section, which includes:

- Trace Plot
- Ergodic Mean Plot
- ACF Plot
- [TBD] More is coming...

```R
library(MASS)
library(base)
library(dplyr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(purrr)
library(tidyr)
library(reshape2)
library(latex2exp)
library(cowplot)
library(ggpubr)
library(ks)

#### Ergodic Mean Plot ####
## Ergodic data
Chain_Para1_erg <-
  sample_array[,,1] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta1") %>%
  as.data.frame()
Chain_Para2_erg <-
  sample_array[,,2] %>%
  as.data.frame() %>%
  mutate(Iteration=seq_len(N)) %>%
  melt(id.vars = "Iteration", measure.vars = paste0("chain:",1:4),variable.name = "Chain",value.name = "Para") %>%
  group_by(Chain) %>%
  mutate(ergodic_mean = cumsum(Para)/seq_len(N),
         Theta = "theta2") %>%
  as.data.frame()
Chain_erg <- 
  Chain_Para1_erg %>% 
  rbind(Chain_Para2_erg) %>%
  mutate(Theta_c = recode(Theta,
                          "theta1" = "theta[1]",
                          "theta2" = "theta[2]")
  )
## Implement ergodic mean plot
Chain_erg %>%
  group_by(Chain) %>%
  # top_n(n = -10000, wt = Iteration) %>%
  ggplot(aes(x=Iteration, y=ergodic_mean, colour=Chain, group=Chain)) +
  facet_grid(rows = vars(Theta_c), labeller = label_parsed) +
  geom_line() +
  theme_minimal(base_size = 14) +
  # geom_hline(aes(yintercept = sqrt(2/pi)), linetype = 5, size = 0.3) +
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.grid = element_blank(),
        legend.key = element_rect(
          color = "white", 
          fill = "white"),
        strip.text.y = element_text(size = 15,angle = 0)) +
  theme(axis.title.x = element_text(margin = margin(r = 20)),
        axis.title.y = element_text(margin = margin(r = 15))) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(35, "pt")) +
  theme(legend.title.align = 0.5) +
  ylab("Ergodic Mean") +
  scale_color_discrete(labels = c(expression(paste(theta^(0)==group("(",list(0.5,0.5),")"))),
                                  expression(paste(theta^(0)==group("(",list(0.5,1.5),")"))),
                                  expression(paste(theta^(0)==group("(",list(1.5,0.5),")"))),
                                  expression(paste(theta^(0)==group("(",list(1.5,1.5),")")))))

```



# R Source Code for 2 and more dimensions

check [Source Code](src.md)

## Two Dimension Case FN(0, ρ = 0.5)

| Type | Plot |
|:----|:----|
| **Trace Plot** |  ![](./fig/media/image2.png) |
| **Ergodic Mean Plot** | ![](./fig/media/image3.png) |
| **ACF Plot** | ![](./fig/media/image4.png)|

## Two Dimension Case FN(0, ρ = 0.95)

| Type | Plot |
|:----|:----|
| **Trace Plot** |  ![](./fig/media/image5.png) |
| **Ergodic Mean Plot** | ![](./fig/media/image6.png) |
| **ACF Plot** | ![](./fig/media/image7.png)|


## Five Dimension Case FN(0,I)

| Type | Plot |
|:----|:----|
| **Trace Plot** |  ![](./fig/media/image8.png) |
| **Ergodic Mean Plot** | ![](./fig/media/image9.png) |
| **ACF Plot** | ![](./fig/media/image10.png)|

# Imputation Algorithm

Check [here](./algorithm2/algorithm2.md) for more details.

# TODO

- [ ] R package: Under development, will be released on CRAN.
- [ ] Dim > 1, re-write in `C++`
  - [x] Dim == 2 in `C++`
  - [ ] Dim >= 3 in `C++`
  - [ ] Imputations
  - [ ] Unitest
  - [ ] Vignettes
- [ ] More documentations
