# Folded Normal RNG

# Figures

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


## R packages
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


## Indicator function of random vector X
```{r}
Indicator <- function(x){
  return(ifelse(FALSE %in% ifelse(x >= 0, TRUE, FALSE), 0, 1))
}
```

## Folded Normal Sampling Funcion
```{r}
FN_Sampling <- function(number, mu, sig){
  sample <- abs(rnorm(number, mean = mu, sd = sqrt(sig)))
  density <- function(x, u, s){
    1/sqrt(2 * pi)/sqrt(s) * (exp(-(x - u)^2/(2 * s)) + exp(-(x + u)^2/(2 * s)) * Indicator(x))
  }
  return(list(sample,density(sample, mu, sig)))
}
```
