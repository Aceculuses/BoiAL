# MAnorm2 Algorithm 


# Within group normalization

# Estimate Size Factor

The example CHIP-seq count table

| Sample 1   | Sample 2    | 
| ---------- | ----------- | 
| 16         |      9      |
| 498        |     477     | 
| 54         |      39     |

```
ref <- apply(counts, 1, function(x){ exp(mean(log(x))) }) #base = e
```
| Sample 1   | Sample 2    |                y            |    ref            |      
|------------|-------------|-----------------------------|-------------------|
| log(16)    | log(9)      |  y1=log(16)+log(9) / 2      |  ref1=exp(y1)     |
| log(498)   | log(477)    |  y2=log(498)+log(477) / 2   |  ref2=exp(y2)     |
| log(54))   | log(39)     |  y3=log(54)+log(39) / 2     |  ref3=exp(y3)     |

```
sizeFactor <- apply(counts, 2, function(x){ median(x / ref, na.rm = TRUE) })
```
| Sample 1    |       y        |    Size factor   |
|-------------|----------------|------------------|
| 16          |  y1=16  / ref1 |                  |
| 498         |  y2=477 / ref2 | median(y1,y2,y3) |
| 54          |  y3=54  / ref3 |                  |

| Sample 2    |    y           |   Size factor    |
|-------------|----------------|------------------|
|      9      |  y1=9 / ref1   |                  |
|     477     |  y2=477 / ref2 |median(y1,y2,y3)  |
|     39      |  y3=39 / ref3  |                  |

Determine which size factor could be used, and the baseline is the chosen size factor, baseline is the index of the chosen size factor, which is refer to baseline samples. eg. baseline <- which.min(abs(log(size.factor))) =1 that means sample1 is the baseline sample.
```
if (all(is.na(size.factor))) {
  stop("Failed to estimate the size factors of samples.
You may specify the baseline sample explicitly")
}
if (length(size.factor) == 2) {
  # To avoid numeric uncertainty when the two size factors
  # are reciprocal to each other
  baseline <- which.min(size.factor)
} else {
  baseline <- which.min(abs(log(size.factor)))
}
base.flag <- TRUE
```
# Estimate Interval Size
There are three choices to set interval size

1.Interval size can be self-define number such 300bp

2.Interval size is set as TRUE. This will calculate each peak size

3.Interval size is set as FALSE. This will set all peak size as 1 for every peak

```
offset <- 0.5
offset <- as.numeric(offset)[1]

#'For condition one
interval.size <- as.numeric(interval.size) / 1000

#'For condition two
if (interval.size[1]) {
    if (!("start" %in% names(x))) stop("Missing \"start\" variable in x")
    if (!("end" %in% names(x))) stop("Missing \"end\" variable in x")
    interval.size <- (x$end - x$start) / 1000
 
 #'For condition three
 interval.size <- 1

#'All the code
if (is.logical(interval.size)) {
  if (interval.size[1]) {
    if (!("start" %in% names(x))) stop("Missing \"start\" variable in x")
    if (!("end" %in% names(x))) stop("Missing \"end\" variable in x")
    interval.size <- (x$end - x$start) / 1000
  } else {
    interval.size <- 1
  }
} else {
  interval.size <- as.numeric(interval.size) / 1000
}
```

# Reads Counts to Signal Intensities

```
Set interval size = 1, offset = 0.5
convert <- function(y){ log(y / interval.size + offset, base = 2) }
```
|   Sample 1  |           y         |
|-------------|---------------------|
| 16          |   log2((16/1)+0.5)  |
| 498         |   log2((498/1)+0.5) |
| 54          |   log2((54/1)+0.5)  |

# MA Correlation Matrix

Select Occupancy of baseline sample 

| Peaks   |     occupancy   |
|---------|-----------------|
|  1      |       0         |
|  2      |       1         |
|  3      |       0         |

0 for no peak, 1 for yes peak

Generate correlation matrix for samples, and only calculate up triangle values
```
j <- length(cnt)
    MA.cor <- matrix(0, nrow = j, ncol = j)
    rownames(MA.cor) <- colnames(MA.cor) <- names(cnt)
    for (i in 1:j) {
        for (k in i:j) {
            if (k == i) next
            flag <- occupancy[[i]] & occupancy[[k]] & common.peak.regions
            MA.cor[i, k] <- MA.pcc(cnt[[i]][flag], cnt[[k]][flag])
        }
    }
```

|           |   Sample 1  |    Sample 2  |
|-----------|-------------|--------------|
|  Sample 1 |             |   Cor value  |
|  Sample 2 |             |              |

MA Pearson Correlation Coefficient ----MA.pcc
```
MA.pcc <- function(x, y) {
    if (length(x) < 2) return(NA_real_)
    a <- x + y
    b <- y - x
    if (sd(a) == 0 || sd(b) == 0) return(0)
    cor(a, b, method = "pearson")
}
```
|  a=x+y  | Sample 1(x)   | Sample 2(y) |    b=y-x   |
|---------| ------------- | ----------- | -----------|
|   25    |      16       |      9      |    -7      |
|   975   |     498       |     477     |    -21     |
|   93    |      54       |      39     |    -15     |
|  sd(a)  |               |             |   sd(b)    |

MA.pcc(a, b)

# MA Normalization Coefficients Matrix
Set flag for common peak regions. 
```
flag <- base.ocupy & occupancy[[i]]
```

| Baseline sample   |   sample    |    flag    |    peaks    |
| ----------------- | ----------- |------------|-------------|
| 0                 |      1      |   FALSE    |   Different |
| 1                 |      1      |   TRUE     |   Common    |
| 0                 |      0      |   FALSE    |   Different |

```
#'The number of common peaks
norm.coef$common.peak.regions[i] <- sum(flag)
```
```
#' Deduce MA Normalization Coefficients
#'
#' @param baseline A numeric vector representing the baseline signal intensity.
#' @param to.norm A numeric vector representing the sample to be normalized.
#' @return \code{c(slope, intercept)}
normCoef <- function(baseline, to.norm) {
    if (length(baseline) < 2) {
        stop("Too few common peak regions to perform the MA normalization", call. = FALSE)
    }
    m1 <- mean(baseline)
    m2 <- mean(to.norm)
    s1 <- sd(baseline)
    s2 <- sd(to.norm)
    if (s1 == 0 && s2 == 0) return(c(1, m1 - m2))
    if (s1 == 0 || s2 == 0) {
        stop("Common peak regions are associated with constant signal intensities in some sample.
Unable to perform the MA normalization", call. = FALSE)
    }
    slope <- s1 / s2
    intercept <- m1 - slope * m2
    c(slope, intercept)
}
```
| Baseline Sample Intensities | to.norm Sample Intensities|      slope   |     intercept                 |
| ----------------------------|---------------------------|--------------|-------------------------------|
| s1=log2((16/1)+0.5)         | x1=log2((9/1)+0.5)        |              |                               |  
| s2=log2((498/1)+0.5)        | x2=log2((477/1)+0.5)      |              |                               |   
| s3=log2((54/1)+0.5)         | x3=log2((39/1)+0.5)       |              |                               |
|mean1                        |    mean2                  |              |   mean1 - (sd1 / sd2) * mean2 |
|sd1                          |    sd2                    |    sd1 / sd2 |                               |            

# Linear Normalization for Signal Intensities. 
```
#' res[1] <- slope, res[2] <- intercept
#'y = ax + b
#'y <- slope * x + intercept
if (i == baseline) next
cnt[[i]] <- cnt[[i]] * res[1] + res[2]
```
| Baseline Sample Intensities | to.norm Sample Intensities|    Slope     |     Intercept                 |      Linear Normalization        |
| ----------------------------|---------------------------|--------------|-------------------------------|----------------------------------|
| s1=log2((16/1)+0.5)         | x1=log2((9/1)+0.5)        |              |                               |      y1=x1 * slope + intercept   |
| s2=log2((498/1)+0.5)        | x2=log2((477/1)+0.5)      |              |                               |      y2=x2 * slope + intercept   |
| s3=log2((54/1)+0.5)         | x3=log2((39/1)+0.5)       |              |                               |      y3=x3 * slope + intercept   |
|mean1                        |    mean2                  |              |   mean1 - (sd1 / sd2) * mean2 |                                  |
|sd1                          |    sd2                    |    sd1 / sd2 |                               |                                  |

# MA plot
M - minus 

A (add) =  (log2(y) + log2(x)) / 2

M (minus, log fold change) = log2(y) - log(x) = log2(y/x) 

x - baseline samole

y - comparison sample

| Baseline Sample Intensities | to.norm Sample Intensities|    Slope     |     Intercept                 |      Linear Normalization      |    M     |     A      |
| ----------------------------|---------------------------|--------------|-------------------------------|--------------------------------|----------|------------| 
| s1=log2((16/1)+0.5)         | x1=log2((9/1)+0.5)        |              |                               |      y1=x1 * slope + intercept |  y1 - s1 | (y1+s1) / 2|
| s2=log2((498/1)+0.5)        | x2=log2((477/1)+0.5)      |              |                               |      y2=x2 * slope + intercept |  y2 - s2 | (y2+s2) / 2|
| s3=log2((54/1)+0.5)         | x3=log2((39/1)+0.5)       |              |                               |      y3=x3 * slope + intercept |  y3 - s3 | (y3+s3) / 2|
|mean1                        |    mean2                  |              |   mean1 - (sd1 / sd2) * mean2 |                                |          |            |
|sd1                          |    sd2                    |    sd1 / sd2 |                               |                                |          |            |


# Between group normalization

# bioCond
Set Occupancy

```
num <- apply(occupancy, 1, function(y){ sum(as.logical(y), na.rm = TRUE) })
```
| Sample1_Occupancy   |     Sample2_Occupancy   |    num   |
|---------------------|-------------------------|----------|
|          0          |           0             |     0    |
|          1          |           1             |     2    |
|          1          |           0             |     1    |

```
occupy.num =  1
length(num) = 3
occupy.num <- as.numeric(occupy.num)
occupy.num <- rep_len(occupy.num, length.out = length(num))
```
|   occupy.num   |      num    |     num >= occupy.num   |    occupancy   |
|----------------|-------------|-------------------------|----------------|
|       1        |      0      |         FALSE           |      FALSE     |
|       1        |      2      |         TRUE            |      TRUE      |
|       1        |      1      |         TRUE            |      TRUE      |


setWeight

```
weight = NULL
strMatrix = NUL
x <- list(name = name, norm.signal = norm.signal, occupancy = occupancy)
setWeight(x, weight = weight, strMatrix = strMatrix)

By default, all the ChIP-seq samples belonging to the bioCond have the same weight for
estimating the mean signal intensities of genomic intervals

weight <- matrix(1, nrow = 1, ncol = m)  # Same weight
```
|      w1    |      w2    |
|------------|------------|
|      1     |      1     |

```
Set inverse Matrix of strMatrix
strMatrix * inv.strMatrix = 1
inv.strMatrix <- lapply(strMatrix, solve)
```
Set scale.var

```
vapply(inv.strMatrix, sum, numeric(1)
```
|    inv.strMatrix| .|    vapply(inv.strMatrix, sum, numeric(1) |
|---------|----------|------------------------------------------|
|     1   |     0    |             2                            |
|     0   |     1    |                                          |

```
scale.var <- 1 / rep_len(vapply(inv.strMatrix, sum, numeric(1)), length.out = n)
```
|  scale.var   |    length.out(n) |
|--------------|------------------|
|    1 / 2     |          1       |
|    1 / 2     |          2       |
|     ...      |         ...      |
|    1 / 2     |          n       |

Set IntervalMeans

Signal1 * weight1 + Signal2 * weight2

Norm.signal

| sample 1 |  sample 2 |
|----------|-----------|
|   4.05   |    3.24   |
|   8.93   |    8.89   |
|   5.76   |    5.30   |

```
coef <- lapply(inv.strMatrix, function(m) {
        y <- colSums(m)
        y / sum(y)
    })
```

|  inv.strMatrix|  .  |  colSums(m) |  sum(y)  |   coef  |
|---------|-----------|-------------|----------|---------|
|     1   |     0     |       1     |     2    |   0.5   |
|     0   |     1     |       1     |          |   0.5   |

```
n <- length(inv.strMatrix) # n = 1 
index <- rep_len(1:n, length.out = nrow(norm.signal))
```
|   index   |
|-----------|
|    1      |
|    1      |
|    1      |

```
x <- norm.signal
vapply(1:nrow(x), function(i) {sum(x[i, ] * coef[[index[i]]])}, numeric(1))
```
coef 
|   sample 1  |   sample 2   |
|-------------|--------------|
|      0.5    |      0.5     |

|   n   |  index |  index[n] |           coef[index[n]]             |    
|-------|--------|-----------|--------------------------------------|
|   1   |    1   |index[1]= 1|   coef[index[1]] = coef[1] = 0.5, 0.5|
|   2   |    1   |index[2]= 1|   coef[index[2]] = coef[1] = 0.5, 0.5|
|   3   |    1   |index[3]= 1|   coef[index[3]] = coef[1] = 0.5, 0.5|

```
weight1 = weight2 = 0.5
IntervalMeans = Signal1 * weight1 + Signal2 * weight2
```

| sample 1 |  sample 2 |   coef_sample1  |  coef_sample2   |     IntervalMeans       |
|----------|-----------|-----------------|-----------------|-------------------------|
|   4.05   |    3.24   |      0.5        |     0.5         |  4.05 * 0.5 + 3.24 * 0.5|
|   8.93   |    8.89   |      0.5        |     0.5         |  8.93 * 0.5 + 8.89 * 0.5|
|   5.76   |    5.30   |      0.5        |     0.5         |  5.76 * 0.5 + 5.30 * 0.5|


Set intervalVars
```
coef <- lapply(inv.strMatrix, function(m) {
        y <- matrix(colSums(m), nrow = 1)
        m - (t(y) %*% y) / sum(y)
    })
```

|     inv.strMatrix| .|
|---------|-----------|
|     1   |     0     |
|     0   |     1     |

|   y     |       |
|---------|-------|
|    1    |   1   |

|  t(y)   |    
|---------|
|    1    |      
|    1    |

|   t(y)*y  |      |    sum(y)  |  
|-----------|------|------------|
|     1     |   1  |     2      |
|     1     |   1  |            |
                                
|     |   [1]  |   [2]  |
|-----|--------|--------|
| [1] |  0.5   |   0.5  |
| [2] |  0.5   |   0.5  |

```
m - (t(y) %*% y) / sum(y)
```
|   |   [1]   |     [2]   |      |    [1]  |  [2]  |      |   [1]   |   [2]   |    
|---|---------|-----------|------|---------|-------|------|---------|---------|
|[1]|     1   |     0     |   -  |    0.5  |   0.5 |   =  |    0.5  |   -0.5  |
|[2]|     0   |     1     |      |    0.5  |   0.5 |      |   -0.5  |   0.5   |

coef 
|      | [1,]  |  [2,]  |
|------|-------|--------|
| [,1] |  0.5  |  -0.5  |
| [,2] | -0.5  |   0.5  |


intervalVars
```
vapply(1:nrow(x), function(i) {
        y <- x[i, , drop = FALSE]
        (y %*% coef[[index[i]]] %*% t(y))[1, 1]
    }, numeric(1)) / (ncol(x) - 1)
```
y; i =1;   s = y %*% coef[[index[i]]]
|    | [,1] |  [,2] |      | [1,]  | [2,] |    |      [1,]                 |             [2,]           |
|----|------|-------|------|-------|------|----|---------------------------|----------------------------|
|[1,]|4.05  | 3.24  |  %*% |  0.5  |  -0.5|  = |  4.05 * 0.5 + 3.24 *(-0.5)|4.05 * (-0.5) + 3.24 * (0.5)|
|[2,]|      |       |      | -0.5  | 0.5  |    |                           |                            | 

z = s %*% t(y)
|    |      [1,]                         |             [2,]                   |   |  [1,]  |   |    [1,]                              |
|----|-----------------------------------|------------------------------------|---|--------|---|--------------------------------------|
|[1,]|  4.05 * 0.5 + 3.24 *(-0.5)= 0.405 |4.05 * (-0.5) + 3.24 * (0.5)= -0.405|%*%|  4.05  | = | 0.405 * 4.05 +(-0.405)*3.24 = 0.32805|
|[2,]|                                   |                                    |   |  3.24  |   |                                      |

df = (ncol(x) - 1)

intervalVars = z / df 

# normBioCond

Signal sample.mean are calculated by biocond which stands for interValmeans. Occupancy is calculated before. 

interValmeans

|    s1     |     s2    |
|-----------|-----------|
|   3.64    |    4.37   |
|   8.91    |    8.84   |  
|   5.53    |    5.28   |

Occupancy

|  o1   |    o2  |
|-------|--------|
|FALSE  |   TRUE |
|TRUE   |   TRUE |
|FALSE  |   TRUE |

```
subset <- apply(occupancy, 1, all)
```
| subset  |
|---------|
|  FALSE  |
|  TRUE   |
|  FALSE  |

```
#Select the overlapping signal, the signals stand for common peaks
#Normalize common peaks between two groups
temp <- signal[subset, , drop = FALSE]
ref <- rowMeans(temp)
```
temp (common peaks)

|  s1  |   s2  |   ref  |
|------|-------|--------|
| 8.91 | 8.84  |  8.87  |
| 8.39 | 8.40  |  8.39  |
| 6.77 | 6.56  |  6.66  |

```
log2.size <- apply(temp, 2, function(x){ median(x - ref) })
size.factor <- 2 ^ log2.size
```

|  s1  |  ref  |    logs.size  |      median    |
|------|-------|---------------|----------------|
| 8.91 | 8.87  |  y1=8.87-8.91 |                |
| 8.39 | 8.39  |  y2=8.39-8.39 |median(y1,y2,y3)|
| 6.77 | 6.66  |  y3=6.66-6.77 |                |

Repeat for sample2

size factor is converted by 2^

baseline is chosen by min(log2.size), which is same as within group normalization process.

```
#implement same normalization algorithm 
#convert is assigned identity which means it is no need to transfer reads count into signal intensity
#common.peak.region is assigned common.peak.region which means we don't need to find common peaks
norm <- normalize(temp, 1:n, (n + 1):(n * 2), baseline = baseline,convert = identity, common.peak.regions = common.peak.regions)
```

Between group normalization is same as within group process. One difference is that it normalize the common peaks between two groups. 
The rest results are remain same. we cannot normalize unique peaks, just only for common peaks between groups.

# fitMeanVarCurve ---- Regression process

# estimateVarRatio

initialize ratio.var

|  condition 1 |   condition 1 |     |   condition2 |   condition2   |
|--------------|---------------|-----|--------------|----------------|
|   rep1       |      rep2     |     |    rep1      |      rep2      |
|   3.38       |      2.52     |     |    4.55      |      4.20      | 
|   8.613      |      8.57     |     |    8.77      |      8.90      |
|   5.21       |      4.72     |     |    5.47      |      5.08      |


```
# n = 2
# m = 2,2
#noRep = FALSE, FALSE
n <- length(conds)
m <- vapply(conds, function(cond){ ncol(cond$norm.signal) }, numeric(1))
noRep <- m <= 1

#var.level = NA, NA
var.level <- rep(NA_real_, n)

#sample var is calculated in biocond step, which is intervalVars
obs.vars <- lapply(conds[!noRep], function(cond){ cond$sample.var })
```
obs.vars

|  s1.Var |   s2.Var  |
|---------|-----------|
|  0.369  |   0.0618  |
| 0.00074 |   0.00776 |
|   0.12  |   0.074   |

```
temp <- lapply(conds[!noRep], function(cond){ cond$occupancy })
subset <- apply(as.data.frame(temp), 1, all)
```

tmp, subset

|   s1.o      |     s2.o    |    subset   |
|-------------|-------------|-------------|
|    FALSE    |    TRUE     |     FALSE   | 
|    TRUE     |    TRUE     |     TRUE    |
|    FALSE    |    TRUE     |     FALSE   |

```
var.level[!noRep] <- estimateSizeFactors(as.data.frame(obs.vars), subset)
```
Each sample size factor is calcualted by estimateSizeFactors() used before.

```
#This step is same as baseline selection.
#base.cond is the index of minimun var.level
base.cond <- which.min(abs(log(var.level)))

#Choose baseline condition
base.index <- base.cond
base.cond <- conds[[base.cond]]
```

```
varRatio <- function(cond1, cond2, invariant = NULL) {
    if (is.null(invariant)) {
        f <- cond1$occupancy & cond2$occupancy
    } else {
        invariant <- as.numeric(invariant)[1]
        if (invariant < 0) stop("invariant must be non-negative")
        f <- abs(cond2$sample.mean - cond1$sample.mean) <= invariant
    }
    if (!any(f)) return(NA_real_)

    med <- median(cond2$sample.var[f] / cond1$sample.var[f], na.rm = TRUE)
    if (is.na(med) || is.infinite(med) || med == 0) return(NA_real_)
    med / qf(0.5, ncol(cond2$norm.signal) - 1, ncol(cond1$norm.signal) - 1)
}
```

|      base.sample.var     |   to.compare.sample.var    |     ratio      |    median        |
|--------------------------|----------------------------|----------------|------------------|
|   i1 = 7.767096e-03      |    j1 = 2.379149e-02       |  r1 = i1 / j1  |                  |
|   i2 =  5.021998e-03     |    j2 = 7.101262e-03       |  r2 = i2 / j2  | median(r1,r2,r3) |
|   i3 =  2.379149e-02     |    j3 = 3.488587e-02       |  r3 = i3 / j3  |                  |

```
After median of variation ratio is calculated, We want to compare the results with F distribution,
when the quantile is 50% to gain the variation ratio. 
As we choose median of the ratio, therefore we want to compare with the middle point of F distribution, therefore, 0.5 is chosen.
Two degree of freedom were determined by the number of sample replicates in each condition. 
df1 = 2 replicates - 1, df2 = 2 replicates - 1
The ratio between those two resutls is defined as variantion ratio between two conditions. 
Therefore, in order to have variation ratio, we need to have at least
two replicates in each condition.
```
```
Why use F distribution ?

F distribution is used in three area: 1. sampling from two normal distribution. 2. Analysis of variance (ANOVA). 3. Regression analysis. 

In this case, F distribution is used to measure the the variance between two conditions.

What is quantile ?
quantile is the cut off value that certian percentage of density will be chosen.
eg. we need to find a value that will split the F prbolility density into two equal parts in given two degree of freedom
qf(0.5,df1,df2) = quantile values or cut-off value.
In this case, the quantile value can been seen as median number of F distribution.
```

# meanVarParaFit

```
fitMeanVarCurve <- function(conds, ratio.var = estimateVarRatio(conds),
                            method = c("parametric fit", "local regression"),
                            occupy.only = TRUE, range.residual = c(1e-4, 15),
                            max.iter = 50, init.coef = NULL, args.lp = list(nn = 0.7),
                            args.locfit = list(), verbose = TRUE)
                            
occupy.only = TRUE. This parameter can be set as FALSE. If TURE, only common peaks is going to be fitted
```

|         |   means  |  1 / 2 ^means  |
|---------|----------|----------------|
| c1.mean |   2.95   |  1 / 2 ^ 2.95  |
|         |   8.59   |  1 / 2 ^ 8.59  |
|         |   4.97   |  1 / 2 ^ 4.97  |
| c2.mean |   4.37   |  1 / 2 ^ 4.37  |
|         |   8.84   |  1 / 2 ^ 8.84  |
|         |   5.28   |  1 / 2 ^ 5.28  |


```
x <- means
y <- vars

fit <- glm(y ~ x, family = Gamma(link = "identity"), subset = good, weights = weight, start = coef)
```

```
Generalize Linear Model (GLM)

This process fit means to vars, glm(vars ~ means), we assume that the data follow exponent distribution ( gamma is one of the exponent distribution family),
One key point in GLM is to deal with errors which is not in normal distribution. 
```

The better explanation can be found here [GLM](https://zhuanlan.zhihu.com/p/110268967)





# Reference
Tu, S., et al., MAnorm2 for quantitatively comparing groups of ChIP-seq samples. bioRxiv, 2020: p. 2020.01.07.896894. https://doi.org/10.1101/2020.01.07.896894.




