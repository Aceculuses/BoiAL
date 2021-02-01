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
|    inv.strMatrix|  |    vapply(inv.strMatrix, sum, numeric(1) |
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
Signal1 * weight1 + Signal2 * weight2
```

| sample 1 |  sample 2 |   coef_sample1  |  coef_sample2   |     IntervalMeans       |
|----------|-----------|-----------------|-----------------|-------------------------|
|   4.05   |    3.24   |      0.5        |     0.5         |  4.05 * 0.5 + 3.24 * 0.5|
|   8.93   |    8.89   |      0.5        |     0.5         |  8.93 * 0.5 + 8.89 * 0.5|
|   5.76   |    5.30   |      0.5        |     0.5         |  5.76 * 0.5 + 5.30 * 0.5|





# Reference
Tu, S., et al., MAnorm2 for quantitatively comparing groups of ChIP-seq samples. bioRxiv, 2020: p. 2020.01.07.896894. https://doi.org/10.1101/2020.01.07.896894.




