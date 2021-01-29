# MAnorm2 algorithm 


# Estimate Size factor

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
# Estimate interval size
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
convert <- function(y){ log(y / interval.size + offset, base = 2) }
```








# Reference












Tu, S., et al., MAnorm2 for quantitatively comparing groups of ChIP-seq samples. bioRxiv, 2020: p. 2020.01.07.896894. https://doi.org/10.1101/2020.01.07.896894.




