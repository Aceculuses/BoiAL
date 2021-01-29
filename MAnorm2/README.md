# MAnorm2 algorithm 


# Estimate Size factor

The example CHIP-seq count table

| Sample 1   | Sample 2    | 
| ---------- | ----------- | 
| 16         |      9      |
| 498        |     477     | 
| 54         |      39     |

```
ref <- apply(counts, 1, function(x){ exp(mean(log(x))) })
```
| Sample 1   | Sample 2    |                y           |    ref      |      
| ---------- | ----------- | ---------------------------|-----------  |
| log(16)    | log(9)      |  -> y=log(16)+log(9) / 2   |  ->  exp(y) |
| log(498)   | log(477)    |  -> y=log(498)+log(477) / 2|  ->  exp(y) |
| log(54))   | log(39)     |  -> y=log(54)+log(39) / 2  |  ->  exp(y) |

```
sizeFactor <- apply(counts, 2, function(x){ median(x / ref, na.rm = TRUE) })
```
| Sample 1    |       y        |    Size factor   |
| ----------  | -----------    | ---------------- |
| 16          |   y1=16 / ref  |                  |
| 498         |  y2=477 / ref  | median(y1,y2,y3) |
| 54          |    y3=54/ref   |                  |

| Sample 2    |    y           |   Size factor    |
| ----------- | ---------------|----------------  |
|      9      |  y1=9 / ref    |                  |
|     477     |  y2=477 / ref  |median(y1,y2,y3)  |
|     39      |  y3=39/ref     |                  |

# Reference

Tu, S., et al., MAnorm2 for quantitatively comparing groups of ChIP-seq samples. bioRxiv, 2020: p. 2020.01.07.896894. https://doi.org/10.1101/2020.01.07.896894.




