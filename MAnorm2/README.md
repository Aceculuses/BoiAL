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
| Sample 1   | Sample 2    |       
| ---------- | ----------- | 
| log(16)    | log(9)      |  -> y=log(16)+log(9) / 2     ->  exp(y)
| log(498)   | log(477)    |  -> y=log(498)+log(477) / 2  ->  exp(y)
| log(54))   | log(39)     |  -> y=log(54)+log(39) / 2    ->  exp(y)


# Reference




