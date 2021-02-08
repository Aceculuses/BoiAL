# TF dynamics 

# Hypergeometric distribution

![image](https://github.com/Aceculuses/BoiAL/blob/main/TFdynamics/venn.hyper.png)

B = 10 genes
C = 8 genes
overlap = 2

To calculate p-value of overlap genes

-C = 8
-B = 6

```
phyper(overlap,C,-C,B) = phyper(2,8,8,10) = 0.003496503
phyper(overlap,B,-B,C) = phyper(2,10,6,8) = 0.003496503
```









