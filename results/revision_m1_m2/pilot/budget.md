# Pilot runtime projection

```
pilot tasks 280, errors 0
ablation (13 configurations): 0.94 h wall
  R0        0.05 h
  S1        0.05 h
  S2        0.05 h
  S3        0.06 h
  S3b       0.06 h
  S4        0.15 h
  S5a       0.08 h
  S5b       0.09 h
  S6        0.05 h
  S7        0.05 h
  R0-BP96F  0.05 h
  R0-BP96   0.06 h
  R0-TEXT   0.15 h
abl2 (R0-ORIG + up to 3 leave-one-out, non-surrogate worst case): 0.34 h
P1 selection + series generation: 0.14 h
headline (R0, R0-BP96, R0-ORIG at 6.8 x R0): 5.56 h
M2 (10200 windows): 0.14 h
TOTAL projected: 7.11 h wall on 4 workers (budget 12 h)
```
