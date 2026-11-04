### Table II. Titration verification (VERIFIED = True)

| check | published | this implementation | pass |
|---|---|---|---|
| V1 logistic r = 3.7 | ~0.75 | 0.637 | True |
| V2 logistic r = 3.575 | ~0.09 | 0.087 | True |
| V3 logistic r = 3.565 (periodic) | ~0 | 0.000 | True |
| V4 bifurcation scan | NL > 0 iff LE > 0 | NL>0 78/78; NL=0 18/18; Spearman 0.85 | True |
| V5 controls (20 each) | NL = 0 | white 0/20 positive; ar2 0/20 positive; limit_cycle 0/20 positive; torus 0/20 positive; logistic_p3 0/20 positive | True |
| V6 flows + Hénon | chaotic > 0, periodic 0 | Lorenz28 0.61, Lorenz160 0.00, Hénon 0.77, MG chaotic 0.33, 0.68, 0.66, MG periodic 0.00, 0.00, 0.00 | True |
