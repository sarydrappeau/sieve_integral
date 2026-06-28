# Computing of sieve integrals using LattE

This repository contains Python 3.11.2 / SageMath 10.7 script files providing methods to compute integrals of the shape

$$ I(P) = \int_P \frac{dt_1 \dotsb dt_k}{t_1 \dotsb t_k} $$

where $P \subset {\mathbb R}_{>0}^k$ is a convex polytope. Also implemented is the variant

$$ I^*(P) = \int_{P \cap H_d} \frac{dt_1 \dotsb dt_{k-1}}{t_1 \dotsb t_k} $

where $H = \{(t_1, \dotsc, t_k) \in R^k, \sum_i t_i = 1\}$.


