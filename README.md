# ExplicitFractionalLaplacianRieszExamples.jl

Companion numerical validation for
[*Explicit fractional Laplacians and Riesz potentials of classical functions*](https://arxiv.org/abs/2311.10896)
by Timon S. Gutleb and Ioannis P. A. Papadopoulos.

This is research companion code, not a registered package. The
displayed Meijer-G forms are evaluated with [MeijerG.jl](https://doi.org/10.5281/zenodo.19430427).
Core formulas are compared with defining-operator quadrature, while derived
representations are compared with operator-validated formulas.

Tests use deterministic seeded random points within the paper's sufficient
parameter ranges: normally 100 per formula, or 20 for the oscillatory formulas
in Table 2. Finite-tail quadratures must agree when the cutoff is doubled.

The numerical solver used for Figure 2 is provided separately in
[FractionalFrames.jl](https://github.com/ioannisPApapadopoulos/FractionalFrames.jl).
