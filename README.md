# Truthful Aggregation of Budget Proposals with Proportionality Guarantees

## A set of Quadratic Programs with Quadratic Constraints to computer worst-case approximation guarantees

We propose a truthful mechanism for the aggregation of budget proposal for $3$ voters. Follow the [companion paper]() for more details.

The code in file loss_computation_3.py computes worst-case upper bounds for the $ell_1$ distance between an aggregated division decided by our mechanism and
the proportional division of the budget, a metric we call the $\ell_1$-loss.

## Usage

Use the following command to compute the maximum loss.

python loss_calculation_3.py --help --optimal --zero --nonzero

-h, --help: shows this help message
-o, --optimal: only programs with optimal solutions printed
-n, --nonzero: computes the upper bound only for the case x_1>0,x_2>0,x_3>0
-z, --zero: computes the upper bound only for the case x_3=0

## Dependancies

The calculator uses the [Gurobi](www.gurobi.com) mathematical optimization solver. If your eligible, you can acquire an academic license, following the instructions [here](https://www.gurobi.com/academia/academic-program-and-licenses/).
