# StructuralTransformations.jl

StructuralTransformations.jl is a package for transforming ModelingToolkit models into structurally better
versions for improved numerical simulation. These transformations are compiler passes on the model form,
such as index reduction, which output new models that are more numerically stable and faster to solve.

## License

StructuralTransformations.jl is free for academic and research use but requires a license for commercial
use. Please see the full License file for more details.

## Example

```julia
using StructuralTransformations
using ModelingToolkit
using OrdinaryDiffEq
using Plots

# Define some variables
@parameters t L g
@variables x(t) y(t) T(t)
@derivatives D'~t

eqs = [D(D(x)) ~ T*x,
       D(D(y)) ~ T*y - g,
       0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, T], [L, g], name=:pendulum)

# Turn into a first order differential equation system
first_order_sys = ModelingToolkit.ode_order_lowering(pendulum2)

# Perform index reduction to get an Index 1 DAE
new_sys = StructuralTransformations.dae_index_lowering(first_order_sys)

u0 = [
  D(x)    => 0.0,
  D(y)    => 0.0,
  x       => 1.0,
  y       => 0.0,
  T       => 0.0
]

p = [
    L => 1.0,
    g => 9.8
]

prob_auto = ODEProblem(new_sys,u0,(0.0,100.0),p)
sol = solve(prob_auto, Rodas5());
plot(sol, vars=(D(x), y))
```

## Methods

- `dae_index_lowering(sys::ODESystem)`: Performs the Pantelides Algorithm to
  transform a higher index DAE to an Index 1 DAE.
