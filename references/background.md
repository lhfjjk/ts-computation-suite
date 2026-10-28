# Background: why the protocol looks like this

This is the conceptual layer behind the five-stage protocol. It is optional
reading for running a calculation, but it explains why the gates use ΔE, why a
constrained point is not a TS, why selectivity is a flux ratio, and what AI can
and cannot contribute.

## 1. Thermodynamics and kinetics are not two independent criteria

For an elementary step, detailed balance gives

```text
k+ / k- = K = exp(-ΔG°/RT)
```

and the net rate can be written in Marcelin–de Donder form:

```text
v = v+ (1 - exp(-A/RT)),    A = -ΔG = RT ln(K/Q)
```

`v+` is kinetic; the bracket is thermodynamic. The sign of `A` sets the
direction, `v+` sets the speed. Near equilibrium this becomes `v ≈ (v+/RT) A`,
the linear irreversible-thermodynamics regime; the entropy production is
`σ = vA/T ≥ 0`.

The barrier itself depends on the driving force (Marcus/Bell–Evans–Polanyi),
so a thermodynamic quantity is already inside the kinetic barrier. The two
"criteria" are projections of one free-energy surface.

## 2. The unified object: F(ξ) plus mobility

In the overdamped limit the probability density obeys

```text
∂t p = ∂ξ [ D exp(-F/kT) ∂ξ ( exp(F/kT) p ) ]
```

`F` is thermodynamics, `D` is kinetics (friction/mobility), and `p(ξ,t)` answers
"whether, how far, and how distributed". Kramers rates, equilibrium
populations, and relaxation times are limits of this equation. For discrete
states, the same statement is a master equation with local detailed balance.

This is why the protocol does not try to replace ΔG or ΔG‡ with one scalar:
the complete object is a landscape plus transport, or a network with rate
constants. A single number cannot encode direction, speed, and spread.

## 3. What the three questions mean in practice

**Can it proceed?** Not a Boolean. Use the Damköhler number `Da = k τ_obs`:
`Da ≫ 1` proceeds, `Da ≪ 1` is frozen, `Da ~ 1` is partial. A negative overall
ΔG only says a downhill direction exists; it does not guarantee a low-barrier
path. A catalyst changes the path and barrier, not the overall ΔG.

**How far?** In a closed system, the equilibrium extent follows `K` (for
A ⇌ B, `x_eq = K/(1+K)`). In an open/flow system there is no equilibrium; the
steady-state conversion depends on residence time and rate. For a catalytic
cycle, the practical unified criterion is the **energy span**:

```text
TOF = (k_BT/h) exp(-δE/RT),    δE = max(TDI - TDI')
```

The most stable intermediate (thermodynamics) lowers `TDI'` and raises the
span; the highest transition state (kinetics) sets `TDI`. This is why a
"thermodynamic sink" can slow a cycle. The degree of rate control
`X_RC,i = ∂ln(TOF)/∂ln k_i` tells you which step actually controls the rate.

**How distributed?** At equilibrium, molecule counts are binomial (Gaussian for
large N), with variance set by `k_BT/G''`; for a macroscopic sample the spread
is ~1/√N. In an open steady state there is no universal Boltzmann form. For
weakly reversible deficiency-zero networks, the stationary distribution is
product-form Poisson and unique; higher deficiency allows bistability. Off
equilibrium, the analogue of free energy is the quasipotential, which depends
on the dynamics and is not a state function.

**Selectivity** is not a Boltzmann product ratio unless the products
equilibrate. It is a flux branching ratio, `J1/J2`. When a common intermediate
is rapidly equilibrated, this reduces to the Curtin–Hammett ratio set by
ΔΔG‡ (the span difference), not by the product ΔG difference. The
"kinetic vs thermodynamic control" dichotomy is a comparison of the observation
time with the equilibration time.

## 4. Why the protocol uses ΔE, constrained points, and IRC

- **ΔE for minima vs saddle:** a saddle's imaginary frequency is excluded from
  the partition function, so its thermal correction is not comparable to a
  minimum's. A sign reversal between ΔE and ΔG is possible.
- **Constrained point ≠ TS:** a constrained optimization is a minimum on the
  constraint surface, not a stationary point of the full surface. It is a
  reference point and a Hessian source, not a barrier.
- **IRC:** a `SADDLE-1` with the right atom pair can still be a different
  channel. Only the two endpoints establish connectivity.
- **Conformer averaging:** a single conformer's R/S energy difference can be
  conformer-dominated; the observed selectivity is a population-weighted flux,
  not one geometry's number.

## 5. What AI can and cannot do here

AI is useful for **candidate generation** (UniTS-Gen, GOAT), **prefiltering**,
**surrogate barriers** (ML potentials/Δ-ML), and **orchestration** (inputs,
parsing, monitoring). It is not a substitute for:

- defining the elementary step and the spin surface;
- checking connectivity and mode identity;
- validating IRC endpoints;
- quantitative selectivity in the 1–3 kcal window for transition metals;
- extrapolating to new catalysts/ligands outside the training distribution.

The practical division is: **AI proposes, physics filters, DFT validates, human
decides.** The evidence and limits are in [ai-role.md](ai-role.md).

## 6. A note on "new standards"

A universal scalar "reaction standard" is unlikely to be both novel and
correct: thermodynamics has no rate, kinetics alone can violate the second law,
and a scalar cannot encode direction, speed, and spread. A useful standard is a
vector or a function:

```text
(A, J, σ, Var)   or   p(ξ,t)
```

or, for a catalytic cycle, `(δE_span, TOF, branching ratios, ⟨S²⟩, uncertainty)`.
That is the language this suite tries to make operational.
