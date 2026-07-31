# Directional Reachable-Set Optimization under Matrix Uncertainty

This repository contains MATLAB illustrations of a directional reachable-set optimization method for linear systems with uncertainty in the state and input matrices. The method selects deterministic matrices \(A^*\) and \(B^*\) from prescribed uncertainty sets so that the reachable set is shaped or extended along a chosen direction of interest.

The repository has three main components:

1. **Illustrations** — applies the optimization method to a two-state example and to the ADMIRE fighter-jet model.
2. **Comparison with CORA** — compares the optimized deterministic propagation with CORA's full uncertain linear reachability computation.
3. **Road scenario** — demonstrates how tighter directional predictions can reduce conservatism in an overtaking scenario with uncertain vehicle dynamics.

> This is research and demonstration code. The examples are intended to illustrate the method and are not a certified safety or control implementation.

## Method overview

Consider the uncertain linear system

\[
\dot{x}(t) = A x(t) + B u(t),
\qquad
A \in \mathcal{A},\quad B \in \mathcal{B},\quad u(t) \in \mathcal{U}.
\]

For a direction of interest \(d\) and time horizon \(T\), the code searches for matrices \(A^*\in\mathcal{A}\) and \(B^*\in\mathcal{B}\) that increase the reachable extent in that direction. The implementation uses the terminal costate relation

\[
p_0(A) = e^{A^\top T}d
\]

and a two-stage procedure:

1. Optimize the state matrix and initial-state vertex to obtain \(A^*\) and \(x_0^*\).
2. With \(A^*\) fixed, optimize the input matrix and control vertex to obtain \(B^*\) and \(u_0^*\).

The selected deterministic system

\[
\dot{x}(t)=A^*x(t)+B^*u(t)
\]

can then be propagated using standard deterministic reachability routines. In the included examples, this produces a directionally targeted reachable set while avoiding the full cost and conservatism of propagating every admissible matrix realization.

## Repository structure

```text
.
├── Illustrations/
│   ├── FullA_BOptimization.m
│   └── AdmireFighterJet_Better.m
│
├── Comparison_with_CORA_uncertainLinear/
│   ├── FullA_BOptimization_CORA.m
│   └── FullA_BOptimization_CORA_ADMIREJET.m
│
└── RoadScenario/
    ├── SPOT_Vs_Me9.m
    ├── straight_window1_green.mp4
    ├── straight_window2_red.mp4
    ├── curved_window1_green.mp4
    └── curved_window2_red.mp4
```

## Requirements

- MATLAB
- [CORA](https://cora.in.tum.de/) for the scripts in:
  - `Comparison_with_CORA_uncertainLinear/`
  - `RoadScenario/`
- CORA must be installed and added to the MATLAB path.

The scripts in `Illustrations/` use standard MATLAB functionality and do not require CORA.

A typical CORA setup is:

```matlab
addpath(genpath('/path/to/CORA'));
savepath;
```

## Running the examples

Run MATLAB from the repository root, or add the repository folders to the MATLAB path.

### 1. Basic two-state illustration

```matlab
run('Illustrations/FullA_BOptimization.m')
```

This example:

- represents both \(A\) and \(B\) as matrix zonotopes;
- optimizes \(A^*\), \(x_0^*\), \(B^*\), and \(u_0^*\);
- constructs reachable-set boundary samples using costate-guided, Hamiltonian-maximizing controls; and
- compares the nominal reachable sets with the optimized reachable sets over time.

The nominal sets are shown in red, the optimized sets in green, and the initial set in blue.

### 2. ADMIRE fighter-jet illustration

```matlab
run('Illustrations/AdmireFighterJet_Better.m')
```

This example applies input-matrix optimization to a linearized ADMIRE fighter-jet model with:

- states \(x=[p,q,r]^\top\), representing roll, pitch, and yaw rates;
- four control surfaces; and
- a selected direction of interest in the state space.

The script computes an optimized input matrix within a Frobenius-norm uncertainty bound and compares the nominal and optimized three-dimensional reachable sets.

## Comparison with CORA

### Two-state comparison

```matlab
run('Comparison_with_CORA_uncertainLinear/FullA_BOptimization_CORA.m')
```

### ADMIRE comparison

```matlab
run('Comparison_with_CORA_uncertainLinear/FullA_BOptimization_CORA_ADMIREJET.m')
```

These scripts compare three systems:

- **Uncertain system:** CORA propagation using `linParamSys` with matrix-zonotope uncertainty in \(A\) and \(B\).
- **Nominal system:** deterministic propagation using the center matrices \(A_c\) and \(B_c\).
- **Optimized system:** deterministic propagation using the selected matrices \(A^*\) and \(B^*\).

The scripts report:

- reachability computation time;
- reachable-set projections;
- state-component bounds over time;
- upper and lower support values along \(d\); and
- directional width at the final time.

The optimized propagation is designed to retain the behavior relevant to the selected direction while replacing the complete uncertain propagation with one deterministic system. This can produce a tighter directional prediction and a lower reachability-computation cost. Exact timing improvements depend on the model, uncertainty size, time horizon, and CORA settings.

## Road scenario

The road example maps reachable tubes from a Frenet frame into straight and curved road geometries. The uncertain vehicle model includes longitudinal and lateral dynamics. Two prediction strategies are visualized:

- **Green:** the directionally optimized deterministic reachable tube;
- **Red:** the full uncertain reachable tube computed with CORA.

The tighter green prediction leaves more usable free space and illustrates how the ego vehicle can follow a less conservative overtaking trajectory. The larger red prediction blocks the corridor and causes the illustrated ego trajectory to wait.

### MATLAB filename note

`RoadScenario/SPOT_Vs_Me9.m` defines the primary function

```matlab
overtake_frenet_animation(roadType)
```

MATLAB normally requires the primary function name and filename to match. Rename the file before running it:

```text
SPOT_Vs_Me9.m  ->  overtake_frenet_animation.m
```

Then run:

```matlab
cd RoadScenario

overtake_frenet_animation('straight')
overtake_frenet_animation('curved')
```

Set the following flag near the top of the function to control video generation:

```matlab
MAKE_VIDEO = true;
```

The generated demonstrations are also included in the repository:

- [Straight road — optimized green prediction](RoadScenario/straight_window1_green.mp4)
- [Straight road — standard red prediction](RoadScenario/straight_window2_red.mp4)
- [Curved road — optimized green prediction](RoadScenario/curved_window1_green.mp4)
- [Curved road — standard red prediction](RoadScenario/curved_window2_red.mp4)

> The ego trajectories in this example are scripted for visualization. The reachable sets illustrate the available planning space but do not directly gate the vehicle motion or implement a closed-loop safety filter.

## Customizing an example

The main parameters are defined near the beginning of each MATLAB file. Common quantities to modify include:

```matlab
Ac       % center of the state-matrix uncertainty set
GA       % state-matrix zonotope generators
Bc       % center of the input-matrix uncertainty set
GB       % input-matrix zonotope generators
V_X0     % vertices of the initial set
Ubound   % control bounds
d        % direction of interest
tf       % reachability horizon
```

For the CORA comparisons, the principal reachability settings are:

```matlab
options.timeStep
options.taylorTerms
options.zonotopeOrder
options.intermediateTerms
```

Runtime and conservatism are sensitive to these parameters, so comparisons should use identical initial sets, input sets, horizons, and compatible numerical settings.

## Interpretation of the results

The optimized set is not intended to reproduce the full uncertain reachable set in every direction. It is a directionally selected reachable set associated with the optimized matrices \(A^*\) and \(B^*\). Its main purpose is to preserve or improve the reachable behavior relevant to a specified planning or control direction while enabling faster deterministic propagation.

Consequently:

- use the full uncertain CORA set when an enclosure of all admissible uncertain trajectories is required;
- use the optimized set when the application is concerned with a particular direction or maneuver and the assumptions of the method are appropriate; and
- validate the resulting set against the uncertainty model before using it in a safety-critical workflow.

## Acknowledgment

The uncertain reachability comparisons and road-scenario visualizations use the CORA toolbox.

## Citation

A publication citation can be added here when the associated paper is available:

```bibtex
@article{das_reachable_set_optimization,
  author  = {Das, Hrishav and Ornik, Melkior},
  title   = {Directional Reachable-Set Optimization under Matrix Uncertainty},
  journal = {To appear},
  year    = {2026}
}
```
