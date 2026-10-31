# Model description

The model computes the excitatory drive that a ganglion cell receives
from an array of diffuse bipolar cells. It has five stages:

1. Place bipolar subunits on a jittered hexagonal mosaic.
2. Filter the stimulus with each subunit's spatial and temporal receptive field.
3. Add noise to each subunit.
4. Couple neighbouring subunits electrically.
5. Pass each subunit through an output nonlinearity and weight it by the
   ganglion cell's receptive field.

`getSubunitResponses` and `subunitModel` run all five stages.
`subunitLinearResponses` and `getCoupling` expose stages 2 and 4
separately, so the model can be assembled step by step (see
`examples/example_moving_bar.m`).

## Units and coordinate conventions

**Pixels and microns.** The stimulus is a square array of size
`N x N x T` (rows x columns x time). Each pixel is `downsample` microns
wide. Every spatial parameter given in microns (`rfDiameter`,
`subunitDiameter`, `subunitSpacing`, `subunitGridNoiseSD`) is divided by
`downsample` inside the functions and used in pixels.

**Coupling length constant.** Subunit positions, and therefore the
pairwise distances used for coupling, are in pixels. `spaceConstant` is
applied to these pixel distances without conversion. To specify a length
constant in microns, pass `spaceConstant = lambda_um / downsample`, as the
examples do.

**Subunit coordinates.** `centers(:,1)` is the x (column) coordinate and
`centers(:,2)` is the y (row) coordinate of each subunit, in pixels, with
the stimulus occupying `[0, N]` along each axis.

**Stimulus size.** `subunitLinearResponses` requires `N` to be divisible
by 8. All functions assume the stimulus is square.

**Time.** There is no explicit time step. The temporal filter must be
sampled at the same rate as the stimulus frames, and the output has one
sample per stimulus frame.

## 1. Subunit mosaic

`simulateBipolarArray(maxRad, spacing, noiseSD, seed)` builds a
triangular (hexagonal) lattice covering `[0, 2*maxRad]` in both
dimensions:

- columns are spaced by `(sqrt(3)/2) * spacing`;
- points within a column are spaced by `spacing`;
- alternate columns are offset by `spacing/2`.

The nearest-neighbour distance is therefore `spacing`. Each centre is then
jittered by independent Gaussian noise with standard deviation
`2*sqrt(noiseSD)` pixels. The jitter is generated after
`rng(seed, 'twister')`, so a given `positionSeed` always gives the same
mosaic. Note that this call also resets MATLAB's global random number
generator.

With the defaults (32 um spacing, 4 um pixels), a 256 x 256 um field
contains about 68 subunits.

## 2. Linear subunit responses

Each subunit `i` has a Gaussian spatial receptive field centred on its
position, with standard deviation `subunitDiameter/2` (in microns), and
normalized to unit sum:

```
z_i(p) = exp(-|p - c_i|^2 / (2 sigma^2)) / sum_p exp(-|p - c_i|^2 / (2 sigma^2))
```

`subunitLinearResponses` optionally subtracts a broader Gaussian surround
(`subunitSurroundDiameter`, `subunitSurroundWeight`; off by default).

The linear response is the stimulus filtered in space by `z_i` and in time
by the temporal filter `f`, using causal convolution truncated to the
stimulus length:

```
r_i(t) = sum_p z_i(p) * (s(p, .) conv f)(t)
```

This is the projection `r_i(t) = F_i' s_t` of the stimulus window `s_t`
onto the subunit's space-time receptive field `F_i = z_i (x) f`.

**Scale factors.** The convolution is divided by 2 in
`getSubunitResponses`, in `subunitModel`, and in the `'parallel'` and
`'gpu'` methods of `subunitLinearResponses`. It is not divided by 2 in the
default `'normal'` method. `getSubunitResponses` then multiplies the
linear responses by 5, to place them on the input range of the output
nonlinearity.

## 3. Noise

When `noise > 0`, each subunit receives an additive offset drawn from a
Gaussian with standard deviation `noise`. Each subunit gets one draw per
repeat, held constant over the stimulus. The model is run `numIter` times
with independent draws, and the responses are averaged before the output
nonlinearity is applied.

- `getSubunitResponses` uses a fixed list of seeds, so its results are exactly
  reproducible.
- `subunitModel` seeds each repeat from the clock (`'shuffle'`), so its results
  vary from run to run.

The default noise level is `exp(-1.947)`, about 0.143.

## 4. Electrical coupling

Coupling is applied once at each time step, to the noisy linear
responses `R_0`. Each subunit moves toward the responses of its
neighbours, by an amount that falls off exponentially with distance:

```
R_i(t) = R_0i(t) + sum_{j != i} g_ij * (R_0j(t) - R_0i(t)) * exp(-d_ij / lambda)
```

- `d_ij` is the distance between subunit centres (pixels).
- `lambda` is `spaceConstant` (pixels).
- `g_ij = gjGainIn` when `R_0j > R_0i`, so current flows into subunit `i`.
- `g_ij = gjGainOut` when `R_0j < R_0i`, so current flows out of subunit `i`.

With `gjGainIn = gjGainOut = g`, this reduces to a single coupling gain.
In that symmetric case coupling conserves the summed response across
subunits; it redistributes signal without creating or removing it.

The coupling step is a single update, not an iteration to steady state.
`getCoupling` implements one time step. `getSubunitResponses` and
`subunitModel` apply the same update inline.

`computeSubunitCoupling.c` is a MEX implementation of the same update with
a single gain, operating on a whole `subunits x time` array:

```matlab
out = computeSubunitCoupling(currents, g * exp(-D / lambda));
```

## 5. Output nonlinearity and pooling

Each coupled subunit response passes through a cumulative-Gaussian
nonlinearity:

```
N(x) = a1 * Phi(a2 * x + a3) + a4
```

`Phi` is the standard normal CDF. The default parameters are
`nlParams = [1.0204 0.4756 0.2307 -0.0204]`.

The ganglion cell's receptive field is a Gaussian with standard deviation
`rfDiameter/2` (in microns), centred on the stimulus. Each subunit's
weight `W_i` is the value of this Gaussian at the subunit's centre; only
subunits with nonzero weight are simulated. The model returns each
subunit's output multiplied by `W_i / sum(W)`, so the ganglion cell's
excitatory drive is the sum over subunits:

```matlab
rgcDrive = squeeze(sum(nlResponse, 3));   % getSubunitResponses: [contrast x gain x time]
rgcDrive = sum(nlResponse, 1);            % subunitModel:        [1 x time]
```

## Default parameters

| Parameter | `getSubunitResponses` | `subunitModel` | Units |
|---|---|---|---|
| `rfDiameter` | 150 | 200 | um |
| `subunitDiameter` | 32 | 32 | um |
| `subunitSpacing` | 32 | 32 | um |
| `subunitGridNoiseSD` | 2 | 2 | um (see section 1) |
| `downsample` | 4 | 4 | um per pixel |
| `spaceConstant` | 60 | 36.4 | pixels |
| `gjGainIn`, `gjGainOut` | swept (see below) | 0.1 | fraction |
| `noise` | exp(-1.947) | exp(-1.947) | response units |
| `numIter` | 24 | 20 | repeats |
| `nlParams` | [1.0204 0.4756 0.2307 -0.0204] | same | |

`getSubunitResponses` sweeps two parameters internally:

- the contrast multiplier `[0.05, 0.1:0.1:1]`;
- the coupling gain `[0.05 0.10 0.15 0.20 0.25 0.30 0.40 0.50]`, applied to both
  `gjGainIn` and `gjGainOut`.

Values of `gjGainIn` and `gjGainOut` passed to it are ignored.
