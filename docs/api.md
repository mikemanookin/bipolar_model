# Function reference

All functions are in `src/coupling`. Spatial parameters are given in
microns unless noted otherwise. See [model.md](model.md) for the
equations and the unit conventions.

---

## subunitLinearResponses

```matlab
[subunitResponse, sCenters] = subunitLinearResponses(stimulus, temporalFilter, Name, Value, ...)
```

Computes the linear response of each bipolar subunit. Each response is
the stimulus filtered by the subunit's Gaussian spatial receptive field
and by the temporal filter. No noise, coupling or nonlinearity is applied.

**Inputs**

| Name | Description |
|---|---|
| `stimulus` | `N x N x T` array (rows x columns x time), or `N x N` for a static image. `N` must be divisible by 8. |
| `temporalFilter` | Temporal filter, one sample per stimulus frame. |

**Name-value parameters**

| Name | Default | Description |
|---|---|---|
| `downsample` | 4 | Microns per pixel. |
| `subunitDiameter` | 32 | Subunit receptive-field diameter (um); the Gaussian SD is half of this. |
| `subunitSurroundDiameter` | 256 | Surround diameter (um). |
| `subunitSurroundWeight` | 0 | Surround strength; 0 disables the surround. |
| `subunitSpacing` | 32 | Mean spacing between subunits (um). |
| `subunitGridNoiseSD` | 2 | Positional jitter parameter (see `simulateBipolarArray`). |
| `positionSeed` | 1 | Seed for the mosaic jitter. |
| `computeMethod` | `'normal'` | `'normal'`, `'parallel'` (uses `parfor`), or `'gpu'`. A 2-D stimulus always uses a single-frame method. |

**Outputs**

| Name | Description |
|---|---|
| `subunitResponse` | `numSubunits x T` linear responses (`numSubunits x 1` for a 2-D stimulus). |
| `sCenters` | `numSubunits x 2` subunit centres in pixels, as `[x y]` (column, row). |

The `'parallel'` and `'gpu'` methods divide the response by 2 and the
`'normal'` method does not, so their outputs differ by a factor of 2.

---

## getCoupling

```matlab
bpCurrents = getCoupling(bpCurrents, bpDistances, spaceConstant, gjGainIn, gjGainOut)
```

Applies one step of electrical coupling to the subunit responses at a
single time point. Each subunit moves toward its neighbours' responses,
weighted by `exp(-d/spaceConstant)`.

| Name | Description |
|---|---|
| `bpCurrents` | `numSubunits x 1` responses at one time point. |
| `bpDistances` | `numSubunits x numSubunits` pairwise distances, e.g. `squareform(pdist(sCenters))`. |
| `spaceConstant` | Length constant, in the same units as `bpDistances`. |
| `gjGainIn` | Gain applied when a neighbour's response is larger (inward current). |
| `gjGainOut` | Gain applied when a neighbour's response is smaller (outward current). |

Returns the coupled responses (`numSubunits x 1`). The trailing
`modelfun` and `nlParams` arguments are accepted but not used.

To couple a whole time series, loop over the columns:

```matlab
D = squareform(pdist(sCenters));
for t = 1:size(R, 2)
    R(:, t) = getCoupling(R(:, t), D, lambdaPix, g, g);
end
```

---

## getSubunitResponses

```matlab
[nlResponse, response, W] = getSubunitResponses(stimulus, temporalFilter, Name, Value, ...)
```

Runs the full model (linear filtering, noise, coupling, nonlinearity and
pooling weights) over a built-in grid of 11 contrast multipliers and 8
coupling gains (see [model.md](model.md#default-parameters)).

**Name-value parameters**

| Name | Default | Description |
|---|---|---|
| `rfDiameter` | 150 | Ganglion cell receptive-field diameter (um). Below 30, only the 12 subunits nearest the centre are kept. |
| `subunitDiameter` | 32 | Subunit receptive-field diameter (um). |
| `downsample` | 4 | Microns per pixel. |
| `nlParams` | `[1.0204 0.4756 0.2307 -0.0204]` | Parameters of the cumulative-Gaussian output nonlinearity. |
| `noise` | `exp(-1.947)` | SD of the additive subunit noise. |
| `spaceConstant` | 60 | Coupling length constant, in **pixels**. |
| `gjGainIn`, `gjGainOut` | 0.1 | Not used; the gain is swept internally. |
| `subunitSpacing` | 32 | Mean subunit spacing (um). |
| `subunitGridNoiseSD` | 2 | Positional jitter parameter. |
| `positionSeed` | 1 | Seed for the mosaic jitter. |
| `numIter` | 24 | Number of noise repeats averaged together. Fixed seeds are used, so results are reproducible. |

**Outputs**

| Name | Size | Description |
|---|---|---|
| `nlResponse` | `11 x 8 x numSubunits x T` | Output after the nonlinearity, multiplied by the normalized pooling weights `W/sum(W)`. Sum over dimension 3 to get the ganglion cell's drive. |
| `response` | `11 x 8 x numSubunits x T` | Coupled responses before the nonlinearity (averaged over noise repeats). |
| `W` | `numSubunits x 1` | Unnormalized ganglion cell pooling weights. |

---

## subunitModel

```matlab
[nlResponse, response] = subunitModel(stimulus, temporalFilter, Name, Value, ...)
```

The same model as `getSubunitResponses`, run at a single contrast and
coupling gain (`gjGainIn`, `gjGainOut`; default 0.1 each). It also
differs in these defaults:

- `rfDiameter` is 200;
- `spaceConstant` is 36.4 pixels;
- `numIter` is 20;
- the noise is seeded from the clock, so results vary between runs.

`nlResponse` is `numSubunits x T` and already weighted by `W/sum(W)`.

This function calls `models.coupling.simulateBipolarArray` and
`models.ln.outputNonlinearity`. These come from the lab's `+models`
package, which is not included in this repository.

---

## simulateBipolarArray

```matlab
bpCenters = simulateBipolarArray(maxRad, bipolarSpacing, noiseSD, positionSeed)
```

Returns the `[x y]` centres (pixels) of a jittered hexagonal mosaic
covering `[0, 2*maxRad]` in each dimension. Neighbouring centres are
`bipolarSpacing` apart before jitter. Each coordinate is jittered by
Gaussian noise with SD `2*sqrt(noiseSD)`. The function calls
`rng(positionSeed, 'twister')`, which resets MATLAB's global random
number generator.

---

## getBipolarActivation

```matlab
[bcWeights, bcActivation] = getBipolarActivation(bpCenters, rgcRFDiam, bpRFDiam, numBars)
```

For a set of bipolar centres, returns two quantities:

- `bcWeights` (`1 x numBipolars`): the overlap of each bipolar's Gaussian
  receptive field with the ganglion cell's receptive field;
- `bcActivation` (`numBars x numBipolars`): the fraction of each bipolar's
  receptive field covered by each of `numBars` equal-width bars spanning
  the ganglion cell's receptive field.

The function uses a 4 um sampling grid.

---

## makeSubunitHexGrid

```matlab
[posX, posY] = makeSubunitHexGrid(stimSize, subunit_RF, subunit_spacing, noiseSD)
```

An alternative hexagonal grid generator. It returns mean-centred,
rounded coordinates. The jitter uses the global random number generator
without a seed.

---

## computeSubunitCoupling (MEX)

```matlab
out = computeSubunitCoupling(currents, couplingWeights)
```

A C implementation of the coupling step applied to every time point at
once.

| Name | Description |
|---|---|
| `currents` | `numSubunits x T` array. |
| `couplingWeights` | `numSubunits x numSubunits` matrix, `g * exp(-D / lambda)`. |

To build it, run `mex computeSubunitCoupling.c` in `src/coupling`. The
repository's `.gitignore` excludes compiled MEX files.
