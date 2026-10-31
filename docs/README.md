# Bipolar cell subunit model

This repository contains the MATLAB model of the diffuse bipolar cell
array that provides excitatory input to parasol and smooth monostratified
ganglion cells in the primate retina. The model consists of a jittered
hexagonal mosaic of bipolar subunits. Each subunit filters the stimulus in
space and time, receives additive noise, is electrically coupled to its
neighbours, and passes its output through a static nonlinearity. The
subunit outputs are then weighted by the ganglion cell's Gaussian
receptive field.

## Contents

| Document | What it covers |
|---|---|
| [model.md](model.md) | The model's stages, the equations each function implements, and its units and coordinate conventions |
| [api.md](api.md) | Reference for every function in `src/coupling`: inputs, outputs and defaults |
| [../examples/README.md](../examples/README.md) | The example scripts and what they show |

## Requirements

- MATLAB R2019b or later.
- Image Processing Toolbox (`fspecial`).
- Statistics and Machine Learning Toolbox (`normcdf`, `pdist`, `squareform`).
- Parallel Computing Toolbox, only for `computeMethod = 'parallel'` or `'gpu'` in
  `subunitLinearResponses`.

## Quick start

```matlab
% From the repository root
addpath(genpath('src'));
run('examples/example_moving_bar.m');       % step-by-step model, ~30 s
run('examples/example_coupling_sweep.m');   % contrast x coupling-gain sweep, 1-2 min
```

The example scripts add `src` to the path themselves, so they can also be
run directly from the `examples` folder.

## Repository layout

```
src/coupling/
    simulateBipolarArray.m      jittered hexagonal mosaic of subunit centres
    subunitLinearResponses.m    linear (spatial x temporal) response of each subunit
    getCoupling.m               one step of electrical coupling between subunits
    getSubunitResponses.m       full model, swept over contrast and coupling gain
    subunitModel.m              full model at a single parameter setting
    getBipolarActivation.m      overlap of bipolar receptive fields with bars
    makeSubunitHexGrid.m        alternative hexagonal grid generator
    computeSubunitCoupling.c    MEX implementation of the coupling step
examples/
    example_moving_bar.m
    example_coupling_sweep.m
docs/
```

## Minimal use

```matlab
addpath(genpath('src'));

umPerPixel = 4;                               % spatial sampling (microns per pixel)
stimulus   = zeros(64, 64, 100);              % [y x time], 256 x 256 um
stimulus(:, 20:30, 30:60) = 1;                % a bar flashed for 30 time steps
tf = [0 0.2 0.6 1 0.6 0.1 -0.2 -0.3 -0.2 -0.1]';   % temporal filter, one tap per time step

% Linear subunit responses (subunits x time) and subunit centres (pixels)
[lin, centers] = subunitLinearResponses(stimulus, tf, 'downsample', umPerPixel);

% Couple neighbouring subunits (length constant 36 um = 9 pixels)
D = squareform(pdist(centers));
coupled = zeros(size(lin));
for t = 1:size(lin, 2)
    coupled(:, t) = getCoupling(lin(:, t), D, 36/umPerPixel, 0.1, 0.1);
end
```

See [model.md](model.md) for the remaining stages (nonlinearity and
pooling), and the examples for complete, plotted versions.

## License

MIT; see [LICENSE](../LICENSE).
