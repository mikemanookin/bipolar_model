# Bipolar cell subunit model

A MATLAB model of the array of diffuse bipolar cells that provides
excitatory input to parasol and smooth monostratified ganglion cells in
the primate retina. The model was used in our study of predictive motion
encoding in primate retinal ganglion cells [CITATION].

## Author and terms of use

This code was written by Mike Manookin, Associate Professor at the
University of Washington.

**This code may not be used in publications without the written approval
of the author.**

## Overview

The model represents a ganglion cell's excitatory input as a mosaic of
bipolar cell subunits. It has five stages:

1. **Subunit mosaic.** Subunits are placed on a jittered hexagonal lattice
   (32 µm mean spacing by default).
2. **Linear filtering.** Each subunit filters the stimulus through a
   Gaussian spatial receptive field and a temporal filter.
3. **Noise.** Each subunit receives additive noise. The model is run
   repeatedly and the responses are averaged.
4. **Electrical coupling.** Each subunit's response moves toward its
   neighbours' responses, by an amount that falls off exponentially with
   distance.
5. **Output and pooling.** Each subunit passes through a static output
   nonlinearity. The outputs are then weighted by the ganglion cell's
   Gaussian receptive field and summed.

The functions can be run as a complete model (`getSubunitResponses`,
`subunitModel`) or assembled one stage at a time
(`subunitLinearResponses`, `getCoupling`).

## Documentation

| Document | Contents |
|---|---|
| [docs/README.md](docs/README.md) | Requirements, quick start, repository layout and minimal usage |
| [docs/model.md](docs/model.md) | Each stage of the model, the equations as implemented, units and coordinate conventions, default parameters |
| [docs/api.md](docs/api.md) | Reference for every function: inputs, outputs and defaults |
| [examples/README.md](examples/README.md) | Description of the example scripts |
| [LICENSE](LICENSE) | License |

## Quick start

Requires MATLAB R2019b or later, the Image Processing Toolbox and the
Statistics and Machine Learning Toolbox.

```matlab
% From the repository root
run('examples/example_moving_bar.m');       % model built step by step (~30 s)
run('examples/example_coupling_sweep.m');   % sweep of contrast and coupling gain (1-2 min)
```

## Repository layout

```
README.md                this file
LICENSE
docs/                    documentation
examples/                runnable example scripts
src/coupling/            model source code
```

## Citation

If you use this code with the author's approval, please cite:

[CITATION]
