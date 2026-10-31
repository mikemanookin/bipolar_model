# Examples

Both scripts add `../src` to the MATLAB path, so they can be run from any
folder. They require the Image Processing Toolbox and the Statistics and
Machine Learning Toolbox.

## example_moving_bar.m

This script builds the model one stage at a time, using a 40 um bar
moving at 400 um/s across a 256 um field. It uses
`subunitLinearResponses` for the linear responses and `getCoupling` for
coupling; the nonlinearity and pooling are written out in the script.

It compares four versions of the model:

- linear bipolar output, without coupling;
- linear bipolar output, with coupling;
- nonlinear bipolar output, without coupling;
- nonlinear bipolar output, with coupling.

The figure shows:

- the subunit mosaic, coloured by each subunit's ganglion cell pooling weight;
- the linear subunit responses, ordered by position along the direction of motion;
- the pooled ganglion cell drive for each version of the model.

Runtime is about 30 s.

Start here to see how the stages fit together, or to modify one stage (for
example, the nonlinearity) without touching the rest.

## example_coupling_sweep.m

This script runs `getSubunitResponses` on the same moving bar. It
simulates the full model, with noise, across the function's built-in grid
of 11 contrasts and 8 coupling gains.

The figure shows:

- the pooled ganglion cell drive over time at full contrast, for each coupling gain;
- the peak pooled drive as a function of contrast, for each coupling gain.

Runtime is one to two minutes.

## Note on the coupling length constant

Distances between subunits are computed in pixels. Both examples
therefore convert the 36 um coupling length constant to pixels
(`36 / umPerPixel`) before passing it to the model. See
[../docs/model.md](../docs/model.md#units-and-coordinate-conventions).
