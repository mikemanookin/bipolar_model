%% Example: sweep of coupling gain and stimulus contrast (getSubunitResponses)
%
% getSubunitResponses runs the full bipolar subunit model (linear
% subunit filtering, additive noise, electrical coupling, output
% nonlinearity and ganglion cell pooling weights) over a built-in grid of
% stimulus contrasts and coupling gains:
%
%   contrast multipliers: 0.05, 0.1, 0.2, ..., 1.0   (11 values)
%   coupling gains:       0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50
%
% The grid is hard-coded inside the function, so the 'gjGainIn' and
% 'gjGainOut' parameters are not used by it.
%
% Noise is drawn from fixed seeds, so repeated runs give identical
% results.
%
% Requirements: MATLAB with the Image Processing and Statistics and
% Machine Learning Toolboxes. Run from any folder; the script adds
% ../src to the path itself. Runtime is one to two minutes on a laptop.

clear; close all;

%% Paths
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'src')));

%% Stimulus: a bar moving across a square field
umPerPixel = 4;      % microns per pixel (the model's 'downsample' parameter)
stimSizeUm = 256;
dt         = 0.01;   % s
numFrames  = 120;

nPix = stimSizeUm / umPerPixel;
barWidthPix = 40 / umPerPixel;
barSpeedPix = 4 / umPerPixel;           % 4 um per time step
stimulus = zeros(nPix, nPix, numFrames);
for t = 1:numFrames
    leftEdge = -barWidthPix + (t-1) * barSpeedPix;
    cols = (1:nPix) > leftEdge & (1:nPix) <= leftEdge + barWidthPix;
    stimulus(:, cols, t) = 1;
end

%% Temporal filter (biphasic)
tFilt = (0:dt:0.3-dt)';
lobe = @(tau) (tFilt/tau).^4 .* exp(-tFilt/tau);
temporalFilter = lobe(0.012)/max(lobe(0.012)) - 0.6*lobe(0.020)/max(lobe(0.020));
temporalFilter = temporalFilter / sum(abs(temporalFilter));

%% Run the model
% Subunit positions and pairwise distances are computed in pixels inside
% getSubunitResponses, so 'spaceConstant' is given in pixels here
% (36 um / 4 um per pixel = 9 pixels).
lambdaUm = 36;

tic;
[nlResponse, linResponse, W] = getSubunitResponses(stimulus, temporalFilter, ...
    'downsample',      umPerPixel, ...
    'rfDiameter',      200, ...            % RGC receptive-field diameter (um)
    'subunitDiameter', 32, ...             % bipolar receptive-field diameter (um)
    'subunitSpacing',  32, ...             % mean bipolar spacing (um)
    'spaceConstant',   lambdaUm / umPerPixel, ...
    'numIter',         8, ...              % noise repeats averaged together
    'positionSeed',    1);
fprintf('getSubunitResponses finished in %.1f s\n', toc);

% nlResponse is [contrast x gain x subunit x time] and is already
% multiplied by the normalized RGC pooling weights, so summing over
% subunits gives the ganglion cell's input.
contrasts = [0.05, 0.1:0.1:1];
gains     = [5 10 15 20 25 30 40 50] / 100;
rgc = squeeze(sum(nlResponse, 3));            % [contrast x gain x time]

% Express each response relative to its value before the bar arrives.
rgc = rgc - rgc(:, :, 1);

time = (0:numFrames-1) * dt * 1000;           % ms

%% Plot
figure('Color', 'w', 'Position', [100 100 1100 420]);
colors = winter(numel(gains));

% Time course at full contrast for each coupling gain
subplot(1, 2, 1); hold on;
for q = 1:numel(gains)
    plot(time, squeeze(rgc(end, q, :)), 'Color', colors(q, :), 'LineWidth', 1.25);
end
xlabel('Time (ms)'); ylabel('Pooled response (re. baseline)');
title('Full contrast');
legend(arrayfun(@(g) sprintf('g = %.2f', g), gains, 'UniformOutput', false), ...
    'Location', 'eastoutside');

% Peak response as a function of contrast for each coupling gain
subplot(1, 2, 2); hold on;
peakResp = max(rgc, [], 3);                   % [contrast x gain]
for q = 1:numel(gains)
    plot(contrasts, peakResp(:, q), '-o', 'Color', colors(q, :), ...
        'MarkerFaceColor', colors(q, :), 'MarkerSize', 4);
end
xlabel('Contrast multiplier'); ylabel('Peak pooled response');
title('Contrast response by coupling gain');
