%% Example: coupled bipolar subunit array responding to a moving bar
%
% This example builds the model step by step, so each stage of the
% computation is visible:
%
%   1. Make a moving-bar stimulus on a square pixel grid.
%   2. Compute the linear response of each bipolar subunit
%      (spatial Gaussian x temporal filter; subunitLinearResponses).
%   3. Apply electrical coupling between neighbouring subunits
%      (getCoupling).
%   4. Pass each subunit through a sigmoidal output nonlinearity.
%   5. Pool the subunit outputs with the ganglion cell's Gaussian
%      receptive-field weights.
%
% It compares four versions of the model: linear vs. nonlinear bipolar
% output, each with and without coupling.
%
% Requirements: MATLAB with the Image Processing and Statistics and
% Machine Learning Toolboxes. Run from any folder; the script adds
% ../src to the path itself.

clear; close all;

%% Paths
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'src')));

%% Spatial and temporal sampling
umPerPixel = 4;      % microns per pixel (passed to the model as 'downsample')
stimSizeUm = 256;    % width and height of the stimulus (microns)
dt         = 0.01;   % time step (s)
numFrames  = 120;    % stimulus duration in time steps

nPix = stimSizeUm / umPerPixel;   % 64 pixels; must be divisible by 8
assert(mod(nPix, 8) == 0, 'Stimulus width (pixels) must be divisible by 8.');

%% 1. Moving-bar stimulus (contrast units; background = 0)
barWidthUm  = 40;    % bar width (microns)
barSpeedUm  = 4;     % bar speed (microns per time step) -> 400 um/s at dt = 10 ms
barContrast = 1;

barWidthPix = barWidthUm / umPerPixel;
startPix    = -barWidthPix;             % bar starts just off the left edge
stimulus    = zeros(nPix, nPix, numFrames);
for t = 1:numFrames
    leftEdge = startPix + (t-1) * barSpeedUm / umPerPixel;
    cols = (1:nPix) > leftEdge & (1:nPix) <= leftEdge + barWidthPix;
    stimulus(:, cols, t) = barContrast;   % bar spans all rows; moves along x (columns)
end

%% Temporal filter (biphasic, Off-lobe after an On-lobe)
tFilt = (0:dt:0.3-dt)';
tau1 = 0.012; tau2 = 0.020; n = 4;
lobe = @(tau) (tFilt/tau).^n .* exp(-tFilt/tau);
temporalFilter = lobe(tau1)/max(lobe(tau1)) - 0.6*lobe(tau2)/max(lobe(tau2));
temporalFilter = temporalFilter / sum(abs(temporalFilter));

%% 2. Linear subunit responses
subunitDiameterUm = 32;   % bipolar receptive-field diameter (microns)
subunitSpacingUm  = 32;   % mean bipolar spacing (microns)

% centers(:,1) is the x (column) coordinate and centers(:,2) the y (row)
% coordinate of each subunit, in pixels.
[linearResp, centers] = subunitLinearResponses(stimulus, temporalFilter, ...
    'downsample',        umPerPixel, ...
    'subunitDiameter',   subunitDiameterUm, ...
    'subunitSpacing',    subunitSpacingUm, ...
    'subunitGridNoiseSD', 2.0, ...
    'positionSeed',      1, ...
    'computeMethod',     'normal');
numSubunits = size(linearResp, 1);

% Scale so that the strongest linear response is 1. This puts the drive
% in the operating range of the output nonlinearity below.
linearResp = linearResp / max(abs(linearResp(:)));

%% Ganglion cell receptive field: Gaussian pooling weights
rgcDiameterUm = 200;                       % RGC receptive-field diameter (microns)
rgcSigmaPix   = rgcDiameterUm / 2 / umPerPixel;
rfCenter      = [nPix nPix] / 2;
d2   = sum((centers - rfCenter).^2, 2);
W    = exp(-d2 / (2 * rgcSigmaPix^2));
W    = W / sum(W);

%% 3. Electrical coupling
% Subunit centres are in pixels, so the length constant must be converted
% from microns to pixels before it is applied to the pairwise distances.
lambdaUm  = 36;                      % coupling length constant (microns)
lambdaPix = lambdaUm / umPerPixel;
gjGain    = 0.1;                     % coupling gain (fraction of the difference shared)
distances = squareform(pdist(centers));

coupledResp = zeros(size(linearResp));
for t = 1:numFrames
    coupledResp(:, t) = getCoupling(linearResp(:, t), distances, lambdaPix, gjGain, gjGain);
end

%% 4. Output nonlinearity (cumulative-Gaussian; same form as getSubunitResponses)
nlParams = [1.0204 0.4756 0.2307 -0.0204];
outputNL = @(a, x) a(1) * normcdf(a(2)*x + a(3), 0, 1) + a(4);
% Scale the unit-normalized drive onto the steep part of the nonlinearity
% (getSubunitResponses applies a similar fixed scale factor internally).
driveGain = 5;
nlUncoupled = outputNL(nlParams, driveGain * linearResp);
nlCoupled   = outputNL(nlParams, driveGain * coupledResp);

% Remove the baseline output to the blank background, so that every
% model reports the change in output.
nlBaseline  = outputNL(nlParams, 0);
nlUncoupled = nlUncoupled - nlBaseline;
nlCoupled   = nlCoupled   - nlBaseline;

%% 5. Pool into the ganglion cell
rgc.linearUncoupled    = W' * linearResp;
rgc.linearCoupled      = W' * coupledResp;
rgc.nonlinearUncoupled = W' * nlUncoupled;
rgc.nonlinearCoupled   = W' * nlCoupled;

%% Plot
time = (0:numFrames-1) * dt * 1000;   % ms

figure('Color', 'w', 'Position', [100 100 1100 750]);

% Subunit mosaic and RGC weights
subplot(2, 2, 1);
scatter(centers(:,1) * umPerPixel, centers(:,2) * umPerPixel, 60, W, 'filled');
hold on;
th = linspace(0, 2*pi, 200);
plot(rfCenter(1)*umPerPixel + rgcDiameterUm/2*cos(th), ...
     rfCenter(2)*umPerPixel + rgcDiameterUm/2*sin(th), 'k--');
axis xy equal; xlim([0 stimSizeUm]); ylim([0 stimSizeUm]);
xlabel('x (\mum)'); ylabel('y (\mum)');
title(sprintf('%d bipolar subunits; colour = RGC weight', numSubunits));
colorbar;

% Linear subunit responses, ordered by position along the motion axis
subplot(2, 2, 2);
[~, order] = sort(centers(:, 1));
imagesc(time, centers(order, 1) * umPerPixel, linearResp(order, :));
axis xy;
xlabel('Time (ms)'); ylabel('Subunit x position (\mum)');
title('Linear subunit responses');
colorbar;

% RGC output: linear bipolar output
subplot(2, 2, 3);
plot(time, rgc.linearUncoupled, 'k', time, rgc.linearCoupled, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Pooled response');
legend({'Uncoupled', 'Coupled'}, 'Location', 'best');
title('Linear bipolar output');

% RGC output: rectified bipolar output
subplot(2, 2, 4);
plot(time, rgc.nonlinearUncoupled, 'k', time, rgc.nonlinearCoupled, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Pooled response');
legend({'Uncoupled', 'Coupled'}, 'Location', 'best');
title('Nonlinear bipolar output');

%% Summary
[~, iLin] = max(rgc.linearUncoupled);
[~, iNL]  = max(rgc.nonlinearUncoupled);
[~, iNLc] = max(rgc.nonlinearCoupled);
fprintf('Peak time, linear output, uncoupled:     %6.0f ms\n', time(iLin));
fprintf('Peak time, nonlinear output, uncoupled:  %6.0f ms\n', time(iNL));
fprintf('Peak time, nonlinear output, coupled:    %6.0f ms\n', time(iNLc));
