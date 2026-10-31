function [subunitResponse, sCenters] = subunitLinearResponses(stimulus, temporalFilter, varargin)

ip = inputParser();
ip.addParameter('subunitDiameter',32.0,@(x)isfloat(x));
ip.addParameter('subunitSurroundDiameter',256,@(x)isfloat(x));
ip.addParameter('subunitSurroundWeight',0.0,@(x)isfloat(x));
ip.addParameter('downsample',4,@(x)isfloat(x));
ip.addParameter('subunitSpacing',32.0,@(x)isfloat(x));
ip.addParameter('subunitGridNoiseSD',2.0,@(x)isfloat(x));
ip.addParameter('positionSeed',1,@(x)isfloat(x) || @(x)isinteger(x));
ip.addParameter('computeMethod','normal',@(x)ischar(x));
ip.parse(varargin{:});

% Pull the results.
downsample = ip.Results.downsample; % Spatial downsample.

% Subunit array parameters.
subunitDiameter = ip.Results.subunitDiameter;
subunitSurroundDiameter = ip.Results.subunitSurroundDiameter;
subunitSurroundWeight = ip.Results.subunitSurroundWeight;
subunitSpacing = ip.Results.subunitSpacing / downsample;
subunitGridNoiseSD = ip.Results.subunitGridNoiseSD / downsample;
positionSeed = round(double(ip.Results.positionSeed));
computeMethod = ip.Results.computeMethod;

% Get the stimulus size.
stimSize = size(stimulus);

% Get then number of stimulus dimensions (third dimension is time).
nd = ndims(stimulus);

% If there is no temporal component, it doesn't make sense to run in
% parallel.
if nd < 3
    computeMethod = '2d';
end

sz = stimSize(1)/2;
% Make sure it's evenly divisible by 2.
if mod(sz/2,2) ~= 0
    disp('Size must be divisible by two.')
    return;
end
% Make the subunit array.
sCenters = round(simulateBipolarArray(sz, subunitSpacing,...
    subunitGridNoiseSD, positionSeed));

% Create the subunit receptive field template.
[x,y] = meshgrid(linspace(-stimSize(1)/2,stimSize(1)/2,stimSize(1)));
a = 1/(2*(subunitDiameter/2/downsample)^2);
as = 1/(2*(subunitSurroundDiameter/2/downsample)^2);

% Define a Gaussian function for RF calculations.
gauss = @(p,x)(p(1)*exp(-(x.^2/(2*p(2)^2))));

% Get the subunit responses.
numSubs = size(sCenters,1);

switch computeMethod
    case 'gpu'
        getConvolutionsGPU();
    case 'parallel'
        getConvolutionsParallel();
    case '2d'
        getConvolutions2D();
    otherwise
        getConvolutionsNoParallel();
end

function getConvolutions2D
    subunitResponse = zeros(numSubs,1);
    for j=1:numSubs
       curX = sCenters(j,1);
       curY = sCenters(j,2);

       % Convolve the stimulus with the filter.
       % Spatial profile.
       zc = exp(-(a*(x-(curX-size(x,1)/2)).^2 + a*(y-(curY-size(x,2)/2)).^2));
       if subunitSurroundWeight > 0
           zs = subunitSurroundWeight*exp(-(as*(x-(curX-size(x,1)/2)).^2 + as*(y-(curY-size(x,2)/2)).^2)) ...
               / (subunitSurroundDiameter/subunitDiameter)^2;
       else
           zs = 0;
       end
       z = (zc - zs) / sum(zc(:)); % Normalize
       
       % Convolution is just the sum of the Hadamard product in the 2-d
       % case.
       tmp = z .* stimulus;
       subunitResponse(j) = sum(tmp(:));
    end
end

function getConvolutionsNoParallel
    subunitResponse = zeros(numSubs, stimSize(3));
    for j=1:numSubs
       curX = sCenters(j,1);
       curY = sCenters(j,2);

       % Convolve the stimulus with the filter.
       % Spatial profile.
       zc = exp(-(a*(x-(curX-size(x,1)/2)).^2 + a*(y-(curY-size(x,2)/2)).^2));
       if subunitSurroundWeight > 0
           zs = subunitSurroundWeight*exp(-(as*(x-(curX-size(x,1)/2)).^2 + as*(y-(curY-size(x,2)/2)).^2)) ...
               / (subunitSurroundDiameter/subunitDiameter)^2;
       else
           zs = 0;
       end
       z = (zc - zs) / sum(zc(:)); % Normalize

       for k = 1 : stimSize(1)
           for m = 1 : stimSize(2)
               tmp = z(k,m)*conv(squeeze(stimulus(k,m,:))',temporalFilter(:)','full');
               subunitResponse(j,:) = subunitResponse(j,:) + tmp(1:stimSize(3));
           end
       end
    end
end

function getConvolutionsParallel
    subunitResponse = zeros(numSubs, stimSize(3));
    parfor j=1:numSubs
       curX = sCenters(j,1);
       curY = sCenters(j,2);

       % Convolve the stimulus with the filter.
       % Spatial profile.
       zc = exp(-(a*(x-(curX-size(x,1)/2)).^2 + a*(y-(curY-size(x,2)/2)).^2));
       if subunitSurroundWeight > 0
           zs = subunitSurroundWeight*exp(-(as*(x-(curX-size(x,1)/2)).^2 + as*(y-(curY-size(x,2)/2)).^2)) ...
               / (subunitSurroundDiameter/subunitDiameter)^2;
       else
           zs = 0;
       end
       z = (zc - zs) / sum(zc(:)); % Normalize

       for k = 1 : stimSize(1)
           for m = 1 : stimSize(2)
               tmp = z(k,m)*conv(squeeze(stimulus(k,m,:))',temporalFilter(:)','full')/2;
               subunitResponse(j,:) = subunitResponse(j,:) + tmp(1:stimSize(3));
           end
       end
    end
end

function getConvolutionsGPU
    subunitResponse = gpuArray.zeros(numSubs, stimSize(3));
    x = gpuArray(x);
    y = gpuArray(y);
    temporalFilter = gpuArray(temporalFilter);
    for j=1:numSubs
       curX = sCenters(j,1);
       curY = sCenters(j,2);

       % Convolve the stimulus with the filter.
       % Spatial profile.
       zc = exp(-(a*(x-(curX-size(x,1)/2)).^2 + a*(y-(curY-size(x,2)/2)).^2));
       if subunitSurroundWeight > 0
           zs = subunitSurroundWeight*exp(-(as*(x-(curX-size(x,1)/2)).^2 + as*(y-(curY-size(x,2)/2)).^2)) ...
               / (subunitSurroundDiameter/subunitDiameter)^2;
       else
           zs = 0;
       end
       z = (zc - zs) / sum(zc(:)); % Normalize

       for k = 1 : stimSize(1)
           for m = 1 : stimSize(2)
               tmp = z(k,m)*conv(squeeze(stimulus(k,m,:))',temporalFilter(:)','full')/2;
               subunitResponse(j,:) = subunitResponse(j,:) + tmp(1:stimSize(3));
           end
       end
    end
    subunitResponse = gather(subunitResponse);
end


end
