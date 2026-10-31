function [nlResponse, response, W] = getSubunitResponses(stimulus, temporalFilter, varargin)

ip = inputParser();
ip.addParameter('rfDiameter',150.0,@(x)isfloat(x));
ip.addParameter('subunitDiameter',32.0,@(x)isfloat(x));
ip.addParameter('downsample',4,@(x)isfloat(x));
ip.addParameter('nlParams',[1.0204 0.4756 0.2307 -0.0204],@(x)isfloat(x));
ip.addParameter('noise',exp(-1.9470),@(x)isfloat(x));
ip.addParameter('spaceConstant',60,@(x)isfloat(x));
ip.addParameter('gjGainOut',0.1,@(x)isfloat(x));
ip.addParameter('gjGainIn',0.1,@(x)isfloat(x));
ip.addParameter('subunitSpacing',32.0,@(x)isfloat(x));
ip.addParameter('subunitGridNoiseSD',2.0,@(x)isfloat(x));
ip.addParameter('positionSeed',1,@(x)isfloat(x) || @(x)isinteger(x));
ip.addParameter('numIter',24,@(x)isfloat(x) || @(x)isinteger(x));
ip.parse(varargin{:});

% Pull the results.
rfDiameter = ip.Results.rfDiameter;
subunitDiameter = ip.Results.subunitDiameter;
downsample = ip.Results.downsample; % Spatial downsample.
nlParams = ip.Results.nlParams;
noise = ip.Results.noise;
spaceConstant = ip.Results.spaceConstant; %space constant for effective gap junction spread (microns)
gjGainOut = ip.Results.gjGainOut; % fraction of signal leaving bipolar via gap junctions
gjGainIn = ip.Results.gjGainIn; % fraction of signal gained from neighboring bipolars via gap junctions

% Subunit array parameters.
subunitSpacing = ip.Results.subunitSpacing / downsample;
subunitGridNoiseSD = ip.Results.subunitGridNoiseSD / downsample;
positionSeed = round(double(ip.Results.positionSeed));
numIter = round(double(ip.Results.numIter)); % Number of iterations for adding noise.

% Use set seeds to make sure that results are exactly replecated.
seeds = [5299,4189,8929,3320,8213,417,1077,5951,5299,4189,3355,6226,4382,7359,5181,5789,6454,9903,8199,4133,8763,8238,8478,3243];
if numIter > length(seeds)
    nseeds = numIter - length(seeds);
    noiseStream = RandStream('mt19937ar', 'Seed', 2048);
    seeds(end+(1:nseeds)) = round(noiseStream.rand(1,nseeds)*10000);
end

% Define the nonlinear function.
nlfun = @(a, x)(a(1)*normcdf(a(2)*x+a(3), 0, 1) + a(4));
% Get the stimulus size.
stimSize = size(stimulus);


% Make the subunit array.
sCenters = round(simulateBipolarArray(stimSize(1)/2, subunitSpacing,...
    subunitGridNoiseSD, positionSeed));

% Determine the number of subunits to keep.
if rfDiameter < 30
    numSubunits = 12;
    % Get the Euclidean distance to the RF center.
    myRad = sqrt((sCenters(:,1)-stimSize(1)/2).^2 + (sCenters(:,2)-stimSize(1)/2).^2);
    % Keep only the cells nearest the RF center to speed processing.
    [~,idx] = sort(myRad);
    sCenters = sCenters(idx(1:numSubunits),:);
end

% pair-wise distances
sDistances = squareform(pdist(sCenters));

% Make the subunit spatial filter.
% sWinSize = 4*subunitDiameter / downsample;
% sGaussWin = fspecial('gaussian', round(sWinSize), subunitDiameter/2/downsample);

% Make the RGC space filter.
RGCgaussWin = fspecial('gaussian', stimSize(1), rfDiameter/2/downsample);

% Create the subunit receptive field template.
[x,y] = meshgrid(linspace(-stimSize(1)/2,stimSize(1)/2,stimSize(1)));
a = 1/(2*(subunitDiameter/2/downsample)^2);

% Get the subunit weights
numSubs = size(sCenters,1);
W = zeros(numSubs, 1);
subunitResponse = zeros(numSubs, stimSize(3));
for i=1:numSubs
   curX = sCenters(i,1);
   curY = sCenters(i,2);
   
   if curX > 0 && curX < stimSize(1)  && curY > 0 && curY < stimSize(2)
       W(i) = RGCgaussWin(curX, curY);
   end
   
   % Convolve the stimulus with the filter.
   if W(i) > 0
       % Spatial profile.
       z = exp(-(a*(x-(curX-size(x,1)/2)).^2 + a*(y-(curY-size(x,2)/2)).^2));
       z = z / sum(z(:)); % Normalize
       
       for k = 1 : stimSize(1)
           for m = 1 : stimSize(2)
               tmp = z(k,m)*conv(squeeze(stimulus(k,m,:))',temporalFilter(:)','full')/2;
               subunitResponse(i,:) = subunitResponse(i,:) + tmp(1:stimSize(3));
           end
       end
   end
end

% Make sure the responses scale properly with the x-axis of the
% nonlinearity.
% max(abs(subunitResponse(:)))
subunitResponse = subunitResponse * 5;

if noise <= 0
    numIter = 1;
end

ctM = [0.05,0.1:0.1:1];
gainValues = [5 10 15 20 25 30 40 50]/100;

response = zeros([length(ctM),length(gainValues),size(subunitResponse)]);
for k = 1 : numIter
    for p = 1 : length(ctM)
        for q = 1 : length(gainValues)
            gjGainIn = gainValues(q);
            gjGainOut = gainValues(q);
            if noise > 0
                % Seed the random number generator with a random seed.
                ns = RandStream('mt19937ar','Seed',seeds(k));
                addedNoise = noise * ns.randn(numSubs,1);
                % Add the noise.
                subunitCurrent = ctM(p)*subunitResponse + addedNoise*ones(1,size(subunitResponse,2));
            else
                subunitCurrent = ctM(p)*subunitResponse;
            end

            % Junctional coupling
            if gjGainOut > 0
                

                for m = 1 : size(subunitCurrent,2)
                    % Pair-wise subunit voltage differences
                    drivingForce = meshgrid(subunitCurrent(:,m)) - meshgrid(subunitCurrent(:,m))';

                    % positive driving forces values
                    posDrivingForce = drivingForce;
                    posDrivingForce(posDrivingForce<0) = 0;

                    % negative driving force values
                    negDrivingForce = drivingForce;
                    negDrivingForce(negDrivingForce>0) = 0;

                    %flow of current into bipolars due to gap junctions
                    gjAdd = sum((posDrivingForce.*exp(-sDistances./spaceConstant)).*gjGainIn, 2);

                    %flow of current out bipolars via gap junctions
                    gjSubtract = sum((negDrivingForce.*exp(-sDistances./spaceConstant)).*gjGainOut, 2);

                    subunitCurrent(:,m) = subunitCurrent(:,m) + gjSubtract + gjAdd;
                end
            end
            response(p,q,:,:) = squeeze(response(p,q,:,:)) + subunitCurrent;
        end
    end
end
response = response / numIter;

% Pipe it through the output nonlinearity.
nlResponse = nlfun(nlParams, response);
% nlResponse = hillFunc(response, 0.3774, 2.5164); %nlfun(nlParams,response);
we = zeros(size(nlResponse));
tmp = (W*ones(1,size(nlResponse,4)))/sum(W);
for p = 1 : length(ctM)
    for q = 1 : length(gainValues)
        we(p,q,:,:) = tmp;
    end
end
% Scale output by the RGC receptive field weights.
nlResponse = nlResponse .* we;