function [bcWeights, bcActivation] = getBipolarActivation(bpCenters,rgcRFDiam,bpRFDiam,numBars)

[x,y] = meshgrid(linspace(-rgcRFDiam/2,rgcRFDiam/2,rgcRFDiam/4));

% Make the RGC RF.
rgcRF = exp(-(sqrt(x.^2 + y.^2).^2/(rgcRFDiam/2)^2));

% Calculate the bar width.
barWidth = size(rgcRF,1)/numBars;

% Create all of the bipolar cell receptive fields.
bcWeights = zeros(1,size(bpCenters,1));
bcActivation = zeros(numBars, size(bpCenters,1));
for k = 1 : size(bpCenters,1)
    xy = bpCenters(k,:) - rgcRFDiam/2;
    bcRF = exp(-(sqrt((x+xy(1)).^2 + (y+xy(2)).^2).^2/(bpRFDiam/2)^2));
    
    % Get the pooling contribution of this bipolar to the RGC RF.
    bcWeights(k) = sum(sum(rgcRF .* bcRF));
    % Get the degree of activation by each bar.
    for m = 1 : numBars
        index = round((m-1)*barWidth + (1 : barWidth));
        bcActivation(m,k) = sum(sum(bcRF(index,:))) / sum(sum(bcRF));
    end
end
