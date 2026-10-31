function bpCenters = simulateBipolarArray(maxRad,bipolarSpacing,noiseSD,positionSeed)

Rad3Over2 = sqrt(3) / 2;
[X, Y] = meshgrid(1:1:2*maxRad);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

%set spacing
X = X * bipolarSpacing;
Y = Y * bipolarSpacing;

Ind = (X <= 2*maxRad & Y <= 2*maxRad);
%Ind = (X-bipolarSpacing_pix*3 <= stimSize & Y-bipolarSpacing_pix*3 <= stimSize);
X = X(Ind); %- bipolarSpacing_pix*1.5;
Y = Y(Ind); %- bipolarSpacing_pix*1.5;


%adding random noise
%noiseVar = 4.0; %microns

noiseToAdd = 2*sqrt(noiseSD);

rng(positionSeed, 'twister');
L = length(X);
X = X + noiseToAdd*randn(L,1);
Y = Y + noiseToAdd*randn(L,1);

%X,Y are now coordinates of bipolar centers

bpCenters = [X, Y]; %add 100 pixels to leave room for RFs