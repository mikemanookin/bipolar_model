function bpCurrents = getCoupling(bpCurrents,bpDistances,spaceConstant,gjGainIn,gjGainOut,modelfun,nlParams)
%
% INPUTS
% bpCurrents
% bpDistances
% spaceConstant
% gjGainIn
% gjGainOut

% pair-wise bipolar voltage differences
drivingForce = meshgrid(bpCurrents) - meshgrid(bpCurrents)';

% drivingForce = drivingForce / (length(bpCurrents)-1);

% positive driving forces values
posDrivingForce = drivingForce;
posDrivingForce(posDrivingForce<0) = 0;

% negative driving force values
negDrivingForce = drivingForce;
negDrivingForce(negDrivingForce>0) = 0;

%flow of current into bipolars due to gap junctions
gjAdd = sum((posDrivingForce.*exp(-bpDistances./spaceConstant)).*gjGainIn, 2);
% gjAdd = sum((posDrivingForce.*exp(-bpDistances./spaceConstant)).*gjGainIn, 2) / (length(bpCurrents)-1);

%flow of current out bipolars via gap junctions
gjSubtract = sum((negDrivingForce.*exp(-bpDistances./spaceConstant)).*gjGainOut, 2);
% gjSubtract = sum((negDrivingForce.*exp(-bpDistances./spaceConstant)).*gjGainOut, 2) / (length(bpCurrents)-1);

% bpCurrents = gjSubtract + gjAdd;
bpCurrents = bpCurrents + gjSubtract + gjAdd;
% 
% bpCurrents = modelfun(nlParams, bpCurrents);