function [Kept, varargout] = FilterQPhaseParticleList(ParticleList,OptionalParam)

% Written by PMM on 2020.07.02

%%% PURPOSE %%%
% To remove objects from a (tracked) particle list due to one or more of 
%   the following issues:
%       1. Not a droplet (e.g. dead pixel)
%       2. Fit parameter uncertainty not available
%       3. Fit (quality) is bad because
%           3.1 Object center is outside of FOV
%           3.2 Object morphology deviates too much from circular
%                   cross-section (e.g. contact pinning, or additional
%                   droplets present in fitting window)
%           3.3 Object is moving significantly during acquisition (e.g.
%                   small and not in complete contact with coverglass)



%%% INPUTS %%%
% ParticleList is the output of e.g. FitQPhaseDropletV3 before or after
%   tracking.
% OptionalParam is a struct with the following fields

%%% OUTPUTS %%%
% FilteredList is an array with the same number of columns as ParticleList
%   containing all those rows of ParticleList that meet the filtering
%   criteria to be kept.
% If a second output argument is specified, then an array is returned 
%   containing al those particles that do NOT meet the criteria for being
%   kept.

%%% UPDATES %%%
% Updated on 2021-06-11 by PMM according to the following change
%   - Moved SkipSolution definition to outside the conditional statement on
%   nargin in Section 4. This fixed the bug whereby SkipSolution was not
%   defined if OP was passed but didn't include the SkipSolution field.
% Updated on 2020-10-21 by PMM according to the following changes
%   - Added OP field 'SkipSolution', because this code section seems buggy.
%   Set = 1 to skip, = 0 to perform filtering of droplets in solution.

%% 1. Not a droplet
% Note: these dead pixels are specific to the loaner QPhase system at the
%   MPI-CBG LMF. If the camera is replaced or the analysis is to be 
%   performed on images acquired on a different system, then the specific
%   coordinates need to be updated.

% ROIs of dead pixels and other fixed imaging issues.
xLB = [177.8,101];
xUB = [179.8,124];
yLB = [437.4,45];
yUB = [439.4,58];
nROI = length(xLB);

tmp = ParticleList;
x0 = tmp(:,1);
y0 = tmp(:,2);

%IndRemove = cell(3,1);
IndRemove{1} = false(length(x0),1);
for i = 1:nROI
    Ind1 = x0 >= xLB(i);
    Ind2 = x0 <= xUB(i);
    Ind3 = y0 >= yLB(i);
    Ind4 = y0 <= yUB(i);
    xInd = and(Ind1,Ind2);
    yInd = and(Ind3,Ind4);
    CurrentInd = and(xInd,yInd);
    IndRemove{1} = or(IndRemove{1},CurrentInd);   % Include newly identified objects
end
IndKeep{1} = ~IndRemove{1};
Kept = tmp(IndKeep{1},:);

RemovedTmp = tmp(IndRemove{1},:);
Removed{1} = [RemovedTmp,1*ones(sum(IndRemove{1}),1)];

%% 2. NaN in all fit parameters
tmp = Kept;

nFitParam = 6;
FitParamErrorColID = [4,6,8,10,12,14];
for i = 1:nFitParam
    CurrentInd = isnan(tmp(:,FitParamErrorColID(i)));
    if i == 1
        IndRemove{2} = CurrentInd;
    else
        IndRemove{2} = and(IndRemove{2},CurrentInd);
    end   
end

IndKeep{2} = ~IndRemove{2};
Kept = tmp(IndKeep{2},:);

RemovedTmp = tmp(IndRemove{2},:);
Removed{2} = [RemovedTmp,2*ones(sum(IndRemove{2}),1)];

%% 3. Droplet centroid is outside FOV
% NaN in xc and yc uncertainty
tmp = Kept;

xcOut = isnan(tmp(:,8));
ycOut = isnan(tmp(:,10));
IndRemove{3} = or(xcOut,ycOut);

IndKeep{3} = ~IndRemove{3};
Kept = tmp(IndKeep{3},:);

RemovedTmp = tmp(IndRemove{3},:);
Removed{3} = [RemovedTmp,3*ones(sum(IndRemove{3}),1)];

%% 4. Droplets in solution (and not in contact with the coverglass)

tmp = Kept;
dR_thresh = 1.0;   % um     Default value
dZ_zstack = 0.2;   % um     Default value
SkipSolution = 1;  % default is to skip

if nargin > 1
    if isfield(OptionalParam,'dR_thresh')
        dR_thresh = OptionalParam.dR_thresh;
    end
    if isfield(OptionalParam,'dZ_zstack')
        dZ_zstack = OptionalParam.dZ_zstack;
    end 
    if isfield(OptionalParam,'SkipSolution')
        SkipSolution = OptionalParam.SkipSolution;
    end
end
ColIndZmax = 13;
ColIndR = 5;
tmp = sortrows(tmp,ColIndR);
tmp = sortrows(tmp,ColIndZmax);
nParticle = length(tmp(:,1));
IndRemove{4} = false(nParticle,1);
if SkipSolution == 0    %   i.e. if false
    if nParticle < 2
        disp('Number of Kept particles after secton 3 is too small to filter for droplets in solution')
        return
    end
    for i = 2:nParticle
        if i == 2
            IndRef = 1;
        end
        dZmax = tmp(i,ColIndZmax) - tmp(IndRef,ColIndZmax);
        dR = tmp(i,ColIndR) - tmp(IndRef,ColIndR);
        Cutoff = -dR_thresh + dZ_zstack*dZmax;
        if dR >= Cutoff
            % Keep particle i, update IndRef
            IndRef = i;
        else
            % Remove particle i, keep IndRef
            IndRemove{4}(i) = true(1);
        end
    end
end
IndKeep{4} = ~IndRemove{4};
Kept = tmp(IndKeep{4},:);

RemovedTmp = tmp(IndRemove{4},:);
Removed{4} = [RemovedTmp,4*ones(sum(IndRemove{4}),1)];

%% Display output summary
nTot = length(ParticleList(:,1));
nKeep = length(Kept(:,1));
nRemove = sum(vertcat(IndRemove{1},IndRemove{2},IndRemove{3},IndRemove{4}),'all');
disp(['Remove ',num2str(nRemove),' of ',num2str(nTot),' objects; Keep ',num2str(nKeep),' of ',num2str(nTot)])

if nargin > 1
    varargout{1} = Removed;
end
