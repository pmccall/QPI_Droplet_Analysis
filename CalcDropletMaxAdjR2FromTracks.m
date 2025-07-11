function MAR2PL = CalcDropletMaxAdjR2FromTracks(Tracks,ColIDArray)

% Written by PMM on 2020.07.17

%%% PURPOSE %%%
% For use in analysis of tracked Qphase z-stack data. This code determines,
%   for each droplet, the z-plane in which the phase surface fit to that
%   droplet was best, and returns a particle list containing one row per
%   tracked object from the oplane of best fit.

%%% INPUTS %%%
% Tracks is a cell array, each element of which contains a sorted list of
%   tracked particles returned by track.m or similar. Different cell
%   elements correspond to different files (e.g. different temperatures or
%   samples)
% ColIDArray is a 1x3 array indicating the column indicies of Tracks that
%   contain the following values: [AdjR2, FOV ID, Track ID]. The default is
%   [15,20,19].

%%% OUTPUTS %%%
% MAR2PL is a cell array with the same dimensions as Tracks. Each element
%   of the cell array contains a subset of the corresponding element of
%   Tracks, now with only one row per particle trajectory. The retained row
%   is that with the maximal value of ColIDArray(1) observed over the
%   trajectory.

%% 1. Interpret inputs
if nargin < 2
    ColID_AdjR2 = 15;
    ColID_FOV = 20;
    ColID_TrackID = 19;
else
    ColID_AdjR2 = ColIDArray(1);
    ColID_FOV = ColIDArray(2);
    ColID_TrackID = ColIDArray(3);
end

[nFile,nFitType] = size(Tracks);

%% 2. Determine MAR2PL

MAR2PL = cell(nFile,2);
for f = 1:nFile
    for i = 1:nFitType
        MaxAdjR2ParticleList = [];
        nFOV = max(Tracks{f,i}(:,ColID_FOV));
        for j = 1:nFOV
            IndCurrentFOV = Tracks{f,i}(:,ColID_FOV) == j;
            CurrentFOV = Tracks{f,i}(IndCurrentFOV,:);
            [~,First,~] = unique(CurrentFOV(:,ColID_TrackID),'first');
            [~,Last,~] = unique(CurrentFOV(:,ColID_TrackID),'last');
            nTrack = length(First);
            
            MaxAdjR2 = zeros(nTrack,1);
            for k = 1:nTrack    
                [~, RelIndMaxAdjR2] = max(CurrentFOV(First(k):Last(k),ColID_AdjR2));
                MaxAdjR2(k) = RelIndMaxAdjR2 + First(k) - 1;
            end
            NewBlock = CurrentFOV(MaxAdjR2,:);
            MaxAdjR2ParticleList = vertcat(MaxAdjR2ParticleList,NewBlock);
        end
        MAR2PL{f,i} = MaxAdjR2ParticleList;
    end
end