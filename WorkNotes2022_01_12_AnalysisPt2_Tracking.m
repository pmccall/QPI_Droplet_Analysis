function WorkNotes2022_01_12_AnalysisPt2_Tracking

%%% Adapted from WorkNotes2021_12_21_Basusree_ControlsTracking.m %%%

%% 1.1 Load data
cd '/Users/mccall/Desktop/AnalysisCode_QPM/SampleData'
tmp = load('WorkNotes2022-01-12_SampleAnalysisPt1_Fitting.mat');
DS = tmp.DS;
FileName = tmp.FileName;
clear tmp

%% 1.2 Restructure data
nFile = length(FileName);
for f = 1:nFile
    nFOV = length(DS{f});
    for j = 1:nFOV
        nZ = length(DS{f}{j});
        tmp = [];
        for i = 1:nZ
            tmp = [tmp;DS{f}{j}{i}];
        end
        DS{f}{j} = tmp;
    end
end
clear tmp

%%% Conclusions %%%
% 1. Ran

%% 2. Track Sphere and Spherical Cap Fits

TrackParam = struct(...
    'mem',[],...        %%%% Set dynamically below %%%%
    'dim',2,...
    'good',2,...
    'quiet',0);
TrackMaxDisp = 12;      %%%%%%%%%% HEURISTIC %%%%%%%%%%

tic
nFile = length(FileName);
for a = 1:2
    Tracks = cell(nFile,1); % {dr off} x {'Sphere','SphericalCap'}
    for f = 1:nFile
        nFOV = length(DS{f});
        nZ = length(DS{f}{1});
        TrackParam.mem = nZ;
        Tracks{f} = [];
        for j = 1:nFOV
            ToDo = 1;
            % Uncomment section below to specify problematic FOVs to skip
            %{
            if f == 6
                if j == 6
                    ToDo = 0;   % Skip
                end
            elseif f == 8
                if j == 16
                    ToDo = 0;
                end
            elseif f == 9
                if j == 9
                    ToDo = 0;
                end
                if j == 16
                    ToDo = 0;
                end
            end
            %}
            if ToDo == 1
                Ind = DS{f}{j}(:,17) == a;       % select fit type
                tmp = DS{f}{j}(Ind,1:18);
                NewTracks = track(tmp,TrackMaxDisp,TrackParam);
                nTrack = length(NewTracks(:,1));
                NewBlock = [NewTracks,j*ones(nTrack,1)];
                Tracks{f} = [Tracks{f}; NewBlock];
            end
        end
    end
    switch a
        case 1
            TrackS = Tracks;
        case 2
            TrackSC = Tracks;
    end
end
toc
clear Tracks

%%% Conclusions %%%
% 1. Indicate FOVs or planes skipped. (none.)
% 2. Elapsed time is .... (0.722287 seconds.)

%% 3.1 Calculate MAR2PL for Sphere and Spherical cap fits
% The ith row in MAR2PLS(C) corresponds to the ith file in DS.
% Each cell in MAR2PLS(C) is an N-by-20 array, where N is the number of
%   tracked particles in the ith file in DS.
% The structure of each row of MAR2PLS(C){i,1} is
% col 1: row of object centroid in image (pixels)
% col 2: column of object centroid in image (pixels)
% col 3: refractive index difference, delta-n (refractive index units)
% col 4: uncertainty in delta-n (refractive index units)
% col 5: object radius, R (micrometers)
% col 6: uncertainty in R (micrometers)
% col 7: xPosition within fitting box (pixels)
% col 8: uncertainty in xPosition within fitting box, (pixels)
% col 9: yPosition within fitting box (pixels)
% col 10: uncertainty in yPosition within fitting box, (pixels)
% col 11: background phase, phi0 (radians)
% col 12: uncertainty in phi0 (radians)
% col 13: height of object equatorial plane above surface, Zeq (micrometers)
% col 14: uncertainty in Zeq (micrometers)
% col 15: Adj R^2 from surface fit
% col 16: RMSE (root-mean-squared error) from surface fit
% col 17: code for fit type (1 = sphere, 2 = spherical cap)
% col 18: Frame number for which AdjR2 is maximal for this object
% col 19: Particle ID
% col 20: Field of view (FOV) ID

Ind = [1];
ColIDArray = [15,20,19];
MAR2PLS = CalcDropletMaxAdjR2FromTracks(TrackS(Ind,1),ColIDArray);
MAR2PLSC = CalcDropletMaxAdjR2FromTracks(TrackSC(Ind,1),ColIDArray);
clear Ind

%%% Conclusions %%%
% 1. Ran

%% 3.2 Set Filter values
AdjR2_thresh = 0.90;    % Typical value is 0.98


%% 4. Plot each FOV in a different color 
nSample = length(MAR2PLSC(:,1));
SampleName = cell(nSample,1);
for i = 1:nSample
    SampleName{i} = FileName(i).Prefix;
end

%%% To specify subplot position for multiple samples, assign PlotInd and
%%% uncomment subplot. For single sample, comment out PlotInd and subplot.

%PlotInd = [1,2,3,6,7,8,4,9,5,10];
figure
for i = 1:nSample
    %subplot(2,5,PlotInd(i)); %%%%%%%%%%% HARDCODED %%%%%%%%%%%%
    tmp = MAR2PLSC{i};
    Ind = tmp(:,15) >= AdjR2_thresh;
    tmp = tmp(Ind,:);
    nFOV = max(tmp(:,20));
    cMap = ametrine(nFOV);
    Leg = cell(nFOV,1);
    for j = 1:nFOV
        Ind = tmp(:,20) == j;
        subset = tmp(Ind,:);
        x = subset(:,5);
        dx = subset(:,6);
        y = subset(:,3);
        dy = subset(:,4);
        errorbar(x,y,-dy,dy,-dx,dx,'o','color',cMap(j,:));
        hold on
        Leg{j} = ['n = ',num2str(length(y))];
    end
    xlabel('Radius, (\mum)');
    ylabel('\Deltan (-)');
    legend(Leg,'location','southeast');
    title(SampleName{i});
end

%%% Conclusions %%%
% 1. No obvious deviations from FOV to FOV.

%% save
cd '/Users/mccall/Desktop/AnalysisCode_QPM/SampleData'
saveas(gcf,'2022-01-12_SampleData_dnVsR_1x1plot.png')

%% 5. Violin plots comparing dn distributions for all samples
PlotInd = [1];  %%% Specify plotting order for multiple samples here.
tmpDS = cell(nSample,1);
for i = 1:nSample
    tmp = MAR2PLSC{i,1};
    tmp = FilterQPhaseParticleList(tmp);
    Ind = tmp(:,15) >= AdjR2_thresh; % Filter on AdjR2
    tmp = tmp(Ind,3);
    TargInd = PlotInd(i);
    tmpDS{TargInd} = tmp;
end

VRA = PrepareArrayForViolinPlot(tmpDS);

% Determine BWparam
ma = max(VRA,[],1);
mi = min(VRA,[],1);
dy = ma-mi;
BWparam = min(dy)*.15;

% Generate ViolinPlot
figure
v = violinplot(VRA,[],'BandWidth',BWparam);
ylabel('\Deltan (-)');
cMap = ametrine(nSample);   %%% Specify colormap here.
for i = 1:nSample
    v(i).ViolinColor = cMap(i,:);
end
a = gca;
a.XTickLabels = SampleName(PlotInd);
title(['Fit Adj. R^2 >= ',num2str(AdjR2_thresh)]);
tmp = ~isnan(VRA);
tmp = sum(tmp,1);
disp(tmp')

%%% Conclusions %%%
% 1. Ran as expected. 
% 2. Looks silly for just one sample.

%% save
cd '/Users/mccall/Desktop/AnalysisCode_QPM/SampleData'
saveas(gcf,'2022-01-12_SampleData_dnViolin.png')

%% 6. Write [dn,ddn,N] and MAR2PLSC as outputs
dn = zeros(nSample,3);  % [dn, ddn, N];
for i = 1:nSample
    tmp = MAR2PLSC{i};
    tmp = FilterQPhaseParticleList(tmp);
    Ind = tmp(:,15) >= AdjR2_thresh; % Weak filter on AdjR2
    tmp = tmp(Ind,3);
    dn(i,1) = mean(tmp);
    dn(i,2) = std(tmp);
    dn(i,3) = length(tmp);
end

%%% Conclusions %%%
% 1. Ran as expected.

%% save dn
save('WorkNotes2022-01-12_AnalysisPt2_Tracking.mat',"dn","MAR2PLSC");