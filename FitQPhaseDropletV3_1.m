function [DS, Filtered, varargout] = FitQPhaseDropletV3_1(FileName, AnalysisOption, PlottingOption)

% Written by PMM on 2019-03-25

%% Description
% Based in part on WorkNotes2019_02_18_QPhaseImageProcessingPipeline.m

% Inputs:
% FileName: name of quantitative phase image file
% AnalysisOption is an optional input. It is a 1xm cell, where each entry 
%   contains a parameter related to the analysis. The code accepts cells 
%   with m = 4,5,6,7,8.
%   -AnalysisOption{1} = FilterDropSizeFromMask can be [] or a scalar. The
%       scalar determines the minimum size (measured as number of included
%       pixels) of a thresholded region in the binary mask that will be
%       kept. Regions with FilterDropSizeFromMask or fewer pixels will be
%       discarded from the analysis. If this parameter is empty (i.e. []),
%       then all regions will be included in the anlysis. A typical value
%       is 4.
%   -AnalysisOption{2} = FilterDropShapeFromMask can be [] or a scalar. The
%       scalar determines the minimum eccentricity of a thresholded region
%       in the binary mask that will be kept. Regions with eccentricies
%       greater than or equal to FilterDropShapeFromMask will be discarded
%       from the analysis. If this parameter is empty, then no regions are
%       discarded on the basis of shape.
%   -AnalysisOption{3} = FitWindowOpt is a string, and takes the values
%       "Fixed" or "SelfAdjusted". If Fixed, then a fitting window of fixed
%       size is used (number is hard-coded). If SelfAdjusted, then the 
%       fitting window size is determined on the basis of the bounding box
%       of the region.
%   -AnalysisOption{4} = FitWindowExtension can be [] or a scalar. The
%       scalar determines the number of pixels by which the bounding box
%       will be extended on either side when a SelfAdjusted fitting window
%       is used. A typical value is 3.
%   -AnalysisOption{5} = lNoise is a scalar >= 0. This is an optional 
%       parameter related to spatial filtering of input data. If 
%       lNoise = 0, then no filtering of the input image is performed. If 
%       non-zero, then a gaussian filter with characteristic lengthscale 
%       lNoise is applied to the images used to identify objects. Note that
%       the fitting is always performed on the unfiltered image.
%   -AnalysisOption{6} = nSig is a scalar > 0. This is an optional 
%       parameter used to determine the number of standard deviations above 
%       the modal value of the pixel intensity histogram at which the 
%       threshold used for  the binary mask should be set. A typical value 
%       is 3. 
%   -AnalysisOption{7} = ThreshAdjR2 is a scalar x, with 0 <= x < 1. This
%       is an optional paramter, used to determine the threshold value of 
%       AdjR2, obtained from the fit of the phase surface, used for 
%       filtering the output DS. If no value is specified by the user, then 
%       the default value of 0.90 is used.

% Outputs:
% DS: An n-by-12 array, where n is the number of identified droplets. The
%   column structure of DS is 
%   (rPeak,cPeak, 4 fit params, 4 fit uncertainty, AdjR2, RMSE)
% Filtered: A k-by-12 array with the same structure as DS, where k is the
%   number of droplets for which AdjR2 > ThreshAdjR2.
% Optional output: FitOut is an n-by-4 cell array of the form
%   {SurfaceFitObj,fGOF,fOut,{hWindow,vWindow,x,y}}, where the first three
%   cells are the outputs of the fit, and the fourth cell includes the
%   necessary information to specificy where in the image the fit is
%   performed.

%% Updates

% Updated to V3_1 by PMM on 2025-01-13 according to the following changes:
% - Added logic in §1 to use 'ResolutionUnit' when determinining pixel
%   size. This is important for data acquired with SophiQ v10 or later, as
%   those .ome.tiff images store XResolution and YResolution in cenimeters
%   instead of microns. If 'ResolutionUnit' = 'None', then units of microns
%   are assumed, consistent with past convention. This is to ensure
%   backwards-compatibility.

% - Note: the two updates below from 2024-10-11 and 2024-12-04 should also
% be considered part of the substantive update to v3_1.

% Updated by PMM on 2024-12-04 according to the following change:
% - Fixed bug that resulted in poor results for samples with contact angles
%   below 90 degrees (high contact angle case). The issue was that the
%   default lower-bound on zMid was 0! The code hadn't been allowed to look for
%   minima corresponding to lower contact angles. Now the mimimum is - the
%   estimate of R.

% Updated by PMM on 2024-10-11 according to the following change:
% - Added conditional statement in fit parameter initialization step that
% adjusts the expected radius when the initial guess for zMid is < 0
% (highly wetting = low contact angle).

% Updated by PMM on 2020-06-25 according to the following changes:
% - Default upper bound on R in fit functions is changed from x(end)/2 to
%       max(x(end),y(end)). This GREATLY improves fit quality for droplets 
%       not fully in the FOV.
% - Default lower bounds on xc and yc is changed from 1 MicronPerPix to 0.
% - Default upper bound on yc is changed from x(end) to y(end).
%   Particularly for droplets at the boundary, the sizes of x and y can
%   differ.
% - Added optional input to set Zmid0 as a fraction of the estimate of R. 
%   Default is 0.75.
% - Improved input determination logic so that AnalysisOptions 11,12,13 can
%   will assign default values if these parameters are passed as [].

% Updated by PMM on 2020-06-22 according to the following changes:
% - Added optional input TolFunTolX to specify the fitting tolerances for
%   each iteration. Default is [1E-6,1E-6];

% Updated by PMM on 2020-06-17 according to the following changes:
%   Things to change:
%{
    - Check type of AnalysisOption{8}: string or cell of strings
    -   If cell, determine nElem (number of fitting steps per drop)
    - Define SurfaceFitType and LengthStartPointVec as arrays
    - Fill in these arrays by looping over AnalysisOption{8} contents
    - 


    - Reduce repetition of FitFunction equations; represent each once
    - For fixed parameters in FitFunctions, add as ,num2str(value) so that the
        value can be passed explicitly.
    - DS output will now have a row for each fit iteration, and therefore
        needs an extra column indicating the iteration
    - Hold DS output dimensions constant, fill in unfit params with NaN
    - Fix output format at
    [x,y,dn,ddn,r,dr,xc,dxc,yc,dyc,zMid,dzMid,phi0,dPhi0,AdjR2,RMSE,Iter,ID,frame]
%}
%   Major update from V2.
% - Added ability to perform sequential fitting, where the fit parameters
%   of a fit to one model are passed to a second model. The fit parameters
%   from only the final round of fitting are returned.
% - Added LowContactAngle function

% Updated by PMM on 2020-06-08 according to the following changes:
%   Major update.
% - Added SphericalCap_tbPhi0_RlessZmidPos, Sphere_Phi0, Sphere_fixPhi0, SphericalCap_fixxPhi0_RlessZmidPos
% - Added AnalysisOption{10} = 1 -> fit only data > Phi0Est.

% Updated by PMM on 2020-05-03/04 according to the following changes:
%   Major update.
% - Added option to specify between two fit different surface fit
%   functions, that of a sphere (the original fit function), and that of a
%   spherical cap. The latter is a better approximation for the shape
%   adopted by a sessile liquid droplet on a solid substrate when the size
%   of the droplet is less than the capillary length given by
%   lc = sqrt(gamma/(drho*g)). The fit function is specified by the new
%   input AnalysisOption{8}, which can be either 'Sphere' or
%   'SphericalCap'. Default is 'Sphere'. 
% Note: 'SphericalCap' includes two additional fit parameters, Zmid and 
%   phi0, which are not present in the fit function 'Sphere'. As a result, 
%   the output DS contains 4 additional rows when 'SphericalCap' is used.
% - In code section 4, added conditional statements to interpret new input
%   and assign default ('Sphere').
% - In code section 4, added switch to specify fit function basd on
%   SurfaceFitType being 'Sphere' or 'SphericalCap'.
% - In code section 5, added switch to specify Startpoint, LowerBound, and
%   UpperBound for all fit parameters based on SurfaceFitType.
%       - phi0 start is phi0Est from code section 2, bounds are ±abs(th).
%       - Zmid start is x(end)/4, bounds are 0 and x(end).
% - In code section 5, update number of columns in DS initialization to
%   reflect the variable number of fit parameters (4 or 6) depending on
%   SurfaceFitType.
% - In code section 6, update DS column ID for AdjR2 filtering based.

% Updated by PMM on 2020-02-29 according to the following changes:
% - Added assignments to varargout in the conditional early termination
%   statements of sections 1 and 3.1
% - In code section 2, changed the criteria for dnSign to be +1 from 
%   FitParams(2) <= median(imFilt(:)); to skewness(imFilt(:)) >= 0. This
%   change helps avoid a bug encountered with very low noise and high dn,
%   particularly in simulated data.

% Updated by PMM on 2020-02-10 according to the following changes:
% - Added varargout to function call: If a third output is present, it
%   contains an n-by-4 cell array of the form
%   {SurfaceFitObj,fGOF,fOut,{hWindow,vWindow,x,y}}, where the first three
%   cells are the outputs of the fit, and the fourth cell includes the
%   necessary information to specificy where in the image the fit is
%   performed.
% Code section 5, added two conditional statements:
%   - At section beginning, initialized an nDrop-by-4 cell array if nargout
%       > 2
%   - Towards end of section, added line writing f, fGOF, fOut, and the 
%       image index information into the varargout cell array for each droplet.

% Updated by PMM on 2019-10-23 according to the following changes:
% Code section 2:
%   -Added conditional statement on FitParams(2) <= DistributionMedian,
%   which is used to specify whether dnSign = +1 or -1.
% Code section 3.1:
%   -Improved procedure so that the no assumptions on the sign or magnitude
%   of th are necessary (i.e. th can now be <0, =0, 0<th<1, 1, or >1).
%   -Added switch on dnSign to toggle whether pixels below or above th are
%   blacked out.
% Code section 5:
%   -Introduced dnFitBounds = [dnInitialGuess,dnUpperBound,dnLowerBound]
%   -Added switch in dnSign which specifies the contents (and sign) of the
%   entries in dnFitBounds.
%   -Increased the default max dn from 0.2 to 0.5.

% Updated by PMM on 2019-09-04 according to the following changes:
%  Inputs:
%   -Added additional (optional) parameter to AnalysisOption used to
%   indicate whether the input image should be gaussian filtered for
%   determining particle locations. AnalysisOption{5} = lNoise. No
%   filtering is done for lNoise = 0. For lNoise > 0, the a gaussian filter
%   with standard deviation lNoise is applied to the images used in
%   sections 2. and 3.1. To enable backwards-compatibility, no filtering is
%   performed if AnalysisOption only contains 4 entries.

%   -Added additional (optional) parameter to AnalysisOption used to
%   determine the number of standard deviations above the modal value of
%   the pixel intensity histogram at which the threshold used for the
%   binary mask should be set. AnalysisOption{6} = nSig. A typical value is 
%   3. To enable To enable backwards-compatibility, the default value of 3 
%   is used if AnalysisOption only contains fewer than 6 entries.

%   -Added additional (optional) parameter to AnalysisOption, used to
%   determine the threshold value of AdjR2 used for filtering the DS
%   output. AnalysisOption{7} = ThreshAdjR2, and the default value is 0.90.
%   To enable backwards-compatibility, the default value is used if
%   AnalysisOption contains fewer than 7 entries.

%  Code section 0:
%   - Set values (in a backwards-compatibile manner) for lNoise and nSig.
%  Code section 2:
%   - Added conditional statement to apply a gaussian filter with standard
%   deviation lNoise to the input image if the input lNoise > 0. Otherwise
%   the "filtered" image is equivalent to the input image.
%   - Pass the filtered image to histogram
%   - Use the parameter nSig to determine the threshold
%  Code section 3:
%   - Make black-and-white mask by thresholding the filtered image
%  Description:
%   - Added description of the AnalysisOption and PlottingOption inputs.
%   - Added description of Filtered output.

% Updated by PMM on 2019-06-02 according to the following changes:
%   Code section 2:
%   -Added conditional statement so that if histogram has less than 3 bins,
%   the code outputs empty arrays, displays a message, then terminates.
%   Code section 3.1:
%   -Added conditional statement so that if no objects are detected, the
%   code outputs empty arrays, displays a message, then terminates.

% Updated by PMM on 2019-04-10 according to the following changes:
%   Code section 5:
%   -Added code to the 'Fixed' window case in Section 5 to truncate the
%   window at the image boundaries. 
%   -Adjusted code to round the centroid coordinates to the nearest integer 
%   -Define meshgrid points after truncation (when necessary).
%   -These fixes combine to resolve the error:
%       Index in position 1 is invalid. Array indices must be positive integers or logical values.
%           Error in FitQPhaseDroplet (line 185)
%           tmp = im(vWindow,hWindow);
%   Code section 2:
%   -Added call to histcounts when plotting is not requested to fix bug
%   where a histogram is plotted on top of whatever plot is currently
%   active

%% 0. Interpret inputs
if nargin >= 2
    nAnalysisParam = length(AnalysisOption);
    FilterDropSizeFromMask = AnalysisOption{1};     % scalar or []
    FilterDropShapeFromMask = AnalysisOption{2};    % scalar or []
    FitWindowOpt = AnalysisOption{3};               % string
    FitWindowExtension = AnalysisOption{4};         % scalar or 0
    if nAnalysisParam >= 5
        lNoise = AnalysisOption{5};                 % scalar >= 0
    else
        lNoise = 0;
    end
    if nAnalysisParam >= 6                          
        nSig = AnalysisOption{6};                   % scalar > 0
    else
        nSig = 3;
    end
    if nAnalysisParam >= 7
        ThreshAdjR2 = AnalysisOption{7};            % scalar > 0, < 1
    else
        ThreshAdjR2 = 0.9;
    end
    if nAnalysisParam >= 8
        if ischar(AnalysisOption{8}) == 1
            InputFitType = {AnalysisOption{8}};     % convert to cell array
        elseif iscell(AnalysisOption{8}) == 1
            InputFitType = AnalysisOption{8};
        else
            disp('Input fit type format unexpected: Check AnalysisOption{8}')
            return
        end
        nFitIter = numel(AnalysisOption{8});
        SurfaceFitType = cell(nFitIter,1);
        LengthStartPointVec = zeros(nFitIter,1);
        IndUserSpecStartPoint = cell(nFitIter,1);
        for i = 1:nFitIter
            if strcmp('Sphere',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 4;
            elseif strcmp('Sphere_Phi0',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 5;
            elseif strcmp('Sphere_fixPhi0',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 5;
            elseif strcmp('SphericalCap',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 6;
            elseif strcmp('SphericalCap_tbPhi0',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 6;
            elseif strcmp('SphericalCap_tbPhi0_RlessZmidPos',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 6;
            elseif strcmp('SphericalCap_fixPhi0_RlessZmidPos',InputFitType{i}) == 1
                SurfaceFitType{i} = InputFitType{i};
                LengthStartPointVec(i) = 6;
            else
                disp('Unrecognized SurfaceFitType. Check AnalysisOption{8}')
                return
            end
            IndUserSpecStartPoint{i} = false(1,LengthStartPointVec(i));
        end
    else
        nFitIter = 1;
        SurfaceFitType{1} = 'Sphere';  % Default
    end
    if nAnalysisParam >= 9
        % Convert inputs to same form
        % Then check input size; disp error as needed
        % Then overwrite IndUserSpecStartPoint as specified
        if iscell(AnalysisOption{9}) == 0   % Assume numerical array
            [nRow,~] = size(AnalysisOption{9});
            UserSpecStartPoint = cell(nRow,1);
            for i = 1:nRow
                UserSpecStartPoint{i} = AnalysisOption{9}(i,:);
            end
        else
            UserSpecStartPoint = AnalysisOption{9};
        end
        nSpecStartPoint = length(UserSpecStartPoint);
        if nSpecStartPoint > nFitIter
            disp('More start points given than fit iterations.')
            disp('Check size of AnalysisOption elements 8 and 9')
        else
            for i = 1:nSpecStartPoint
                if length(UserSpecStartPoint{i}) == LengthStartPointVec(i)
                    IndUserSpecStartPoint{i} = ~isnan(UserSpecStartPoint{i});
                else
                    disp('User-supplied StartPointVector wrong size for selected SurfaceFitType.')
                    disp('Check AnalysisOption elements 8 and 9')
                    return
                end
            end
        end
    end
    if nAnalysisParam >= 10
        RestrictFittingDomain = AnalysisOption{10};
    else
        RestrictFittingDomain = 0;  % Default, no restriction
    end
    if nAnalysisParam >= 11
        if ~isempty(AnalysisOption{11})
            ZmidPenalty = AnalysisOption{11};
        else
            ZmidPenalty = 10;   % Default, Applies to RlessZmidPos fits
        end
    end
    if nAnalysisParam >= 12
        if ~isempty(AnalysisOption{12})
            TolFunTolX = AnalysisOption{12};
        else
            TolFunTolX = ones(nFitIter,1)*[1E-6,1E-6];
        end
    end
    if nAnalysisParam >= 13
        if ~isempty(AnalysisOption{13})
            Zmid0_Default = AnalysisOption{13};
        else
            Zmid0_Default = 0.75;
        end
    end
else
    FilterDropSizeFromMask = [];
    FilterDropShapeFromMask = [];
    FitWindowOpt = 'Fixed';
    FitWindowExtension = 0;
    lNoise = 0;
    nSig = 3;
    ThreshAdjR2 = 0.9;
    SurfaceFitType{1} = 'Sphere';  % Default
    nFitIter = 1;
    RestrictFittingDomain = 0;  % Default, no restriction
    ZmidPenalty = 10;           % Default, Applies to RlessZmidPos fits
    TolFunTolX = ones(nFitIter,1)*[1E-6,1E-6];
    Zmid0_Default = 0.75;
end

if nargin >= 3
    FigOpt1 = PlottingOption(1);
    FigOpt2 = PlottingOption(2);
else
    FigOpt1 = 1;
    FigOpt2 = 1;
end

%% 1. Load image
if ischar(FileName) == 1
    % Assume file name
    imf = imfinfo(FileName);
    ResUnit = imf(1).ResolutionUnit;
    switch ResUnit
        case 'None' % Assume units of microns. This was the default prior to SophiQ v10 (i.e. for data acquired by PMM on or before 2023.12.04)
            PixPerMicron = imf(1).XResolution;
        case 'Centimeter' % Convert from centimeters to microns
            PixPerMicron = imf(1).XResolution/1000;
        case 'Micrometer' % This case doesn't appear yet at time of writing, but I anticipate it may be relevent in the future
            PixPerMicron = imf(1).XResolution;
    end
    MicronPerPix = 1/PixPerMicron;
    im = double(imread(FileName));
else
    % Assume array
    im = FileName;
    PixPerMicron = 6.3767;  % Default for 40x Objective on Qphase. Valid for G1 system.
    MicronPerPix = 1/PixPerMicron;
    imf(1).Width = length(im(1,:));
    imf(1).Height = length(im(:,1));
end

if FigOpt1 == 1
    figure
        subplot(2,2,1)
        imagesc(im)
        hold on
end

%% 2. Determine threshold
if lNoise > 0
    imFilt = imgaussfilt(im,lNoise);
else
    imFilt = im;
end

if FigOpt1 == 1
    subplot(2,2,2)
    Histo = histogram(imFilt);
    BinMid = Histo.BinEdges(1:(end-1)) + Histo.BinWidth/2;
    HistoValues = Histo.Values;
else
    [HistoValues,BinEdges] = histcounts(imFilt);
    BinWidth = BinEdges(2)-BinEdges(1);
    BinMid = BinEdges(1:(end-1)) + BinWidth/2;
end

if length(BinMid) < 3  % 3 points is minimum for fitting the histogram with a gaussian
    DS = [];
    Filtered = [];
    if nargout > 2
        varargout{1} = [];
    end
    disp([FileName,' is likely a blank image'])
    return
end

f = fit(BinMid',HistoValues','Gauss1');

if FigOpt1 == 1
    hold on
    set(gca,'YScale','log');
    plot(f)
    ylim([.5,1E5]);
    xlabel('Phase Delay (Radians)');
    ylabel('Count');
end

FitParams = coeffvalues(f); % [Amplitude, x_center, sqrt(2)*sigma]
skew = skewness(imFilt(:));
if skew >= 0
    dnSign = 1; % Phase objects are bright in image
else
    dnSign = -1;    % Phase objects are dark in image
end
th = FitParams(2) + dnSign*nSig*FitParams(3);   % ++++++++++++ HEURISTIC ++++++++++++
phi0Est = FitParams(2);
dPhi0 = FitParams(3)/sqrt(2);

%% 3.1 Identify droplets
bw = imFilt;

% Convert distribution to two delta functions positioned outside domain
BiggestValue = max(imFilt(:));
LowestValue = min(imFilt(:));
bw(bw>th) = 2*BiggestValue;
bw(bw<th) = -2*abs(LowestValue);
bw(bw==th) = 2*BiggestValue;
switch dnSign
    case 1  % Assign above threshold to 1
        bw(bw==2*BiggestValue) = 1;
        bw(bw==-2*abs(LowestValue)) = 0; 
    case -1 % Assign below threshold to 1
        bw(bw==2*BiggestValue) = 0;
        bw(bw==-2*abs(LowestValue)) = 1;              
end

cc = bwconncomp(bw);
reg = regionprops(cc,'area','centroid','boundingbox','eccentricity');
nDropTot = length(reg);

if isequal(nDropTot,0)
    DS = [];
    Filtered = [];
    if nargout > 2
        varargout{1} = [];
    end
    disp(['No objects detected in ',FileName])
    return
end

DropList = zeros(nDropTot,8);
for i = 1:nDropTot
    DropList(i,1) = reg(i).Area;
    DropList(i,2:3) = reg(i).Centroid;
    DropList(i,4:7) = reg(i).BoundingBox;
    DropList(i,8) = reg(i).Eccentricity;
end
if FigOpt1 == 1
    subplot(2,2,3)
        imagesc(bw)    
    subplot(2,2,4)
        histogram(DropList(:,1))
        xlabel('Region Area (pixels)')
        ylabel('Count');        
    subplot(2,2,1)
        plot(DropList(:,2),DropList(:,3),'oy');
end
%% 3.2 Filter droplets by size
if isempty(FilterDropSizeFromMask) == 0
    Area = DropList(:,1);
    Ind = Area > FilterDropSizeFromMask;
    if FigOpt1 == 1
        aInd = logical(1 - Ind);
        subplot(2,2,1)
        plot(DropList(aInd,2),DropList(aInd,3),'xr');
    end
    DropList = DropList(Ind,:);
    nDrop = length(DropList(:,1));
else
    nDrop = nDropTot;
end

%% 3.3 Filter droplets by shape (eccentricity)
if isempty(FilterDropShapeFromMask) == 0
    Eccentricity = DropList(:,8);
    Ind = Eccentricity < FilterDropShapeFromMask;
    if FigOpt1 == 1
        aInd = logical(1 - Ind);
        subplot(2,2,1)
        plot(DropList(aInd,2),DropList(aInd,3),'xg');
    end
    DropList = DropList(Ind,:);
    nDrop = length(DropList(:,1));
end

if FigOpt1 == 1
    subplot(2,2,1)
    plot(DropList(:,2),DropList(:,3),'ob');
end
%% 4. Set up fitting
CheckParam = {'phi0','Zmid'};
IndexPos = [11,13];

FitCoefficientName = cell(nFitIter,1);
nFitCoeff = zeros(nFitIter,1);
surfaceFitFunction = cell(nFitIter,1);
ColIndDS = cell(nFitIter,1);
for i = 1:nFitIter
    switch SurfaceFitType{i}
        case 'Sphere'
            FitCoefficientName = {'dn','R','xc','yc'};            
            FitExpression = ...
                '(4*pi/0.65)*dn*sqrt(R^2 - (x-xc)^2 - (y-yc)^2)*heaviside(R^2 - (x-xc)^2 - (y-yc)^2)';
        case 'Sphere_Phi0'
            FitCoefficientName = {'dn','R','xc','yc','phi0'};
            FitExpression = ...
                '(4*pi/0.65)*dn*sqrt(R^2 - (x-xc)^2 - (y-yc)^2)*heaviside(R^2 - (x-xc)^2 - (y-yc)^2) + phi0';
            Phi0Bounds = 'Default';
        case 'Sphere_fixPhi0'
            FitCoefficientName = {'dn','R','xc','yc','phi0'};
            FitExpression = ...
                '(4*pi/0.65)*dn*sqrt(R^2 - (x-xc)^2 - (y-yc)^2)*heaviside(R^2 - (x-xc)^2 - (y-yc)^2) + phi0';
            Phi0Bounds = 'Fixed';
        case 'SphericalCap'
            FitCoefficientName = {'dn','R','xc','yc','phi0','Zmid'};
            FitExpression = ['(2*pi/0.65)*dn',...
                '.*( sqrt(R.^2-(x-xc).^2-(y-yc).^2)',...
                '.*(1+heaviside(Zmid.^2+(x-xc).^2+(y-yc).^2-R.^2))',...
                '.*heaviside(R.^2-(x-xc).^2-(y-yc).^2)',...
                '+ Zmid.*heaviside(R.^2-Zmid.^2-(x-xc).^2-(y-yc).^2) )',...
                '+ phi0'];
            Phi0Bounds = 'Default';
        case 'SphericalCap_tbPhi0'
            FitCoefficientName = {'dn','R','xc','yc','phi0','Zmid'};
            FitExpression = ['(2*pi/0.65)*dn',...
                '.*( sqrt(R.^2-(x-xc).^2-(y-yc).^2)',...
                '.*(1+heaviside(Zmid.^2+(x-xc).^2+(y-yc).^2-R.^2))',...
                '.*heaviside(R.^2-(x-xc).^2-(y-yc).^2)',...
                '+ Zmid.*heaviside(R.^2-Zmid.^2-(x-xc).^2-(y-yc).^2) )',...
                '+ phi0'];
            Phi0Bounds = 'TightBounds';
        case 'SphericalCap_tbPhi0_RlessZmidPos'
            FitCoefficientName = {'dn','R','xc','yc','phi0','Zmid'};
            FitExpression = ['(2*pi/0.65)*dn',...
                '.*( sqrt(R.^2-(x-xc).^2-(y-yc).^2)',...
                '.*(1+heaviside(Zmid.^2+(x-xc).^2+(y-yc).^2-R.^2))',...
                '.*heaviside(R.^2-(x-xc).^2-(y-yc).^2)',...
                '+ Zmid.*heaviside(R.^2-Zmid.^2-(x-xc).^2-(y-yc).^2) )',...
                '+ phi0',...
                '+ ',num2str(ZmidPenalty),'*((Zmid-R).^2)*heaviside(Zmid-R)'];
            Phi0Bounds = 'TightBounds';
        case 'SphericalCap_fixPhi0_RlessZmidPos'
            FitCoefficientName = {'dn','R','xc','yc','phi0','Zmid'};
            FitExpression = ['(2*pi/0.65)*dn',...
                '.*( sqrt(R.^2-(x-xc).^2-(y-yc).^2)',...
                '.*(1+heaviside(Zmid.^2+(x-xc).^2+(y-yc).^2-R.^2))',...
                '.*heaviside(R.^2-(x-xc).^2-(y-yc).^2)',...
                '+ Zmid.*heaviside(R.^2-Zmid.^2-(x-xc).^2-(y-yc).^2) )',...
                '+ phi0',...
                '+ ',num2str(ZmidPenalty),'*((Zmid-R).^2)*heaviside(Zmid-R)'];
            Phi0Bounds = 'Fixed';
    end
    nFitCoeff(i) = length(FitCoefficientName);
    surfaceFitFunction{i} = fittype(FitExpression,...
        'independent',{'x','y'},...
        'dependent','z',...
        'coefficients',FitCoefficientName);
    ColIndDS{i} = [3,5,7,9];
    for j = 1:numel(CheckParam)
        if find(contains(FitCoefficientName,CheckParam{j}))
            ColIndDS{i} = [ColIndDS{i},IndexPos(j)];
        end
    end
end
[nFitCoeffMax,~] = max(nFitCoeff);
%[nFitCoeffMax,IndFitIterMaxParam] = max(nFitCoeff);
FieldNameGOF = {'adjrsquare','rmse'};
nFieldNameGOF = length(FieldNameGOF);
colAdjR2 = 2 + 2*nFitCoeffMax + 1;       % Assumes FieldNamesGOF{1} = 'adjrsquare'
DSnCol = 2 + 2*nFitCoeffMax + nFieldNameGOF + 1;
%{
ParamColInd = 2 + (1:2:(2*nFitCoeffMax-1));
ParamList = coeffnames(surfaceFitFunction{IndFitIterMaxParam});
ColID = cell(2,nFitCoeffMax);
for i = 1:nFitCoeffMax
    ColID{1,i} = ParamList{i};
    ColID{2,i} = ParamColInd(i);
end
%}
%fOptions = fitoptions(surfaceFitFunction)
%% 5. Loop over all peaks and fit
DS = nan(nDrop*nFitIter,DSnCol);
% Columns in DS (rPeak,cPeak, 6 fit params, 6 fit uncertainty, AdjR2, RMSE, FitIteration) = 17 
if nargout > 2
    varargout{1} = cell(nDrop*nFitIter,4);   % {f,fGOF,fOUT,{hWindow,vWindow,x,y}}
end
for i = 1:nDrop
    %%%%%%%%%% Determine FitWindow %%%%%%%%%%
    switch FitWindowOpt
        case 'SelfAdjusted' % Use boundingbox with Optional Extension
            LHS = ceil(DropList(i,4)) - FitWindowExtension;
            RHS = ceil(DropList(i,4)+DropList(i,6)) + FitWindowExtension;
            TopSide = ceil(DropList(i,5)) - FitWindowExtension;
            BottomSide = ceil(DropList(i,5)+DropList(i,7)) + FitWindowExtension;
            if LHS < 1
                LHS = 1;
            end
            if TopSide < 1
                TopSide = 1;
            end
            if RHS > imf(1).Width
                RHS = imf(1).Width;
            end
            if BottomSide > imf(1).Height
                BottomSide = imf(1).Height;
            end
            hWindow = (LHS:RHS);
            vWindow = (TopSide:BottomSide);
            x = MicronPerPix*(1:length(hWindow))';
            y = MicronPerPix*(1:length(vWindow))';
        case 'Fixed'
            TypicalDropletSize = 6;%13;    % ++++++++++++ HEURISTIC ++++++++++++
            FitWindowSize = TypicalDropletSize*2+1; % ++++++++++++ HEURISTIC ++++++++++++
            rFitWindow = (round(FitWindowSize-1))/2;
            xc = round(reg(i).Centroid(2));
            yc = round(reg(i).Centroid(1));
            vWindow = (yc-rFitWindow):(yc+rFitWindow);
            hWindow = (xc-rFitWindow):(xc+rFitWindow);
            vWindow = vWindow(vWindow >= 1);
            vWindow = vWindow(vWindow <= imf(1).Height);
            hWindow = hWindow(hWindow >= 1);
            hWindow = hWindow(hWindow <= imf(1).Width);
            x = MicronPerPix*(1:length(hWindow))';
            y = MicronPerPix*(1:length(vWindow))';
    end
    tmp = im(vWindow,hWindow);
    [X,Y,Z] = prepareSurfaceData(x,y,tmp);
    if RestrictFittingDomain == 1
        switch dnSign
            case 1
                IndFitDomain = Z >= phi0Est;
            case -1
                IndFitDomain = Z <= phi0Est;
        end
        X = X(IndFitDomain);
        Y = Y(IndFitDomain);
        Z = Z(IndFitDomain);
        %disp('restricted domain for fitting')          %%%% DISPLAY %%%%
    end
    %%%%%%%%%% Initialization %%%%%%%%%%
    switch dnSign
        case 1
            dnFitBounds = [.1,.5,.00001];
        case -1
            dnFitBounds = [-0.1,-0.00001,-0.5];
    end
    CheckParamDefaultStartPoint = [phi0Est,x(end)/4];   % Default initial guesses
    nCheckParam = numel(CheckParamDefaultStartPoint);
    for k = 1:nFitIter
        RowInd = nFitIter*(i-1)+k;
        DS(RowInd,1:2) = DropList(i,2:3);
        DS(RowInd,end) = k;
        StartPoint = nan(1,nFitCoeff(k));
        UB = nan(1,nFitCoeff(k));
        LB = nan(1,nFitCoeff(k));
        UB(1:4) = [dnFitBounds(2),max(x(end),y(end)),x(end),y(end)];
        LB(1:4) = [dnFitBounds(3),MicronPerPix,0,0];
        if k == 1
            StartPoint(1:4) = [dnFitBounds(1),5,x(end)/2,x(end)/2];     % Default initial guesses
        else
            ParamFromLastFit = DS(RowInd-1,ColIndDS{k-1});
            StartPoint(1:4) = ParamFromLastFit(1:4);
            for a = 5:(4+nCheckParam)
                if nFitCoeff(k) >= a
                    if numel(ParamFromLastFit) >= a
                        StartPoint(a) = ParamFromLastFit(a);
                    else
                        if a == 6   % If Zmid not fit previously, estimate from fit of R param
                            StartPoint(a) = Zmid0_Default*StartPoint(2);     %%%% User-defined as of 2020.06.25 %%%%
                            if Zmid0_Default < 0
                                StartPoint(2) = sqrt(ParamFromLastFit(2)^2 + StartPoint(6)^2);      %%% Added 2024.10.11 to handle low contact angles
                            end
                        end
                    end
                end
            end
        end
        for a = 5:(4+nCheckParam)
            if nFitCoeff(k) >= a
                if isnan(StartPoint(a))     % Initialize needed parameters if not done already above
                    StartPoint(a) = CheckParamDefaultStartPoint(a-4);
                end
                switch a
                    case 5                  % Phi0
                        switch Phi0Bounds
                            case 'Default'
                                UB(a) = abs(th);
                                LB(a) = -abs(th);
                            case 'TightBounds'
                                UB(a) = phi0Est+dPhi0;
                                LB(a) = phi0Est-dPhi0;
                            case 'Fixed'
                                UB(a) = phi0Est;
                                LB(a) = phi0Est;
                        end
                    case 6                  % Zmid
                        UB(a) = x(end);
                        %LB(a) = 0;                 %%% Previous value, prior to 2024.12.04
                        LB(a) = -StartPoint(2);     %%% Added 2024.12.04 to allow for low contact angles 
                end
            end
        end
        if exist('IndUserSpecStartPoint','var')
            StartPoint(IndUserSpecStartPoint{k}) = UserSpecStartPoint{k}(IndUserSpecStartPoint{k});
        end
        %disp(StartPoint)                            %%%% DISPLAY %%%%
        %%%%%%%%%% Fit %%%%%%%%%%
        [f, fGOF, fOutput] = fit([X,Y],Z,surfaceFitFunction{k},...
            'StartPoint',StartPoint,...
            'Upper',UB,...
            'Lower',LB,...
            'TolFun',TolFunTolX(k,1),...
            'TolX',TolFunTolX(k,2));
        %disp(fOutput)                              %%%% DISPLAY %%%%
        CoeffValue = coeffvalues(f);
        ConfInt = confint(f);
        Uncertainty = CoeffValue - ConfInt(1,:);
        for j = 1:nFitCoeff(k)
            ColInd = ColIndDS{k}(j);
            DS(RowInd,ColInd) = CoeffValue(j);
            DS(RowInd,ColInd+1) = Uncertainty(j);
        end
        for j = 1:nFieldNameGOF
            ColInd = 2+2*nFitCoeffMax+j;
            DS(RowInd,ColInd) = fGOF.(FieldNameGOF{j});
        end
        if nargout > 2
            varargout{1}{RowInd,1} = f;
            varargout{1}{RowInd,2} = fGOF;
            varargout{1}{RowInd,3} = fOutput;
            varargout{1}{RowInd,4} = {hWindow,vWindow,x,y};
        end
    end
    clear X Y Z f fGOF fOutput CoefValue ConfInt Uncertainty
    if rem(i,100) == 1
        disp(['Processing ',num2str(i),' of ',num2str(nDrop),' Droplets'])
    end
end

%% 6. Post-processing and filtering
Ind = DS(:,colAdjR2) > ThreshAdjR2;   % Include only dropf for which Adj-R^2 > 0.9
Filtered = DS(Ind,:);

if FigOpt2 == 1
    figure
    for j = 1:2
        switch j
            case 1
                DSplot = DS;
            case 2
                DSplot = Filtered;
        end
        subplot(2,4,(j-1)*4+1)
        histogram(DSplot(:,3))
        xlabel('\Deltan (-)');
        ylabel('counts');

        subplot(2,4,(j-1)*4+2)
        histogram(DSplot(:,5));
        xlabel('R, \mum');
        ylabel('counts');

        dn = DSplot(:,3);
        R = DSplot(:,5);
        dnErr = DSplot(:,4);
        dRerr = DSplot(:,6);
        AdjR2 = DSplot(:,colAdjR2);

        subplot(2,4,(j-1)*4+3)
        errorbar(R,dn,-dnErr,+dnErr,-dRerr,+dRerr,'ok')
        xlabel('R, \mum');
        ylabel('\Deltan (-)');

        subplot(2,4,(j-1)*4+4)
        errorbar(AdjR2,dn,-dnErr,+dnErr,'ok')
        xlabel('Adj. R^2 (-)');
        ylabel('\Deltan (-)');
        title(['N = ',num2str(length(DSplot(:,1)))]);
    end
end