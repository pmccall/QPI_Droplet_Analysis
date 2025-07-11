# QPI_Droplet_Analysis
Extract the shape and refractive index difference of droplets from quantitative phase images (QPI)

% Listing of the code dependencies for the analysis of Quantitative Phase
%   Imaging (QPI) images of biomolecular condensates. Code written by
%   Patrick M McCall except 
%       track.m (Daniel Blair & Eric Dufresne https://site.physics.georgetown.edu/matlab/)
%       ametrine.m (Geissbuehler)
%       violinplot.m (Bastian Bechtold)

% AnalysisPt1_Fitting should be run first. The output is loaded in the
%   beginning of AnalysisPt2_Tracking.

% Sample Data available on Zenodo.

%% List of individual .m files and their MATLAB toolbox dependencies
%{
    WorkNotes_AnalysisPt1_Fitting:
        GenCompPhaseMDAFileNameList.m   -> Utility function to
            automatically generate a list of Compensated Phase filenames
            assuming QPM images from Telight's SophiQ software.
            No further dependencies
        FitQPhaseDropletV3.m    -> Main fitting code.
            imgaussfilt.m (Image processing toolbox)
            fit.m (Curve fitting toolbox)
            skewness.m (statistics toolbox?)
            bwconncomp.m (Image processing toolbox)
            regionprops.m (Image processing toolbox)
            fittype.m (Curve fitting toolbox)
            prepareSurfaceData.m (Curve fitting toolbox)
            Heaviside.m (Symbolic math toolbox)
        parfor.m    (Parallel computing toolbox)
            
    WorkNotes_AnalysisPt2_Tracking:
        track.m (Blair and Dufresne)    -> Main tracking code.
        CalcDropletMaxAdjR2FromTracks.m     -> Identify best parameter 
            estimates for each tracked object.
            No further dependencies
        FilterQPhaseParticleList.m      -> Filter objects.
            No futher dependencies
        ametrine.m (Geissbuehler Optics Express 2013)   -> Colormap
        PrepareArrayForViolinPlot.m     -> Utility function
            No further dependencies
        violinplot.m (Bastian Bechtold)     -> Make pretty violin plot.
%}
