function WorkNotes2022_01_12_AnalysisPt1_Fitting

%%% Adapted from WorkNotes2021_12_15_Basusree_LIPS6.m %%%

%% 1.0 Setup directories

FileName(1).Dir = '/Users/mccall/ownCloud/_SharedFolders/2024-01-16_PoL_Master/Ln3_a_zstack_dz0-5um_1FOV_20C_40x';

nFOVvec = [1];
nZvec =   [41];
nTvec =   [1];
Prefix = {'S0188g_FOV1'};

nFile = length(FileName);
for i = 1:nFile
    FileName(i).Prefix = Prefix{i};     %FileName(i).Dir(53:58);
    FileName(i).nFOV = nFOVvec(i);   %FileName(i).Dir(113);    
    FileName(i).nZ = nZvec(i);    %FileName(i).Dir(104:105);
    FileName(i).nT = nTvec(i);
    FileName(i).MG = 'G1';
end

%% 2.0 Run fits to 2 models over all images in zstacks - domain restriction, write-after-each-plane-structure

% Run with parfor on nFOV
tic
DS = cell(nFile,1);
for f = 1:nFile
    cd(FileName(f).Dir);
    nFOV = FileName(f).nFOV;
    nZ = FileName(f).nZ;
    nT = FileName(f).nT;
    MicroscopeGeneration = FileName(f).MG;
    FileNameList = GenCompPhaseMDAFileNameList(nFOV,nZ,nT,{[],[],[]},MicroscopeGeneration);
    tmp = cell(nFOV,1);
    parfor i = 1:nFOV
        disp(['Processing FOV ',num2str(i),' of ',num2str(nFOV)])
        DSV3_tmp = cell(nZ,1);
        for j = 1:nZ 
            disp(['Processing plane ',num2str(j),' of ',num2str(nZ)])
            FileNameListInd = nZ*(i-1) + j;
            DataFileName = FileNameList{FileNameListInd};
            
            % Domain Restriction off
            AnalysisOption = {[4],[],'SelfAdjusted',3,0,5,0.7,...
                {'Sphere','SphericalCap_fixPhi0_RlessZmidPos'},...
                {nan(1,4),nan(1,6)},...
                0,1E5,[],[]};
            PlottingOption = [0,0];
            [DSV3,~,~] = FitQPhaseDropletV3(DataFileName,AnalysisOption,PlottingOption);
            nParticle = length(DSV3(:,1));
            DSV3_tmp{j} = [DSV3,j*ones(nParticle,1),i*ones(nParticle,1)];
        end
        tmp{i,1} = DSV3_tmp;
    end
    DS{f,1} = tmp;
end
toc
%%% Conclusions %%%
% 1. Elapsed time is 352.935903 seconds.
% 2. Output format is DS{f}{FOV}{Z}

%% save data
cd '/Users/mccall/Desktop/AnalysisCode_QPM/SampleData'
%Type save command in command line
%save('WorkNotes2022-01-12_SampleAnalysisPt1_Fitting.mat')
