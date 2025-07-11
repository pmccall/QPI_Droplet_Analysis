function ViolinReadyArray = PrepareArrayForViolinPlot(DataSet)

% Written by PMM on 2019.11.25

%%% PURPOSE %%%
% Generatate an array, the columns of which represent the different
%   distributions to be compared in a violin plot, with NaN padding such
%   that they are all the same length.

%%% INPUTS %%%

% DataSet is an n-element cell array. Each element is an m-by-1 vector,
%   where m may be different for each element.

%%% OUTPUTS %%%
% ViolinReadyArray is a k-by-n array, where k is the lenth of the longest
%   vector in DataSet. The ith column of ViolinReadyArray contains the 
%   contents of the ith elements of DataSet, with NaN padding at the end
%   (as necessary) to bring the total number of elements in the column up
%   to k.

%% Code Start

% Check input type
if iscell(DataSet) == 0
    disp('DataSet not of expected to by of type cell')
    ViolinReadyArray = [];
    return
end

% Determine max column length
nCol = length(DataSet(:));
ColLength = zeros(nCol,1);
for i = 1:nCol
    ColLength(i) = length(DataSet{i});
end
k = max(ColLength);

% Generate Output array
ViolinReadyArray = NaN(k,nCol);
for i = 1:nCol
    ViolinReadyArray(1:ColLength(i),i) = DataSet{i};
end
end