function violins = violinplot(data, cats, varargin)
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the
%   cell-of-strings CATEGORYNAMES.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA.
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Defaults to 0.3.
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false
%     'GroupOrder'   Cell of category names in order to be plotted.
%                    Defaults to alphabetical ordering

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license:

%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met:

%   1. Redistributions of source code must retain the above copyright 
%   notice, this list of conditions and the following disclaimer.

%   2. Redistributions in binary form must reproduce the above copyright 
%   notice, this list of conditions and the following disclaimer in the 
%   documentation and/or other materials provided with the distribution.

%   3. Neither the name of the copyright holder nor the names of its 
%   contributors may be used to endorse or promote products derived from 
%   this software without specific prior written permission.

%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
%   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
%   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
%   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
%   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%   Description of BSD 3-clause license obtained from 
%   https://opensource.org/licenses/BSD-3-Clause on 2022.01.12

%   Original source code available from
%   https://github.com/bastibe/Violinplot-Matlab

%%
    hascategories = exist('cats','var') && not(isempty(cats));
    
    %parse the optional grouporder argument 
    %if it exists parse the categories order 
    % but also delete it from the arguments passed to Violin
    grouporder = {};
    idx=find(strcmp(varargin, 'GroupOrder'));
    if ~isempty(idx) && numel(varargin)>idx
        if iscell(varargin{idx+1})
            grouporder = varargin{idx+1};
            varargin(idx:idx+1)=[];
        else
            error('Second argument of ''GroupOrder'' optional arg must be a cell of category names')
        end
    end

    % tabular data
    if isa(data, 'dataset') || isstruct(data) || istable(data)
        if isa(data, 'dataset')
            colnames = data.Properties.VarNames;
        elseif istable(data)
            colnames = data.Properties.VariableNames;
        elseif isstruct(data)
            colnames = fieldnames(data);
        end
        catnames = {};
        for n=1:length(colnames)
            if isnumeric(data.(colnames{n}))
                catnames = [catnames colnames{n}];
            end
        end
        for n=1:length(catnames)
            thisData = data.(catnames{n});
            violins(n) = Violin(thisData, n, varargin{:});
        end
        set(gca, 'xtick', 1:length(catnames), 'xticklabels', catnames);

    % 1D data, one category for each data point
    elseif hascategories && numel(data) == numel(cats)
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end

        catnames = categories(cats);
        for n=1:length(catnames)
            thisCat = catnames{n};
            thisData = data(cats == thisCat);
            violins(n) = Violin(thisData, n, varargin{:});
        end
        set(gca, 'xtick', 1:length(catnames), 'xticklabels', catnames);

    % 1D data, no categories
    elseif not(hascategories) && isvector(data)
        violins = Violin(data, 1, varargin{:});
        set(gca, 'xtick', 1);

    % 2D data with or without categories
    elseif ismatrix(data)
        for n=1:size(data, 2)
            thisData = data(:, n);
            violins(n) = Violin(thisData, n, varargin{:});
        end
        set(gca, 'xtick', 1:size(data, 2));
        if hascategories && length(cats) == size(data, 2)
            set(gca, 'xticklabels', cats);
        end

    end

end
