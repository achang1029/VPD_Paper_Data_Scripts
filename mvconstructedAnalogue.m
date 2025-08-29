% flow_analogue flow analogue search
%  [DIST, IDOY, IYR, ID] = flow_analogue(X, Lon, Lat, indYrRef, nFA) 
%  returns analogue day index and distance between analogue and observation.
%  Lon and Lat are the longitudes and latitudes in a vector; X is the input
%  circulation pattern variable array (e.g. geopotential) with the size of
%  [N(lon), N(lat), N(day-of-year), N(year)]; indYrRef is the reference
%  year index for searching analogues, e.g. 1:32; nFA is the number of flow
%  analogue returned to output. Ouput analogue indices include ID, IDOY, 
%  and IYR: IDOY and IYR are index of DOY (day of year) and year, and ID  
%  equals IDOY + (IYR-1)*N(day-of-year). All output has the same size of 
%  [N(day-of-year), N(year), nFA, N(ensemble)].
%  
%  [DIST, IDOY, IYR, ID] = flow_analogue(..., OPTIONAL_ARGUMENT, ARG_VALUE)
%  Optional arguments include:
%    nd: half window size for analogue search, window size is nd*2+1.
%      Default is 30;
%
%    distFunc: method for calculating distance, include "Pearson", 
%      "Spearman", and"Euclidean" (default)
%
%    optMov: option to perform moving average, true (default) or false. 
%    nMov: moving average size. Default is 5.
%
%    nEns: flow analogue ensemble size; for each ensemble, analogue is 
%      search every nEns days. Default is 5.
%
%    optCrop: option to crop certain region, true or false (default).
%    lonb: bounding box for longitude, ie. [-160,-80]
%    latb: bounding range for latitude, ie. [20, 60]
%
%    optEof: option to perform EOF and then use reconstructed data as 
%      input, true (default) or false.
%    thEof: cutoff explained variance (%) for EOF. Default is 90.

%    indYrRef: index of years to search for analogues, default is 1:N(year)

%    optNormalize: option to normalize data, include "none", "anomaly", 
%      "standardized anomaly", "smooth anomaly", and "smooth standardized 
%      anomaly" (default).
%    indClm: index of years for calculating climatology. Default is the 
%      same as the reference period (indYrRef)
%    nSmClm: window size for calculating smooth climatology. Default is 31.


%  Yizhou Zhuang (zhuangyz@atmos.ucla.edu)
%  References:
%    Zhuang et al. 2021 (preparing)
%    Yiou et al. 2007 (GRL)
%
function [dist, iDoy, iYr, id, b_cca] = mvconstructedAnalogue(X, lon, lat, nFA, varargin)
%% check input size
if ndims(X) ~= 5
    error('Input X must be a 5-dimension array!');
end
nlon = size(X, 1); nlat = size(X, 2);
if nlon ~= length(lon) || nlat ~= length(lat)
    error('Size mismatch between input X and lon/lat!');
end
nVar = size(X, 3);
nDoy = size(X, 4);  % number of days in year, must be <= 365
if nDoy > 365
    error('Days of a year exceed 365 days!');
end
nYr = size(X, 5);   % number of years
%% parameters with default value
nd = 30;  % calendar window for analogue search is 2*nd+1 days centered on current day

optMov = true; % option to perform moving average on input data X
nMov = 5; % moving average length

nEns = 5; % search analogue every nEns days, creating nEns ensembles

optCrop = false; % option to use certain region for flow analogue search
lonb = bounds(lon); % bounding range for longitude, ie. [-160, -80]
latb = bounds(lat); % bounding range for latitude, ie. [20, 60]

optEof = true;  % perform EOF on data X and then use reconstructed X as input
thEof = 90; % cutoff explained variance for EOF

% indYrRef = 1 : 32;   % use 1979-2010 (input is 1979-2020)
indYrRef = 1 : nYr;  % use all years

optNormalize = "smooth standardized anomaly";
indClm = indYrRef;  % year index for calculating climatology; default is the same as reference period
nSmClm = 31; % window size for calculating smooth climatology

distFunc = "Euclidean"; % "Pearson", "Spearman", "Euclidean"
%% override the default parameters if specified
if nargin < 4
    error('At least 4 input arguments are needed: X, lon, lat, nFA.');
elseif nargin > 4 && mod(nargin,2)==0
    i = 1;
    while i < nargin-4
        switch lower(varargin{i})
            case "nd"
                nd = varargin{i+1};
            case "optmov"
                optMov = varargin{i+1};
            case "nmov"
                nMov = varargin{i+1};
            case "nens"
                nEns = varargin{i+1};
            case "optcrop"
                optCrop = varargin{i+1};
            case "lonb"
                lonb = varargin{i+1};
            case "latb"
                latb = varargin{i+1};
            case "opteof"
                optEof = varargin{i+1};
            case "theof"
                thEof = varargin{i+1};
            case "indyrref"
                indYrRef = varargin{i+1};
            case "optnormalize"
                optNormalize = varargin{i+1};
            case "indclm"
                indClm = varargin{i+1};
            case "nsmclm"
                nSmClm = varargin{i+1};
            case "distfunc"
                distFunc = varargin{i+1};
        end
        i = i + 2;
    end
end
%% Optional pre-processing
% perform moving average for input X
if optMov
    if nDoy == 365
        X = movmean(reshape(X, nlon, nlat, nVar, []), nMov, 4, ...
            'includenan', 'endpoints', 'fill');
    elseif nDoy < 365
        X = reshape(movmean(X, nMov, 4, 'includenan', 'endpoints', 'fill'), ...
            nlon, nlat, nVar, []);
    end
end
% crop certain region of X for searching flow analogue
if optCrop
    [X, lon, lat] = crop(X, lon, lat, lonb, latb);
    nlon = length(lon); nlat = length(lat);
end
% Option to normalize X data
if optNormalize ~= "none"
    X = normalizeData(X);
end
% perform EOF reconstruction
if optEof
    [X, P] = eofReconstruct(X);
    nCCA = size(P, 2); % use EOF number to repsent CCA number
    if nFA < nCCA
        disp('Warning! Analogue number is less than EOF number.');
    end
else
    nCCA = nFA;
end
%% Data preparation and output initialization
% prepare input array X for flow analogue searching
X = reshape(X, nlon*nlat*nVar, nDoy, []); X0 = reshape(X, nlon*nlat*nVar, []);
X = cat(2, cat(3, nan(nlon*nlat*nVar,nd,1), X(:,nDoy-nd+1:end,1:end-1)), ...
    X, cat(3, X(:,1:nd,2:end), nan(nlon*nlat*nVar,nd,1)));
if optEof
    nEof = size(P, 2);
    P = reshape(P', nEof, nDoy, nYr); P0 = reshape(P, nEof, []);
    P = cat(2, cat(3, nan(nEof,nd,1), P(:,nDoy-nd+1:end,1:end-1)), ...
        P, cat(3, P(:,1:nd,2:end), nan(nEof,nd,1)));
end
% create a temporary time index array
ti = reshape(1 : nDoy*nYr, nDoy, []);
if nDoy == 365
    ti = [[nan(nd,1), ti(nDoy-nd+1:end,1:end-1)]; ti; ...
        [ti(1:nd,2:end), nan(nd,1)]];
else
    ti = [nan(nd,nYr), ti, nan(nd,nYr)];
end
% initialize distance array
dist = nan(nFA, nDoy, nYr, nEns);  
id = dist; iDoy = dist; iYr = dist;
b_cca = zeros(nFA, nDoy, nYr, nEns);
%% search flow analogue for each calendar day
for i = 1 : nDoy
    % temporary time index array for analogue ID assignment
    ti1 = ti(i:i+2*nd,:); ti1 = ti1(:);
    % calculate distance between current day and all possible analogue days
    switch lower(distFunc)
        case "pearson"
            X1 = reshape(X(:,i:i+2*nd,:), nlon*nlat*nVar, []);
            dist0 = 1 - corr(X1, squeeze(X(:,i+nd,:)));
        case "spearman"
            X1 = reshape(X(:,i:i+2*nd,:), nlon*nlat*nVar, []);
            dist0 = 1 - corr(X1, squeeze(X(:,i+nd,:)), 'type', 'spearman');
        case "euclidean"
            if optEof
                dist0 = pdist2(reshape(P(:,i:i+2*nd,:), nEof, [])', ...
                    squeeze(P(:,i+nd,:))');
            else
                dist0 = pdist2(reshape(X(:,i:i+2*nd,:), nlon*nlat*nVar, [])', ...
                    squeeze(X(:,i+nd,:))');
            end
    end
    % sort distance to current calendar day of all years to search analogue
    [dist0_sorted, dist0_ind] = sort(dist0, 1, 'ascend', 'missingplacement', 'last');
    % get index of year
    dist0_indYr = ceil(dist0_ind / (2*nd+1));
    % current year and years not in the reference period will not be considered
    FAindYr = dist0_indYr ~= repmat(1:nYr,(2*nd+1)*nYr,1) & ...
        ismember(dist0_indYr, indYrRef); 
    % loop for all years
    for j = 1 : nYr
%         if j==6;
%             1;
%         end
        dist0 = dist0_sorted(FAindYr(:,j), j);
        id0 = ti1(dist0_ind(FAindYr(:,j), j));
        if all(isnan(dist0))
            dist(:,i,j,:) = nan;
            id(:,i,j,:) = nan;
            b_cca(:,i,j,:) = nan;
        else
            for k = 1 : nEns
                ind = find(mod(id0,nEns)==k-1);
                dist(:,i,j,k) = dist0(ind(1:nFA));
                id(:,i,j,k) = id0(ind(1:nFA));
                % CCA
                nCCA1 = nCCA;
				%{
                if optEof
                    while rank(X0(:,id(1:nCCA1,i,j,k))) < nCCA
                        nCCA1 = nCCA1 + 1;
                    end
                end
				%}
                B = regress(X0(:,(j-1)*nDoy+i), X0(:,id(1:nCCA1,i,j,k)));
                b_cca(1:nCCA1,i,j,k) = B;
            end
        end
    end
    iDoy = mod(id-1, nDoy) + 1;
    iYr = ceil(id/nDoy); 
end

%% Functions 
    function X = normalizeData(X)
        if contains(lower(optNormalize), "smooth")
            Xm = reshape(movmean(X, nSmClm, 4, 'includenan', 'endpoints', 'fill'), ...
                nlon, nlat, nVar, nDoy, []);
            Xm = nanmean(Xm(:,:,:,:,indClm), 5);
            Xs = reshape(movmean(X.^2, nSmClm, 4, 'includenan', 'endpoints', 'fill'), ...
                nlon, nlat, nVar, nDoy, []);
            Xs = nanmean(Xs(:,:,:,:,indClm), 5);
            n1 = length(indClm) * nSmClm;
            Xs = sqrt((Xs - Xm.^2) * n1/(n1-1));
            X = reshape(X, nlon, nlat, nVar, nDoy, []);
        else
            X = reshape(X, nlon, nlat, nVar, nDoy, []);
            Xm = mean(X(:,:,:,:,ind_clm), 5);
            Xs = nanstd(X(:,:,:,:,ind_clm), 0, 5);
        end
        switch lower(optNormalize)
            case {"anomaly", "smooth anomaly"}
                X = (X - Xm);
            case {"standardized anomaly", "smooth standardized anomaly"}
                X = (X - Xm) ./ Xs;
        end
    end
    % Perfomr EOF reconstruction, output reconstructed X and PCs
    function [X, P] = eofReconstruct(X)
        [yq, ~] = meshgrid(lat, lon);
        w = sqrt(cosd(yq));  % grid weighting
        w = repmat(w, [1,1,nVar]);
        [E, P] = eof(reshape(X, nlon*nlat*nVar, []) .* w(:), thEof);
        % reconstruct input array X with EOFs and PCs
        X = E * P' ./ w(:);
    end
end
%% More functions
% Calculate EOF and retain first few modes according to set minimum
% explained variance
function [eof, pc] = eof(data, thEof)
[np, nt] = size(data);
% remove time instant at which data are all NaN
indt_nan = all(isnan(data), 1);
data = data(:, ~indt_nan);
% remove grid points at which data are all NaN
indp_nan = all(isnan(data), 2);
data = data(~indp_nan, :);
% perform EOF using PCA function
[eof0, pc0, ~, ~, expl] = pca(data', 'centered', false);
ne = find(cumsum(expl) >= thEof, 1, 'first'); % number of mode to retain
eof = nan(np, ne);
eof(~indp_nan,:) = eof0(:, 1:ne);
pc = nan(nt, ne);
pc(~indt_nan,:) = pc0(:, 1:ne);
end

% Cropping input data according to bounding box of longitude and latitude
function [X, lon, lat] = crop(X, lon, lat, lonb, latb)
% change scale of longitude bound to match longitude data
if any(lon<0) && all(lonb>=180)
    lonb = lonb - 360;
elseif all(lon>=0) && all(lonb<0)
    lonb = lonb + 360;
end
ind_lon = lon>=lonb(1) & lon<=lonb(2);
ind_lat = lat>=latb(1) & lat<=latb(2);
lon = lon(ind_lon); 
lat = lat(ind_lat); 
X = X(ind_lon, ind_lat, :,:);
end