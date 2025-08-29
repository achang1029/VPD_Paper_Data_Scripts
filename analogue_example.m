clear; clc;
warning off;
datasets = 'era5'; 
yrbstr = '19500101-20220115';

% 3 bounding boxes from large to small
lonbs = {[-160, -80], [-150, -90], [-140, -100]}; nBd = length(lonbs);
latbs = {[20, 60], [25, 55], [30, 50]};
yrb = [1950,2021];  % bound for years
yrb_clm = [1950, 1999]; % year bound for climatology
yrb_ref = [1950, 1999]; % year bound for analogue reference

nFA = 20;           % analogue number saved 5/10/20/...
optEof = false;     % use EOF for circulation true/false
thEof = 90;         % maximum % variance maintained for EOF, if optEof=true
nEns = 5;           % number of flow analogue ensemble
nMov = 5;           % moving average days 
func = {'euclidean', 'pearson', 'spearman'}; % distance function: Euc/CorPS/CorSM
nFun = length(func);           
nSmoothClm = 31;    % parameters for calculating smooth climatology
dimName = {'Analog Rank (20)', 'Day of Year (365)', 'Year (72)', ...
    'Spatial Range (3, L/M/S)', 'Distance Function (3, Euc/CorPS/CorSM)', 'Pentad Center (5, 1-5)'};  % dimension name

load(sprintf('/data/rong4/Data/ERA5/%s.z500-minusgm.%s.mat', datasets, yrbstr)); % now is data %var = z_mgm;   % use Z500 with global mean removed

% removed 2/29 data first
V = datevec(time);
ind_t = V(:,1)>=yrb(1) & V(:,1)<=yrb(2) & (V(:,2)~=2 | V(:,3)~=29);
time = time(ind_t); V = datevec(time);
yr = min(V(:,1)) : max(V(:,1)); ny = length(yr);
indClm = yr>=yrb_clm(1) & yr<=yrb_clm(2);
indRef = yr>=yrb_ref(1) & yr<=yrb_ref(2);
z = data(:,:,ind_t); [nlon, nlat, nt] = size(z);

% fill nan data at the end and reshape to (nlon,nlat,365,nyear)
if mod(nt, 365) > 0
    z = cat(3, z, nan(nlon,nlat,365-mod(nt,365)));
end
z = reshape(z, nlon, nlat, 365, []);

% start computing
id = nan(nFA, 365, ny, nBd, nFun, nEns); b_cca = nan(nFA, 365, ny, nBd, nFun, nEns);
for ib = 1 : nBd % spatial ranges
    for iF = 1 : nFun % distance function
        disp([ib, iF]);
        [~, ~, ~, id0,b_cca0] = constructedAnalogue(z, lon, lat, nFA, 'indYrRef', find(indRef), 'nEns', nEns, ...
            'optCrop', true, 'lonb', lonbs{ib}, 'latb', latbs{ib}, 'distFunc', func{iF}, ...
            'optMov', true, 'nMov', nMov, 'optNormalize', 'Smooth Standardized Anomaly', ...
            'indClm', find(indClm), 'nSmClm', nSmoothClm, 'optEof', optEof, 'thEof', thEof); %%%optEof true,def:thEof:90
        id(:,:,:,ib,iF,:) = id0; b_cca(:,:,:,ib,iF,:) = b_cca0;
    end
end

outfile = 'test.mat';
save(outfile, 'id', 'b_cca', 'yr', 'yrbstr', 'dimName');


%% prepare surface anomalies
clear; clc;
nMov = 5; nSmoothClm = 31; 
load test.mat

[nFA, ~, ~, nBd, nFun, nEns] = size(id); 

load era5.vpd.19500101-20210619.125W-100W.30N-50N.3htodaily.mat
V = datevec(time); ind = V(:,1)<=2020 & (V(:,2)~=2 | V(:,3)~=29); nlon = length(lon); nlat = length(lat);
data = (data(1:end-1,1:end-1,:) + data(1:end-1,2:end,:) + data(2:end,1:end-1,:) + data(2:end,2:end,:))/4;
vars = reshape(data(:,end:-1:1,ind), [], sum(ind));


nR = 2;
indst = {[1], [2]};
h = matfile('uswest_coast_intermountain.mat'); mask = h.inds; lonci = h.lon; latci = h.lat;% lonr latr
A = get_ind_lonlat(lonci, latci, [-125,-102], [30,50]);
mask = mask(A.ind_lon, A.ind_lat);
[yy, xx] = meshgrid(lat, lon);
wg = cosd(yy);

vars1 = [];
for j = 1 : nR
    ind = ismember(mask(:),indst{j});
    vars1(:,j) = nansum(vars(ind,:) .* wg(ind), 1) ./ sum(wg(ind));
end

indClm = 1 : 50;
nadd = 365 - mod(size(vars1,1)-1, 365) - 1; 
vars1 = cat(1, vars1, nan(nadd, nR));
vO = movmean(vars1, nMov, 1, 'omitnan', 'endpoints', 'fill');      % 5-day moving average
vOm = reshape(movmean(vO, nSmoothClm, 1, 'omitnan', 'endpoints', 'fill'), 365, [], nR);
vOm = nanmean(vOm(:,indClm,:), 2);
vOa = reshape(reshape(vO, 365, [], nR) - vOm, [], nR); 

% % id for days in 2020 which dont' have data are set to NaN
% for example, End date: ERA5 (12-31); MERRA2 (11-30); CFSR (10-05); JRA55 (10-31)
% iEnd = [0, 31, 87, 61]; nhMov = floor(nMov/2);
iEnd = 0;
nhMov = floor(nMov/2);
id = reshape(id, nFA, 365*length(yr), nBd, nFun, nEns);  %
id(:,[1:nhMov,end-iEnd-nhMov+1:end],:,:,:) = nan;  % 
b_cca = reshape(b_cca, nFA, 365*length(yr), nBd, nFun, nEns);
b_cca(:,[1:nhMov,end-iEnd-nhMov+1:end],:,:,:) = nan;  %

%% calculate analogue value
nYr = size(id, 2)/365;
vFMd = nan(365*nYr, nBd, nFun, nEns, nR);  % median analogue
vCCA = vFMd;  % constructed analogue
for iR = 1 : nR     % region
    tmp1 = vOa(:,iR);
    disp(iR);
    for iB = 1 : nBd       % bounding box
        for iF = 1 : nFun  %%%%%%%%%%% distance function
            for iE = 1 : nEns  % analogue ensemble
                ind = any(isnan(id(:,:,iB,iF,iE)), 1); %%%%%%%%%%%%%%%%
                tmp2 = tmp1(id(:,~ind,iB,iF,iE));
                B = b_cca(:,~ind,iB,iF,iE);
                vCCA(~ind,iB,iF,iE,iR) = sum(B .* tmp2, 1);
                vFMd(~ind,iB,iF,iE,iR) = median(tmp2(:,:), 1, 'omitnan');
            end
        end
    end
end
vO = single(vO); vOa = single(vOa); vOm = single(vOm);
vFMd = single(vFMd); 
vCCA = single(vCCA);