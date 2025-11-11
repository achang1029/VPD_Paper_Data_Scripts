function out = qdm_z(obs, mod_cp, var, frq, trace, jitter_factor)
% This function performs QDM bias correction and was combined and modified 
% from ClimQMBC (Matlab) and ClimDown (R)
% The codes with moving projected model windows are adapted from ClimQMBC
% The core code for non-parametric distribution fitting and correction are
% adapted from ClimDown. (while ClimQMBC uses parametric fitting)
%
% Input:
%   obs = observed data (nlon x nlat x nt_c); nt_c: length of historical
%   period (reference for correction)
%
%   mod_cp = modeled data (nlon x nlat x nt_cp); nt_cp = nt_c+nt_p, length
%   of model period, the first nt_c should represent the same time period
%   as observation
%
%   var = A flag that identifies if data type.
%         This flag tells the getDist function if it has to discard
%         distribution functions that allow negative numbers, and if the 
%         terms in the correction equations are multiplied/divided or
%         added/subtracted.
%         Temperature (non-bounded): var = 0, correction is done
%         added/substracted
%         Precipitation (bounded, positive, ...): var = 1, correction is
%         done multiplied/divided
%
%   NOTE: This routine considers that obs and mod series start in january
%         of the first year and ends in december of the last year.
%
%   frq = A string specifying if the input is annual or monthly data. If
%         not specified, it is set monthly as default.
%         Monthly:    frq = 'M'
%         annual:     frq = 'A'
% 
%   trace = A float indicating the threshold to consider physically 
%                  null precipitation values. (when var == 1)
%
%   jitter_factor = A float indicating the maximum value of the random values
%               that added to original data to acommadate ties.
%
%   rel_change_th = A float indicating the maximum scaling factor (Equation
%                   4 of Cannon et al. (2015)) when the denominator is
%                   below inv_mod_th.
%
%   inv_mod_th = A float indicating the upper threshold of the denominator
%                of the scaling factor (Equation 4 of Cannon et al. (2015))
%                to truncate the scaling factor. This parameter is defined
%                as default as the pp_threshold parameter, described above.
%
% Output:
%   out = QDM corrected modeled data (same size as mod_cp)
%
% References:
%   Cannon, A. J., et al., 2015: Bias correction of
%   GCM precipitation by quantile mapping: How well do methods preserve
%   changes in quantiles and extremes? J. Climate, 28(17), 6938-6959,
%   https://doi.org/10.1175/JCLI-D-14-00754.1
%
%   ClimDown: https://github.com/pacificclimate/ClimDown
%   ClimQMBC: https://github.com/saedoquililongo/climQMBC

% Yizhou Zhuang, UCLA
% May 15, 2024
%%
% Define default values for arguments
if ~exist('trace','var')
    trace = 0.05;
end

if ~exist('jitter_factor','var')
    jitter_factor = 0.01;
end

if ~exist('frq', 'var')
    frq = 'M';
end

if strcmp(frq, 'M')
    I = 12;
elseif strcmp(frq, 'Y')
    I = 1;
elseif strcmp(frq, 'D')
    I = 365;
elseif strcmp(frq, 'X') % for 360-day calendar models
    I = 360;
else
    error('undefined frq option!');
end
%% 
% Check if input (observation or model) is all NaN
if all(isnan(mod_cp), 'all') || all(isnan(obs), 'all')
    out = nan(size(mod_cp));
    return;
end

% Check if spatial size is matched between observation and models
N_obs = size(obs); N_mod = size(mod_cp); 
if N_obs(1)~=N_mod(1) || N_obs(2)~=N_mod(2)  
    error('spatial dimension not matched!');
end

% Remove grid point with all NaN data
n_ll = N_obs(1)*N_obs(2);
mask = ~(all(isnan(obs), 3:ndims(obs)) | all(isnan(mod_cp), 3:ndims(mod_cp)));
n_llm = sum(mask(:)); N = n_llm * I;
obs = reshape(single(obs), n_ll, I, []);
mod_cp = reshape(single(mod_cp), n_ll, I, []);
obs = reshape(obs(mask(:), :, :), N, []); ny_obs = size(obs, 2);
mod_cp = reshape(mod_cp(mask(:), :, :), N, []); ny_mod = size(mod_cp, 2);

% Apply jitter to accommodate ties and measurement precision
obs = obs + jitter_factor * rand(size(obs));
mod_cp = mod_cp + jitter_factor * rand(size(mod_cp));

% Handle ratio data for zero values as left censored
if var == 1
    epsilon = eps('double');
    ind = obs < trace & ~isnan(obs);
    obs(ind) = epsilon + (trace-epsilon) * rand(sum(ind(:)), 1);
    ind = mod_cp < trace & ~isnan(mod_cp);
    mod_cp(ind) = epsilon + (trace-epsilon) * rand(sum(ind(:)), 1);
end

% Calculate quantiles of observation and models during reference period
tau = (1 : ny_obs)/(ny_obs + 1);
quant_obs = quantile(obs, tau, 2); %iosr.statistics.
mod_c = mod_cp(:,1:ny_obs);
quant_mod_c = quantile(mod_c, tau, 2); %iosr.statistics.

% Bias correction for model reference period
mhat_c = nan(N, ny_obs);
parfor i = 1 : N
    %Everything inside this for loop is from ChatGPT to address
    %uniqueness problem with griddedInterpolant
    x = quant_mod_c(i,:);
    [x_unique, ia] = unique(x, 'stable');
    y_unique = quant_obs(i, ia);

    if numel(x_unique) < 2 || any(isnan(x_unique)) || any(isnan(y_unique))
        mhat_c(i,:) = nan;
    else
        F = griddedInterpolant(x_unique, y_unique, 'linear', 'nearest');
        mhat_c(i,:) = F(mod_c(i,:));
    end
end

% Bias correction for model projected period 
if var == 1
    F1 = griddedInterpolant(tau', quant_mod_c', 'linear', 'nearest');
    F2 = griddedInterpolant(tau', quant_obs', 'linear', 'nearest');
else
    FD = griddedInterpolant(tau', quant_obs'-quant_mod_c', 'linear', 'nearest');
end
% construct projected windows with the same length as the current period
ind = (1:ny_obs)'+(1:ny_mod-ny_obs);
mod_P = reshape(mod_cp(:,ind(:)), N, ny_obs, []);
quant_mod_P = quantile(mod_P, tau, 2);  %iosr.statistics.
mod_P1 = squeeze(mod_P(:,end,:));

adj = nan(N, ny_mod-ny_obs);
if var == 1
    parfor j = 1 : ny_mod-ny_obs
        F_mod_p1 = nan(N,1);
        for i = 1 : N
            %Everything inside this for loop is from ChatGPT to address
            %uniqueness problem with griddedInterpolant
            x = quant_mod_P(i,:,j);
            [y, ia] = unique(x, 'stable');  % Ensure uniqueness
            
            if numel(y) < 2 || any(isnan(y)) || any(isnan(tau(ia)))
                F_mod_p1(i) = nan;
            else
                F = griddedInterpolant(y, tau(ia), 'linear', 'nearest');
                F_mod_p1(i) = F(mod_P1(i,j));
            end
        end
        adj(:,j) = diag(F2(F_mod_p1)) ./ diag(F1(F_mod_p1));
    end
    mhat_p = mod_P1 .* adj;
else
    parfor j = 1 : ny_mod-ny_obs
        F_mod_p1 = nan(N,1);
        for i = 1 : N
            %Everything inside this for loop is from ChatGPT to address
            %uniqueness problem with griddedInterpolant
            x = quant_mod_P(i,:,j);
            [y, ia] = unique(x, 'stable');  % Ensure uniqueness
            
            if numel(y) < 2 || any(isnan(y)) || any(isnan(tau(ia)))
                F_mod_p1(i) = nan;
            else
                F = griddedInterpolant(y, tau(ia), 'linear', 'nearest');
                F_mod_p1(i) = F(mod_P1(i,j));
            end
        end
        adj(:,j) = diag(FD(F_mod_p1));
    end
    mhat_p = mod_P1 + adj;
end

mhat = [mhat_c, mhat_p];
if var == 1
    mhat(mhat < trace) = 0;
end
% mhat = reshape(mhat, N_mod);

out = nan(n_ll, I, ny_mod);
out(mask(:),:,:) = reshape(mhat, n_llm, I, []);
out = reshape(out, N_mod);
end