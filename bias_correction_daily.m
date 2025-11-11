projectdir = '/home/achang/CMIP6';
% dinfo = dir(fullfile(projectdir));
% dinfo([dinfo.isdir]) = [];     %get rid of all directories including . and ..
% nfiles = length(dinfo);

models = load(append(projectdir, '/model_list.mat')).models;
experiments = ["ssp126", "ssp245", "ssp585"];
variables = ["hurs", "tasmin", "tasmax"];

t2_daily = load(append(projectdir, '/t2_daily.mat')).t2;
d2_daily = load(append(projectdir, '/d2_daily.mat')).d2;
t2max_daily = load(append(projectdir, '/t2max_daily.mat')).t2max;
t2min_daily = load(append(projectdir, '/t2min_daily.mat')).t2min;

hurs_daily = d2m_to_hurs(t2_daily, d2_daily);
for h = 1:length(variables)
    for i = 1:length(models)
        for k = 1:length(experiments)
            try
                disp(variables(h))
                disp(models(i)) 
                disp(experiments(k))
                
                if strcmp(models(i),'HadGEM3-GC31-LL')
                    ens = 'r1i1p1f3';
                    per_year = 360;
                elseif strcmp(models(i), 'KACE-1-0-G')
                    ens = 'r1i1p1f1';
                    per_year = 360;
                else
                    ens = 'r1i1p1f1';
                    per_year = 365;
                end

                fname = string(append(projectdir, '/regridded/r1i1p1f1_daily/', variables(h), '_day_', models(i), '_', experiments(k), '_', ens , '_regridded.mat'));
                fname_historical = string(append(projectdir, '/regridded/r1i1p1f1_daily/', variables(h), '_day_', models(i), '_historical_', ens , '_regridded.mat'));
                experim = load(fname);
                historical = load(fname_historical);
                experim = experim.(string(fieldnames(experim)));
                historical = historical.(string(fieldnames(historical)));

                if variables(h) == 'hurs'
                    if isfile(string(append(projectdir, '/bias_corrected/r1i1p1f1_daily/ea_day_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')))
                        disp("Already bias corrected. Continuing to next.")
                        continue
                    end
                    obs = d2m_to_svp(d2_daily);

                    fname_tasmax_historical = string(append(projectdir, '/regridded/r1i1p1f1_daily/tasmax_day_', models(i), '_historical_', ens , '_regridded.mat'));
                    fname_tasmin_historical = string(append(projectdir, '/regridded/r1i1p1f1_daily/tasmin_day_', models(i), '_historical_', ens , '_regridded.mat'));
                    tasmax_historical = load(fname_tasmax_historical);
                    tasmin_historical = load(fname_tasmin_historical);
                    tasmax_historical = tasmax_historical.(string(fieldnames(tasmax_historical)));
                    tasmin_historical = tasmin_historical.(string(fieldnames(tasmin_historical)));

                    fname_tasmax_experim = string(append(projectdir, '/regridded/r1i1p1f1_daily/tasmax_day_', models(i), '_', experiments(k), '_', ens , '_regridded.mat'));
                    fname_tasmin_experim = string(append(projectdir, '/regridded/r1i1p1f1_daily/tasmin_day_', models(i), '_', experiments(k), '_', ens , '_regridded.mat'));
                    tasmax_experim = load(fname_tasmax_experim);
                    tasmin_experim = load(fname_tasmin_experim);
                    tasmax_experim = tasmax_experim.(string(fieldnames(tasmax_experim)));
                    tasmin_experim = tasmin_experim.(string(fieldnames(tasmin_experim)));                     

                    % historical_max = hurs_to_d2m(tasmax_historical, historical);
                    % experim_max = hurs_to_d2m(tasmax_experim, experim);
                    % 
                    % historical_min = hurs_to_d2m(tasmin_historical, historical);
                    % experim_min = hurs_to_d2m(tasmin_experim, experim);

                    historical = (d2m_to_svp(hurs_to_d2m(tasmax_historical, historical)) + d2m_to_svp(hurs_to_d2m(tasmin_historical, historical))) ./ 2;
                    experim = (d2m_to_svp(hurs_to_d2m(tasmax_experim, experim)) + d2m_to_svp(hurs_to_d2m(tasmin_experim, experim))) ./ 2;

                    var = 0;
                elseif variables(h) == 'tasmin'
                    if isfile(string(append(projectdir, '/bias_corrected/r1i1p1f1_daily/', variables(h), '_day_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')))
                        disp("Already bias corrected. Continuing to next.")
                        continue
                    end
                    obs = t2min_daily;
                    var = 0;
                elseif variables(h) == 'tasmax'
                    if isfile(string(append(projectdir, '/bias_corrected/r1i1p1f1_daily/', variables(h), '_day_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')))
                        disp("Already bias corrected. Continuing to next.")
                        continue
                    end
                    obs = t2max_daily;
                    var = 0;
                end
                %Reference period is 1979-2000. obs from 1979-2023,
                %historical from 1850-2014, future 2015 onwards
                obs_shortened = obs(:,:,1:((2001-1979)*per_year));
                historical_experim = cat(3, historical(:,:,(1979-1850)*per_year+1:(2015-1850)*per_year), experim);

                if per_year == 360
                    frq = 'X';
                else
                    frq = 'D';
                end
                
                if variables(h) == 'hurs'
                %     corrected_max = qdm_z(obs_shortened, historical_experim, experim_max), var);
                %     save(string(append(projectdir, '/bias_corrected/d2max_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')), 'corrected_max')                    
                %     corrected_min = qdm_z(obs_shortened, historical_experim, experim_min), var);
                %     save(string(append(projectdir, '/bias_corrected/d2min_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')), 'corrected_min')                    
                    corrected = qdm_z(obs_shortened, historical_experim, var, frq);
                    save(string(append(projectdir, '/bias_corrected/r1i1p1f1_daily/ea_day_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')), 'corrected')
                else
                    corrected = qdm_z(obs_shortened, historical_experim, var, frq);
                    save(string(append(projectdir, '/bias_corrected/r1i1p1f1_daily/', variables(h), '_day_', models(i), '_', experiments(k), '_', ens , '_corrected.mat')), 'corrected')
                end
                disp("Complete. Continuing to next...")
                % disp("")
            catch exception
                msgtext = getReport(exception);
                disp(msgtext)
                disp("File not found. Continuing...")
                % disp("")
                continue
            end
        end
    end
end