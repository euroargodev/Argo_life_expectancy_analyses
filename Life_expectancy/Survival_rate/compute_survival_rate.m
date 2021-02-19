% plot_survival_rate
% plots survival rates for a list of floats given in number of cycles, 
% vertical_km and float age
%
% Input
% - list_dir: path to floats list. Floats list should contains at least WMO
%       and STATUS (list from jcommops web site)
% - dac_dir: path to gdac
% - output_folder: path to output folder
%
% Output
% - Figures with survival rates depending on number of cycles, vertical km 
%       and float age 
% - Figures with sample histograms of number of cycles, vertical km and
%       float age
%
% Auxiliary functions:
%    read_csv
%    get_floats_filespath
%    get_vertical_km_multiprof
%    format_data_for_plotting
%
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/13



% add paths (packages and auxiliary functions)
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/M_map
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/MyTools
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/seawater
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/addaxis % more than 2 y axis in a figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/legendflex % more than 1 column in legend
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/setgetpos_V1.2
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/ezyfit/ezyfit % fit data to a curve
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/aux_functions
% clear variables
close all
clear variables

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_dir = '/home1/datahome/co_arg/larduini/Lists/Comparison_EA_Intern/Intern_Iridium_open_ocean.csv'
%list_dir = '/home1/datahome/co_arg/agarciaj/life_expentancy/lists/ARVOR_globalocean_from2018.csv'
dac_dir = '/home/ref-argo/gdac/'
output_folder = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% get floats from list
[Floats_list] = read_csv(list_dir,',');

% number of floats
Floats_list.WMO = cellstr(Floats_list.REF);

% when there is a not valid WMO in list
index_wmo = cellfun(@(x) length(x)~=7, Floats_list.WMO);
Floats_list_field = fieldnames(Floats_list);
for ifield = 1: length(Floats_list_field)
    Floats_list.(Floats_list_field{ifield})(index_wmo,:) = [];
end

% look for dac path
% get dac
[Floats_paths] = get_floats_filespath(dac_dir, Floats_list.WMO);
% n_floats = length(Floats_list.WMO.data);
index_dac = ~ismember(cellstr(Floats_list.WMO), Floats_paths.WMO);
for inames = 1 : length(Floats_list_field)
    Floats_list.(Floats_list_field{inames})(index_dac,:) = [];
end


% sort by wmo before merge structs
[Floats_list.WMO, index_sort] = sort(Floats_list.WMO);
Floats_list_field(contains(Floats_list_field, 'WMO')) = [];
for inames = 1 : length(Floats_list_field)
    Floats_list.(Floats_list_field{inames}) = Floats_list.(Floats_list_field{inames})(index_sort,:);
end
[Floats_paths.WMO, index_sort] = sort(Floats_paths.WMO);
Floats_list.DAC = Floats_paths.DAC(index_sort)';

%[CTD_points] = calculate_CTDPoints(Floats_list,'CONFIG_ParkPressure_dbar',0);
%Floats_list.WMO = Floats_paths.WMO;

%% get vertical km and float age

% vertical km
[Analysis.vertical_km, Analysis.vertical_km_mean,Analysis.float_age, Analysis.last_cycle] = get_verticalkm_age_multiprof(Floats_list, dac_dir);

% format
index_nan = isnan(Analysis.float_age);
Analysis.vertical_km(index_nan) = [];
Analysis.vertical_km_mean(index_nan) = [];
Analysis.float_age(index_nan) = [];
Analysis.last_cycle(index_nan) = [];
Floats_list.STATUS(index_nan,:) = [];
Floats_list.DEPLOYMENTDATE(index_nan,:) = [];

idx_2019 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2019');
idx_2018 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2018');
idx_2017 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2017');
idx_2016 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2016');
idx_2015 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2015');
idx_2014 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2014');
idx_2013 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2013');
idx_2012 = contains(cellstr(Floats_list.DEPLOYMENTDATE), '2012');

Analysis_cycle_2018 = Analysis.last_cycle(idx_2018);
Analysis_cycle_2017 = Analysis.last_cycle(idx_2017);
Analysis_cycle_2016 = Analysis.last_cycle(idx_2016);

% only death floats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index_death = (contains(cellstr(Floats_list.STATUS), 'INACTIVE') | contains(cellstr(Floats_list.STATUS), 'CLOSED'));
% Analysis.vertical_km(~index_death) = [];
% Analysis.vertical_km_mean(~index_death) = [];
% Analysis.float_age(~index_death) = [];
% Analysis.last_cycle(~index_death) = [];
% Floats_list.STATUS(~index_death,:) = [];

n_floats = length(Analysis.vertical_km);

%% get survival rate

% get biggest cycles number
max_cycle = max(Analysis.last_cycle);
max_verticalkm = max(Analysis.vertical_km);
max_float_age = max(Analysis.float_age);
 

    death_floats = sum(contains(cellstr(Floats_list.STATUS), 'INACTIVE')) + sum(contains(cellstr(Floats_list.STATUS), 'CLOSED'));
    alive_index = ~(contains(cellstr(Floats_list.STATUS), 'INACTIVE') | contains(cellstr(Floats_list.STATUS), 'CLOSED'));
    death_index = (contains(cellstr(Floats_list.STATUS), 'INACTIVE') | contains(cellstr(Floats_list.STATUS), 'CLOSED'));
    
    % loop all cycles and counting
    sample_s = NaN(1,max_cycle);
    floats_survived = NaN(1,max_cycle);
    floats_not_survived = NaN(1,max_cycle);
    plot_data_cycle = NaN(1,max_cycle);
    plot_data_cycle_DEAD = NaN(1,max_cycle);
    plot_data_cycle_ALIVE = NaN(1,max_cycle);
    for icycle = 1:max_cycle
        sample_s(icycle) =  sum(Analysis.last_cycle(~alive_index) <= icycle) + sum(Analysis.last_cycle >= icycle);
        floats_survived(icycle) = (sum(Analysis.last_cycle >= icycle)/n_floats)*100;
        floats_not_survived(icycle) = (sum(Analysis.last_cycle <= icycle)/n_floats)*100;
        plot_data_cycle(icycle) = sum(Analysis.last_cycle >= icycle)/sample_s(icycle)*100;
        plot_data_cycle_formula2(icycle) = sum(Analysis.last_cycle >= icycle)/n_floats*100; 
        plot_data_cycle_DEAD(icycle) = sum(Analysis.last_cycle(death_index) >= icycle)/ sample_s(icycle)*100;
        plot_data_cycle_ALIVE(icycle) = sum(Analysis.last_cycle(alive_index) >= icycle)/ sample_s(icycle)*100;
    end

    
    % loop vertical km
    verticalkm_intervals = 0:max_verticalkm;
    sample_s =  NaN(1,length(verticalkm_intervals));
    plot_data_verticalkm = NaN(1,length(verticalkm_intervals));
    for ivkm = 1:length(verticalkm_intervals)
        sample_s(icycle) =  sum(Analysis.vertical_km(~alive_index) <= verticalkm_intervals(ivkm)) + ... 
                                   sum(Analysis.vertical_km >= verticalkm_intervals(ivkm));
        plot_data_verticalkm(ivkm) = sum(Analysis.vertical_km > verticalkm_intervals(ivkm))/sample_s(icycle)*100;
    end
    
    % loop float age
    floatage_intervals = 0:0.1:max_float_age;
    sample_s =  NaN(1,length(floatage_intervals));
    plot_data_floatage = NaN(1,length(floatage_intervals));
    for iage = 1:length(floatage_intervals)
        sample_s(iage) =  sum(Analysis.float_age(~alive_index) <= floatage_intervals(iage)) + ... 
                                   sum(Analysis.float_age >= floatage_intervals(iage));
        plot_data_floatage(iage) = sum(Analysis.float_age > floatage_intervals(iage))/sample_s(iage)*100;
    end
    
    


%% Import  extern variables for comparisons of life expectancies   
% data_AIC = readtable('/home1/datahome/co_arg/larduini/Lists/Survival Rates_2_data_1month_freq.txt');
% AIC_interval = data_AIC.Step;
% AIC_surv_rate = data_AIC.All;
% 
% surv_rate_vkm_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/surv_rate_vkm_ea_arvor_I_openocean.mat')
% data_open_surv_rate = surv_rate_vkm_openocean.plot_data_verticalkm;
% intervals_vkm_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/intervals_vkm_ea_arvor_I_openocean.mat')
% data_open_intervals = intervals_vkm_openocean.verticalkm_intervals;
% 
% surv_rate_vkm_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_vkm_marginalseas_surv_rate.mat')
% data_marginal_surv_rate = surv_rate_vkm_marginal.plot_data_verticalkm;
% intervals_vkm_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_vkm_marginalseas_intervals.mat')
% data_marginal_intervals = intervals_vkm_marginal.verticalkm_intervals;
% 
% surv_rate_vkm_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/Intern_array_surv_rate.mat')
% data_intern_surv_rate = surv_rate_vkm_intern.plot_data_verticalkm;
% intervals_vkm_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/Intern_array_intervals.mat')
% data_intern_intervals = intervals_vkm_intern.verticalkm_intervals;

%% Output folder

working_date = datestr(now,'yyyymmdd'); 

% check if output folder exits
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

%% Save variables
save([output_folder '/floatage_intervals_intern_open.mat'], 'floatage_intervals')
save([output_folder '/verticalkm_intervals_intern_open.mat'], 'verticalkm_intervals')
save([output_folder '/plot_data_cycle_intern_open.mat'], 'plot_data_cycle')
save([output_folder '/plot_data_floatage_intern_open.mat'], 'plot_data_floatage')
save([output_folder '/plot_data_verticalkm_intern_open.mat'], 'plot_data_verticalkm')


%% Analysis
% sample mean cycles
cycles_mean = mean(Analysis.last_cycle)
% sample mean age
age_mean = mean(Analysis.float_age)
% sample mean vertical km
vkm_mean = mean(Analysis.vertical_km)
% percentage of dead floats
per_dead = death_floats/n_floats*100

%% Plotting
disp(' ')
disp('Plotting ...')

% number of cycles
figure(1)
title_label = '2 different survival rate computations for European Iridium floats'
% full screen
% set(gcf, 'Position', get(0, 'Screensize'));

plot(plot_data_cycle, 'b','LineWidth',2, 'DisplayName', 'Survival rate from formula (1)');
hold on
plot(plot_data_cycle_formula2,'r', 'LineWidth',2, 'DisplayName', 'Survival rate from formula (2)');
% plot(floats_survived,'r','LineWidth',2, 'DisplayName', 'Floats that survived to the n cycle');
% hold on
% plot(plot_data_cycle_DEAD, 'g', 'LineWidth',2, 'DisplayName', 'Dead survival rate computations (dead floats only)');
% hold on
% plot(plot_data_cycle_ALIVE, 'k', 'LineWidth',2, 'DisplayName', 'Alive survival rate computations (alive floats only)');

% format
lgd = legend('Location','NorthEast');
lgd.FontSize = 18;
title(title_label, 'Interpreter', 'none', 'FontSize',20)
xlim([0 400])
ylabel('Survival rate (%)','FontSize',14)
xlabel('Number of cycles','FontSize',14)
set(gcf,'color','w')

%legend({['all floats (' num2str(n_floats) ')'], ['dead floats (' num2str(n_floats_DEAD) ')']})
%legend({['from 2008 (' num2str(n_floats) ')'], ['from 2014 (' num2str(n_floats_2014) ')'], ['from 2018 (' num2str(n_floats_2018) ')']})


% fig_path = [output_folder '/' 'SurvivalRate_cycles_DEAD_' working_date '.png'];
% export_fig(fig_path)



% % vertical km
% figure(2)
% title_label2 = 'Survival rate per vertical distance traveled for the International and European arrays'
% % full screen
% % set(gcf, 'Position', get(0, 'Screensize'));
% 
% 
% plot(data_intern_intervals, data_intern_surv_rate,'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the International array');
% hold on
% plot(verticalkm_intervals, plot_data_verticalkm,'m', 'LineWidth',4, 'DisplayName', 'Survival rate for all European Arvor-I floats');
% hold on
% plot(data_open_intervals, data_open_surv_rate,'b', 'LineWidth',2, 'DisplayName', 'Survival rate for European, open Ocean Arvor-I sample');
% hold on
% plot(data_marginal_intervals, data_marginal_surv_rate,'r', 'LineWidth',2, 'DisplayName', 'Survival rate for European, marginal Arvor-I sample');
% 
% %plot(verticalkm_intervals_2014, plot_data_verticalkm_2014,'LineWidth',2);
% %plot(verticalkm_intervals_2018, plot_data_verticalkm_2018,'LineWidth',2);
% % format
% title(title_label2, 'Interpreter', 'none', 'FontSize',18)
% lgd = legend('Location','NorthEast')
% lgd.FontSize = 18
% xlim([0 1200])
% ylabel('Floats surviving/ floats deployed (%)','FontSize',14)
% xlabel('Vertical distance (km)','FontSize',14)
% set(gcf,'color','w')
% %legend({['all floats (' num2str(n_floats) ')'], ['dead floats (' num2str(n_floats_DEAD) ')']})
% %legend({['from 2008 (' num2str(n_floats) ')'], ['from 2014 (' num2str(n_floats_2014) ')'], ['from 2018 (' num2str(n_floats_2018) ')']})
% 
% % fig_path = [output_folder '/' 'SurvivalRate_verticalkm_DEAD_' working_date '.png'];
% % export_fig(fig_path)
% 
% 
% 
% % float age
% figure(3)
% title_label3 = 'Comparison between Survival rates computed from Matlab Script and the Ocean Ops tool'
% titla_label3 = [title_label3  newline 'for European Arvor-I floats deployed in open Ocean']
% % full screen
% % set(gcf, 'Position', get(0, 'Screensize'));
% 
% plot(floatage_intervals, plot_data_floatage,'r', 'LineWidth',2, 'DisplayName', 'Survival rate computed from MatlabScript');
% hold on
% plot(AIC_interval/12, AIC_surv_rate, 'b', 'LineWidth',2,'DisplayName', 'Survival rate computed from the OceanOps tool');
% %plot(floatage_intervals_2014, plot_data_floatage_2014,'LineWidth',2);
% %plot(floatage_intervals_2018, plot_data_floatage_2018,'LineWidth',2);
% 
% % format
% title(title_label3, 'Interpreter', 'none')
% legend('Location','East')
% ylabel('Floats surviving/ floats deployed (%)','FontSize',12)
% xlabel('Floats age (years)','FontSize',12)
% set(gcf,'color','w')
% %legend({['all floats (' num2str(n_floats) ')'], ['dead floats (' num2str(n_floats_DEAD) ')']})
% %legend({['from 2008 (' num2str(n_floats) ')'], ['from 2014 (' num2str(n_floats_2014) ')'], ['from 2018 (' num2str(n_floats_2018) ')']})
% 
% fig_path = [output_folder '/' 'SurvivalRate_floatage_DEAD_' working_date '.png'];
% export_fig(fig_path)
% 
% figure(4)
% title_label4 = 'Characteriztion of the European Arvor-Iridium fleet'
% % full screen
% % set(gcf, 'Position', get(0, 'Screensize'));
% title(title_label4, 'Interpreter', 'none')
% 
% subplot(3,1,1)
% histogram(Analysis.last_cycle)
% title('Number of cycles in sample')
% xlabel('Number of cycles')
% ylabel('Number of floats')
% 
% subplot(3,1,2)
% histogram(Analysis.vertical_km)
% title('Vertical km in sample')
% xlabel('Vertical distance (km)')
% ylabel('Number of floats')
% 
% subplot(3,1,3)
% histogram(Analysis.float_age)
% title('Float age in sample')
% xlabel('Float age (years)')
% ylabel('Number of floats')
% 
% fig_path = [output_folder '/' 'SurvivalRate_HIST_DEAD_' working_date '.png'];
% export_fig(fig_path)


%%
% n_floats_DEAD = n_floats;
% plot_data_cycle_DEAD = plot_data_cycle;
% plot_data_verticalkm_DEAD = plot_data_verticalkm;
% plot_data_floatage_DEAD = plot_data_floatage;
% floatage_intervals_DEAD = floatage_intervals;

% n_floats_2014 = n_floats;
% plot_data_cycle_2014 = plot_data_cycle;
% plot_data_verticalkm_2014 = plot_data_verticalkm;
% verticalkm_intervals_2014 = verticalkm_intervals;
% plot_data_floatage_2014 = plot_data_floatage;
% floatage_intervals_2014 = floatage_intervals;

% n_floats_2018 = n_floats;
% plot_data_cycle_2018 = plot_data_cycle;
% plot_data_verticalkm_2018 = plot_data_verticalkm;
% verticalkm_intervals_2018 = verticalkm_intervals;
% plot_data_floatage_2018 = plot_data_floatage;
% floatage_intervals_2018 = floatage_intervals;