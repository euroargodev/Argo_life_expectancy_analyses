%% Script description
% This script permits to plot different survival rates computations stored as variables (.mat). It permits comparison between different floats samples, models, etc… 
% The survival rates of the floats samples, differentiation between models, etc. is done in an auxiliary script (to rename): /home1/datahome/co_arg/larduini/Andrea/Life_expectancy/plot_survival_rate.m
% The survival rates are computed according to different x-axis: number of cycles, float age and vertical distance traveled (in km)

% AUTHOR: Luca Arduini Plaisant, Euro-Argo ERIC
%         (luca.arduini.plaisant@euro-argo.eu)
%
% Modified on 2020/12/14
%%

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

clear variables
close all

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of variables to compare/plot
% var_name{1} = 'ALL_ALTO';
% var_name{2} = 'ALL_APEX';
% var_name{1} = 'APEX_BGC';
% var_name{1} = 'ALL_APEX_D';
% var_name{1} = 'ALL_ARVOR';
var_name{1} = 'ALL_ARVOR_A';
% var_name{6} = 'DEAD_ARVOR_I';
% var_name{2} = 'ALL_ARVOR_I_DO';
% var_name{1} = 'ALL_ARVOR_D';
% var_name{2} = 'DEAD_ARVOR_D';
% var_name{6} = 'ALL_HM2000';
% var_name{7} = 'ALL_NAVIS_A';
% var_name{8} = 'ALL_NAVIS_EBR';
% var_name{2} = 'NAVIS_A_BGC';
% var_name{3} = 'NAVIS_EBR_BGC';
% var_name{9} = 'ALL_NEMO';
% var_name{3} = 'ALL_NINJA_D';
% var_name{11} = 'ALL_NOVA';
% var_name{12} = 'ALL_PROVOR';
% var_name{13} = 'ALL_PROVOR_II';
% var_name{4} = 'ALL_PROVOR_III';
% var_name{5} = 'ALL_PROVOR_IV';
% var_name{6} = 'ALL_PROVOR_V';
% var_name{17} = 'ALL_S2A';
% var_name{18} = 'ALL_S2X';
% var_name{4} = 'ALL_SOLO_D';
% var_name{5} = 'ALL_SOLO_D_MRV';
% var_name{21} = 'ALL_SOLO_II';

% var_name{1} = 'ARCTIC';
% var_name{2} = 'ATLANTIC';
% var_name{3} = 'PACIFIC';
% var_name{4} = 'INDIAN';
% var_name{5} = 'SOUTHERN';
% var_name{5} = 'MEDITERRANEAN';
% var_name{6} = 'ARVOR_I_Open_Ocean';
% var_name{7} = 'BLACK';
% var_name{8} = 'BALTIC';
% var_name{9} = 'BANDA_SEA';
% var_name{10} = 'CARIBBEAN';
% var_name{11} = 'MEXICO_GULF';
% var_name{1} = 'APEX_D_JAPAN';
% var_name{2} = 'APEX_D_UK';
% var_name{3} = 'APEX_D_USA';

% var_name{1} = 'INTERNATIONAL_GLOBAL_IRIDIUM_ARRAY';
% var_name{2} = 'INTERNATIONAL_Marginal_Seas_IRIDIUM_ARRAY';
% var_name{3} = 'INTERNATIONAL_Open_Ocean_IRIDIUM_ARRAY';

% var_name{4} = 'EUROPEAN_GLOBAL_IRIDIUM_ARRAY';
% var_name{5} = 'EUROPEAN_Marginal_Seas_IRIDIUM_ARRAY';
% var_name{6} = 'EUROPEAN_Open_Ocean_IRIDIUM_ARRAY';


% filepath to where the .mat from the "compute_survival_rate.m" are stored
% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ALTO/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/APEX/';
% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/APEX_D/';
data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ARVOR_A/';
% data_folder{6} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ARVOR_I_DO/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ARVOR_I_DO/';
% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ARVOR_D/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/ARVOR_D/';
% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/APEX_D/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/APEX_D/';
% data_folder{3} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/APEX_D/';
% data_folder{6} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/HM2000/';
% data_folder{7} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/NAVIS_A/';
% data_folder{8} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/NAVIS_EBR/';
% data_folder{9} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/NEMO/';
% data_folder{3} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/NINJA_D/';
% data_folder{11} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/NOVA/';
% data_folder{12} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/PROVOR/';
% data_folder{13} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/PROVOR_II/';
% data_folder{14} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/PROVOR_III/';
% data_folder{15} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/PROVOR_IV/';
% data_folder{16} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/PROVOR_V/';
% data_folder{17} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/S2A/';
% data_folder{18} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/S2X/';
% data_folder{4}= '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/SOLO_D/';
% data_folder{5} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/SOLO_D_MRV/';
% data_folder{21} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/SOLO_II/';


% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_ARCTIC/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_ATLANTIC/';
% data_folder{3}= '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_PACIFIC/';
% data_folder{4} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_INDIAN/';
% data_folder{5} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_SOUTHERN/';
% data_folder{5} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_MEDITERRANEAN/';
% data_folder{6} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/Open_Ocean_ARVOR_I/';
% data_folder{7} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_BLACK/';
% data_folder{8} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_BALTIC/';
% data_folder{9} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_BANDA_SEA/';
% data_folder{10} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_CARIBBEAN/';
% data_folder{11} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_basins/ARVOR_I_BASINS/ARVOR_I_MEXICO_GULF/';

% data_folder{1} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';
% data_folder{2} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';
% data_folder{3} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';

% data_folder{4} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';
% data_folder{5} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';
% data_folder{6} = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_EA_Intern/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_folder = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/Comparison_models/Figs/';
fig_path = [output_folder '/' 'SurvivalRate_comparison_EA_Intern_floats_01.png'];
isexport = 0;

super_title = ['Survival rates comparison between EA and International arrays'];

sample_limit = 10 %Min number of floats necessary for a model to be ploted in the final output

%Loading the .mat files
for i= 1:length(var_name)
    
    if isempty(var_name{i})
        continue
    else
        % VKms
        surv_rate_vkm = matfile([data_folder{i} 'plot_data_verticalkm_' var_name{i} '.mat']);
        data_surv_rate_vkm{i} = surv_rate_vkm.plot_data_verticalkm;
        intervals_vkm = matfile([data_folder{i} 'verticalkm_intervals_' var_name{i} '.mat']);
        data_intervals_vkm{i} = intervals_vkm.verticalkm_intervals;

    %     Float age
        surv_rate_age = matfile([data_folder{i} 'plot_data_floatage_' var_name{i} '.mat']);
        data_surv_rate_age{i} = surv_rate_age.plot_data_floatage;
        intervals_age = matfile([data_folder{i} 'floatage_intervals_' var_name{i} '.mat']);
        data_intervals_age{i} = intervals_age.floatage_intervals;
    % 
    %     Cycles
        surv_rate_cycle = matfile([data_folder{i} 'plot_data_cycle_' var_name{i} '.mat']);
        data_surv_rate_cycle{i} = surv_rate_cycle.plot_data_cycle;
        n_floats{i} = surv_rate_cycle.n_floats;
    end

end
    

%% Plotting
disp(' ')
disp('Plotting ...')

% Creation of a mutli-color array
% Color = {'g','r','b','m','c', 'y','k', [0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],...
%         [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.8 0 0.8], [0.2 0.4 0.2], [0.4 0 0], [0.2 0.4 0.6],...
%         [0.5 0.6 0.8], [1 0.4 0.8], [0.4 0.4 0.6]};

Color = {[0.4940 0.1840 0.5560],[0 0.4470 0.7410],[0.8500 0.3250 0.0980], 'g','r','b',[0.9290 0.6940 0.1250],...
        [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.8 0 0.8], [0.2 0.4 0.2], [0.4 0 0], [0.2 0.4 0.6],...
        [0.5 0.6 0.8], [1 0.4 0.8], [0.4 0.4 0.6]};


% vertical km
figure(1)
%full screen
set(gcf, 'Position', get(0, 'Screensize'));
%super_title

annotation('textbox',    [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 16, 'FitBoxToText', 'on', 'Interpreter', 'none')

subplot(2,2,1)
title_label1 = 'Survival rate per cycles reached';

for i= 1:length(var_name)
    if n_floats{i} < sample_limit
        continue
    else
        plot(data_surv_rate_cycle{i}, 'LineWidth',2, 'color',Color{i}, 'DisplayName', ['Survival rate for ' var_name{i} ' floats (' num2str(n_floats{i}) ')']);
%         xline(400,'--',{'Theoretical lifetime', 'Marginal Seas'},'color','r')
%         xline(270,'--',{'Theoretical lifetime', 'Open Ocean'},'color','b')

    end
    hold on
    
end

% format
title(title_label1, 'Interpreter', 'none', 'FontSize',12)
% lgd = legend('Location','NorthEast','Interpreter', 'none');
% lgd.FontSize = 11;
xlim([0 450])
ylabel('Survival rate (%)','FontSize',10)
xlabel('Cycle number reached','FontSize',10)
set(gcf,'color','w')


subplot(2,2,2)
title_label2 = 'Survival rate per age reached';

for i= 1:length(var_name)
    if n_floats{i} < sample_limit
        continue
    else
        plot(data_intervals_age{i}, data_surv_rate_age{i}, 'LineWidth',2, 'color', Color{i}, 'DisplayName', ['Survival rate for ' var_name{i} ' floats (' num2str(n_floats{i}) ')']);
    end
    hold on
end

% format
title(title_label2, 'Interpreter', 'none', 'FontSize',12)
% lgd = legend('Location','NorthEast','Interpreter', 'none');
% lgd.FontSize = 11;
xlim([0 6])
ylabel('Survival rate (%)','FontSize',10)
xlabel('Age reached (years)','FontSize',10)
set(gcf,'color','w')


subplot(2,2,3:4)
title_label3 = 'Survival rate per vertical distance traveled';

for i= 1:length(var_name)
    if n_floats{i} < sample_limit
        continue
    else
        plot(data_intervals_vkm{i}, data_surv_rate_vkm{i}, 'LineWidth',2, 'color', Color{i},'DisplayName', ['Survival rate for ' var_name{i} ' floats (' num2str(n_floats{i}) ')']);
%         xline(1080,'LineWidth',2, 'color','r','LineStyle','--','DisplayName','coucou')
    end
    hold on
end

% format
title(title_label3, 'Interpreter', 'none', 'FontSize',12)
lgd = legend('Location','eastoutside','Interpreter', 'none');
lgd.FontSize = 11;
xlim([0 1600])
ylabel('Survival rate (%)','FontSize',10)
xlabel('Vertical distance (km)','FontSize',10)
set(gcf,'color','w')


if isexport == 1
    export_fig(fig_path)
end

