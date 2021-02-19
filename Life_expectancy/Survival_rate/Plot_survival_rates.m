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
output_folder = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/Figs'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Import  extern variables for comparisons of life expectancies   

%%% OceanOps survival rate computation
data_AIC = readtable('/home1/datahome/co_arg/larduini/Lists/Survival Rates_2_data_1month_freq.txt');
AIC_interval = data_AIC.Step;
AIC_surv_rate = data_AIC.All;


%%%European all floats
%VKms
surv_rate_vkm_EA_all = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ALL/EA_all_vkm_surv_rate.mat');
data_EA_ALL_surv_rate_vkm = surv_rate_vkm_EA_all.plot_data_verticalkm;
intervals_vkm_EA_all = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ALL/EA_all_vkm_intervals.mat');
data_EA_ALL_intervals_vkm = intervals_vkm_EA_all.verticalkm_intervals;

%Float age
surv_rate_age_EA_all = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ALL/EA_all_floatage_surv_rate.mat');
data_EA_ALL_surv_rate_age = surv_rate_age_EA_all.plot_data_floatage;
intervals_age_EA_all = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ALL/EA_all_floatage_intervals.mat');
data_EA_ALL_intervals_age = intervals_age_EA_all.floatage_intervals;

%Cycles
surv_rate_cycle_EA_all = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ALL/EA_all_cycle_surv_rate.mat');
data_EA_ALL_surv_rate_cycle = surv_rate_cycle_EA_all.plot_data_cycle;


%%% European Arvor_I Open Ocean
%VKms
surv_rate_vkm_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/surv_rate_vkm_ea_arvor_I_openocean.mat');
data_EA_Arvor_I_surv_rate_vkm_open = surv_rate_vkm_openocean.plot_data_verticalkm;
intervals_vkm_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/intervals_vkm_ea_arvor_I_openocean.mat');
data_EA_Arvor_I_intervals_vkm_open = intervals_vkm_openocean.verticalkm_intervals;

%Float age
surv_rate_age_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/surv_rate_age_ea_arvor_I_openocean.mat');
data_EA_Arvor_I_surv_rate_age_open = surv_rate_age_openocean.plot_data_floatage;
intervals_age_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/intervals_age_ea_arvor_I_openocean.mat');
data_EA_Arvor_I_intervals_age_open = intervals_age_openocean.floatage_intervals;

%Cycles
surv_rate_cycle_openocean = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_openocean/surv_rate_cycle_ea_arvor_I_openocean.mat');
data_EA_Arvor_I_surv_rate_cycle_open = surv_rate_cycle_openocean.plot_data_cycle;

%%% European Arvor_I Marginal Seas
%VKms
surv_rate_vkm_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_vkm_marginalseas_surv_rate.mat');
data_EA_Arvor_I_surv_rate_vkm_marginal = surv_rate_vkm_marginal.plot_data_verticalkm;
intervals_vkm_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_vkm_marginalseas_intervals.mat');
data_EA_Arvor_I_intervals_vkm_marginal = intervals_vkm_marginal.verticalkm_intervals;

%Float age
surv_rate_age_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_age_marginalseas_surv_rate.mat');
data_EA_Arvor_I_surv_rate_age_marginal = surv_rate_age_marginal.plot_data_floatage;
intervals_age_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_age_marginalseas_intervals.mat');
data_EA_Arvor_I_intervals_age_marginal = intervals_age_marginal.floatage_intervals;

%Cycles
surv_rate_cycle_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_marginal_seas/EA_ARVOR_I_cycle_marginalseas_surv_rate.mat');
data_EA_Arvor_I_surv_rate_cycle_marginal = surv_rate_cycle_marginal.plot_data_cycle;


%%%European all Arvor-I
%VKms
surv_rate_vkm_EA_all_Arvor_I = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_ALL/EA_ARVOR_I_ALL_surv_rate_vkm.mat');
data_EA_ALL_Arvor_I_surv_rate_vkm = surv_rate_vkm_EA_all_Arvor_I.plot_data_verticalkm;
intervals_vkm_EA_all_Arvor_I = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_ALL/EA_ARVOR_I_ALL_intervals_vkm.mat');
data_EA_ALL_Arvor_I_intervals_vkm = intervals_vkm_EA_all_Arvor_I.verticalkm_intervals;

%Float age
surv_rate_age_EA_all_Arvor_I = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_ALL/EA_ARVOR_I_ALL_surv_rate_age.mat');
data_EA_ALL_Arvor_I_surv_rate_age = surv_rate_age_EA_all_Arvor_I.plot_data_floatage;
intervals_age_EA_all_Arvor_I = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_ALL/EA_ARVOR_I_ALL_intervals_age.mat');
data_EA_ALL_Arvor_I_intervals_age = intervals_age_EA_all_Arvor_I.floatage_intervals;

%Cycles
surv_rate_cycle_EA_all_Arvor_I = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/EA_ARVOR_I_ALL/EA_ARVOR_I_ALL_surv_rate_cycle.mat');
data_EA_ALL_Arvor_I_surv_rate_cycle = surv_rate_cycle_EA_all_Arvor_I.plot_data_cycle;

%%%International all
%VKms
surv_rate_vkm_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/intern_all_floats_surv_rate_vkm.mat');
data_intern_surv_rate_vkm = surv_rate_vkm_intern.plot_data_verticalkm;
intervals_vkm_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/intern_all_floats_vkm_intervals.mat');
data_intern_intervals_vkm = intervals_vkm_intern.verticalkm_intervals;

%Float age
surv_rate_age_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/intern_all_floats_surv_rate_age.mat');
data_intern_surv_rate_age = surv_rate_age_intern.plot_data_floatage;
intervals_age_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/intern_all_floats_age_intervals.mat');
data_intern_intervals_age = intervals_age_intern.floatage_intervals;

%Cycles
surv_rate_cycle_intern = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/intern_all_floats_surv_rate_cycle.mat');
data_intern_surv_rate_cycle = surv_rate_cycle_intern.plot_data_cycle;

%%%International open ocean
%VKms
surv_rate_vkm_intern_open = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_verticalkm_intern_open.mat');
data_intern_surv_rate_vkm_open = surv_rate_vkm_intern_open.plot_data_verticalkm;
intervals_vkm_intern_open = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/verticalkm_intervals_intern_open.mat');
data_intern_intervals_vkm_open = intervals_vkm_intern_open.verticalkm_intervals;

%Float age
surv_rate_age_intern_open = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_floatage_intern_open.mat');
data_intern_surv_rate_age_open = surv_rate_age_intern_open.plot_data_floatage;
intervals_age_intern_open = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/floatage_intervals_intern_open.mat');
data_intern_intervals_age_open = intervals_age_intern_open.floatage_intervals;

%Cycles
surv_rate_cycle_intern_open = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_cycle_intern_open.mat');
data_intern_surv_rate_cycle_open = surv_rate_cycle_intern_open.plot_data_cycle;

%%%International marginal seas
%VKms
surv_rate_vkm_intern_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_verticalkm_intern_marginal.mat');
data_intern_surv_rate_vkm_marginal = surv_rate_vkm_intern_marginal.plot_data_verticalkm;
intervals_vkm_intern_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/verticalkm_intervals_intern_marginal.mat');
data_intern_intervals_vkm_marginal = intervals_vkm_intern_marginal.verticalkm_intervals;

%Float age
surv_rate_age_intern_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_floatage_intern_marginal.mat');
data_intern_surv_rate_age_marginal = surv_rate_age_intern_marginal.plot_data_floatage;
intervals_age_intern_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/floatage_intervals_intern_marginal.mat');
data_intern_intervals_age_marginal = intervals_age_intern_marginal.floatage_intervals;

%Cycles
surv_rate_cycle_intern_marginal = matfile('/home1/datahome/co_arg/larduini/Exports/Survival_rates/EA_ARVOR_I/INTERN_ALL/plot_data_cycle_intern_marginal.mat');
data_intern_surv_rate_cycle_marginal = surv_rate_cycle_intern_marginal.plot_data_cycle;




%% Plotting
disp(' ')
disp('Plotting ...')

% vertical km
figure(1)
%full screen
set(gcf, 'Position', get(0, 'Screensize'));

title_label1 = 'Survival rate per vertical distance traveled for the European Arvor-I and the global European array'
% plot(data_intern_intervals_vkm, data_intern_surv_rate_vkm,'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the International array');
% hold on
% plot(data_open_intervals_vkm, data_open_surv_rate_vkm, 'r','LineWidth',4, 'DisplayName', 'Survival rate for the open Ocean European Iridium floats');
% hold on
% plot(data_intern_intervals_vkm, data_intern_surv_rate_vkm,'m', 'LineWidth',4, 'DisplayName', 'Survival rate for all International Iridium floats');
% hold on
% plot(data_intern_intervals_vkm_open, data_intern_surv_rate_vkm_open,'b', 'LineWidth',4, 'DisplayName', 'Survival rate for open Ocean International Iridium floats');
% hold on
% plot(data_intern_intervals_vkm_marginal, data_intern_surv_rate_vkm_marginal,'r', 'LineWidth',2, 'DisplayName', 'Survival rate for International Iridium floats deployed marginal Seas');
plot(data_EA_ALL_intervals_vkm, data_EA_ALL_surv_rate_vkm, 'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Iridium array');
hold on
plot(data_EA_ALL_Arvor_I_intervals_vkm, data_EA_ALL_Arvor_I_surv_rate_vkm, 'm', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Arvor-I floats');
hold on
plot(data_EA_Arvor_I_intervals_vkm_open, data_EA_Arvor_I_surv_rate_vkm_open, 'b', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in open Ocean');
hold on
plot(data_EA_Arvor_I_intervals_vkm_marginal, data_EA_Arvor_I_surv_rate_vkm_marginal, 'r', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in marginal Seas');


% format
title(title_label1, 'Interpreter', 'none', 'FontSize',18)
lgd = legend('Location','NorthEast');
lgd.FontSize = 18;
xlim([0 1200])
ylabel('Survival rate (%)','FontSize',14)
xlabel('Vertical distance (km)','FontSize',14)
set(gcf,'color','w')

fig_path = [output_folder '/' 'SurvivalRate_EA_Arvor_I_vkm.png'];
export_fig(fig_path)

% float age
figure(2)
%full screen
set(gcf, 'Position', get(0, 'Screensize'));

title_label2 = 'Survival rate per age reached for the global European Iridium array and the European Arvor-I floats'
% plot(data_intern_intervals_age, data_intern_surv_rate_age,'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the International array');
% hold on
% plot(data_open_intervals_age, data_open_surv_rate_age, 'r','LineWidth',4, 'DisplayName', 'Survival rate for the open Ocean European Iridium floats');
% hold on
% plot(data_intern_intervals_age, data_intern_surv_rate_age,'m', 'LineWidth',4, 'DisplayName', 'Survival rate for all International Iridium floats');
% hold on
% plot(data_intern_intervals_age_open, data_intern_surv_rate_age_open,'b', 'LineWidth',4, 'DisplayName', 'Survival rate for open Ocean International Iridium floats');
% hold on
% plot(data_intern_intervals_age_marginal, data_intern_surv_rate_age_marginal,'r', 'LineWidth',2, 'DisplayName', 'Survival rate for International Iridium floats deployed marginal Seas');
plot(data_EA_ALL_intervals_age, data_EA_ALL_surv_rate_age, 'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Iridium array');
hold on
plot(data_EA_ALL_Arvor_I_intervals_age, data_EA_ALL_Arvor_I_surv_rate_age, 'm', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Arvor-I floats');
hold on
plot(data_EA_Arvor_I_intervals_age_open, data_EA_Arvor_I_surv_rate_age_open, 'b', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in open Ocean');
hold on
plot(data_EA_Arvor_I_intervals_age_marginal, data_EA_Arvor_I_surv_rate_age_marginal, 'r', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in marginal Seas');

% format
title(title_label2, 'Interpreter', 'none', 'FontSize',18)
lgd = legend('Location','NorthEast');
lgd.FontSize = 18;
xlim([0 8])
ylabel('Survival rate (%)','FontSize',14)
xlabel('Age reached (years)','FontSize',14)
set(gcf,'color','w')

fig_path = [output_folder '/' 'SurvivalRate_EA_Arvor_I_age.png'];
export_fig(fig_path)

% number of cycles
figure(3)
%full screen
set(gcf, 'Position', get(0, 'Screensize'));


title_label3 = 'Survival rate per cycles reached for the global European Iridium array and the European Arvor-I floats'

% plot(data_intern_surv_rate_cycle,'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the International array');
% hold on
% plot(data_open_surv_rate_cycle, 'r','LineWidth',4, 'DisplayName', 'Survival rate for the European Iridium floats deployed the open Ocean');
% hold on
% plot(data_intern_surv_rate_cycle,'m', 'LineWidth',4, 'DisplayName', 'Survival rate for all International Iridium floats');
% hold on
% plot(data_intern_surv_rate_cycle_open,'b', 'LineWidth',4, 'DisplayName', 'Survival rate for open Ocean International Iridium floats');
% hold on
% plot(data_intern_surv_rate_cycle_marginal,'r', 'LineWidth',2, 'DisplayName', 'Survival rate for International Iridium floats deployed in marginal Seas');
plot(data_EA_ALL_surv_rate_cycle, 'g', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Iridium array');
hold on
plot(data_EA_ALL_Arvor_I_surv_rate_cycle, 'm', 'LineWidth',4, 'DisplayName', 'Survival rate for the European Arvor-I floats');
hold on
plot(data_EA_Arvor_I_surv_rate_cycle_open, 'b', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in open Ocean');
hold on
plot(data_EA_Arvor_I_surv_rate_cycle_marginal, 'r', 'LineWidth',2, 'DisplayName', 'Survival rate for the European Arvor-I floats deployed in marginal Seas');

% format
title(title_label3, 'Interpreter', 'none', 'FontSize',18)
lgd = legend('Location','NorthEast');
lgd.FontSize = 18;
xlim([0 500])
ylabel('Survival rate (%)','FontSize',14)
xlabel('Cycle number reached','FontSize',14)
set(gcf,'color','w')


fig_path = [output_folder '/' 'SurvivalRate_EA_Arvor_I_cycle.png'];
export_fig(fig_path)






