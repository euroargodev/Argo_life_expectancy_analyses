% CONFIG_fleet_status
% plots the status of a list of floats regarding a configuration parameter 
% and spliting the results depending on the country, deployment date and 
% float model
%


%
% Input (below)
% - list_dir: path to .csv or .txt file with floats list, including country
%       deployment date and model information. Format should be like file
%       deployed from JCOMMOPS website 
%       (exemple: ARVOR_PROVOR_from2008_noBGC_2010927.csv).
% - dac_dir: path to floats nc file (folder hierarchy should be the same as
%       in gdac, including index files)
% - config_param: configuration parameter to be analysed
% - units: configuration parameter units
% - output_folder: path to output folder
%
% Output
% - Figure config param changed?
% - Figure config param value for floats which did not change configuration
% - Figure config param value per cycle. 
% * All figures are saved in output_folder
%
% Auxiliary functions (from Toolbox)
%    read_csv
%    get_floats_filespath
%    get_floats_data_gdac
%    calculate_CyclePeriod
%    calculate_CTDPoints
%    format_data_for_plotting
%    get_matrix_barplot
%    suptitle
%
% NOTES:
% (1) meta index file is needed for looking for floats paths
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/30

addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/seawater
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/M_Map 
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/MyTools % all functions I developed
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/addaxis % more than 2 y axis in a figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/legendflex % more than 1 column in legend
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/setgetpos_V1.2
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/ezyfit/ezyfit % fit data to a curve
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox



clear variables

close all




% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list from JCOMMOPS with fields name = REF, STATUS, COUNTRY, DEPLOYMENT DATE, MODEL
%list_dir = '/home1/datahome/co_arg/agarciaj/life_expentancy/lists/ARVOR_Iridium_standardmission_from2008_globalocean.csv'
list_dir = '/home1/datahome/co_arg/larduini/Lists/Comparison_EA_Intern/ALL_EA_iridium_floats_since_2008.csv'
%list_dir = '/home1/datahome/co_arg/agarciaj/Other_scripts/life_expentancy/lists/MOCCA_list_jcommops.csv'
dac_dir = '/home/ref-argo/gdac/'
% configuration parameter one want to investigate and monitor its status
config_param = {'CONFIG_ProfilePressure_dbar'}
% The units variable is only used for confugration parameters exprssed with a time unit, such as Cycle_time,etc...
units = 'hours'
output_folder = '/home1/datahome/co_arg/larduini/Exports/Fleet_status/EA_Intern_comparison'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% get floats from list

% read list of floats
[Floats_list] = read_csv(list_dir,',');
Floats_list.WMO = cellstr(Floats_list.REF);


% when there is a not valid WMO in list
index_wmo = cellfun(@(x) length(x)~=7, Floats_list.WMO);
floats_list_fnames = fieldnames(Floats_list);
for inames = 1 : length(floats_list_fnames)
    Floats_list.(floats_list_fnames{inames})(index_wmo,:) = [];
end
Floats_list = rmfield(Floats_list,'REF');
floats_list_fnames = fieldnames(Floats_list);


% look for dac path
[Floats] = get_floats_filespath(dac_dir, cellstr(Floats_list.WMO));
index_dac = ~ismember(cellstr(Floats_list.WMO), Floats.WMO);
for inames = 1 : length(floats_list_fnames)
    Floats_list.(floats_list_fnames{inames})(index_dac,:) = [];
end


% sort by wmo before merge structs
[~, index_sort] = sort(Floats_list.WMO);
for inames = 1 : length(floats_list_fnames)
    Floats_list.(floats_list_fnames{inames}) = Floats_list.(floats_list_fnames{inames})(index_sort,:);
end
[~,index_sort] = sort(Floats.WMO);
Floats_list.DAC = Floats.DAC(index_sort)';
Floats_list.WMO = char(Floats_list.WMO);


%% get data from gdac
disp(' ')
disp('Getting data from gdac...')


% if CycleTime param
if contains(config_param,'CycleTime')
    calcul_param = config_param;
    config_param = {'CONFIG_CycleTime_days', 'CONFIG_CycleTime_hours', 'CONFIG_CycleTime_minutes'};
    where_file = {'config','config', 'config'};
elseif contains(config_param,'CTDPoints')
    % use also CONFIG_CycleTime_days
    calcul_param = config_param;
    config_param = {'CONFIG_CycleTime_hours', 'CONFIG_CycleTime_days',...
                    'CONFIG_ParkSamplingPeriod_hours', 'CONFIG_ProfilePressure_dbar',...
                    'CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar',...
                    'CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar',...
                    'CONFIG_ProfileSurfaceSlicesThickness_dbar',...
                    'CONFIG_ProfileIntermediateSlicesThickness_dbar',...
                    'CONFIG_ProfileBottomSlicesThickness_dbar'};
    where_file = {'config','config','config','config','config','config',...
                  'config','config','config'};
else
    calcul_param = config_param;
    where_file = {'config'};
end

[Floats] = get_floats_data_gdac_v3_FINAL(Floats_list, config_param, where_file, dac_dir);


%% Calculations
disp(' ')
disp('Calculations...')

if contains(calcul_param,'CycleTime')
    [Floats] = calculate_CyclePeriod(Floats,units);

elseif contains(calcul_param,'CTDPoints')
    [Floats] = calculate_CyclePeriod(Floats, 'hours');
    config_param(contains(config_param,'CONFIG_CycleTime_days')) = [];
    [Floats] = calculate_CTDPoints(Floats,config_param,0);
    
end


%% make format

disp(' ')
disp('Formatting...')

remove_cycles = [0,1];
[Floats,notused_floats] = format_data_for_plotting(Floats, remove_cycles);


n_floats = length(Floats.WMO.data);
n_cycles = cellfun(@(x) length(x), Floats.(calcul_param{1}).data);

% same floats for country, model and date
index_data = ismember(cellstr(Floats_list.WMO), notused_floats.WMO);
for inames = 1 : length(floats_list_fnames)
    Floats_list.(floats_list_fnames{inames})(index_data,:) = [];
end


%% Data analysis
disp(' ')
disp('Analysing data...')


Floats_list.DEPLOYMENTDATE = Floats_list.DEPLOYMENTDATE(:,1:4);

%%%%%%% Has the chosen config parameter changed durgin float lifetime ? %%%%%%%
CONFIG_changed = cellfun(@(x) length(unique(x(~isnan(x)))), Floats.(calcul_param{1}).data);


[barmatrix_changes_country, DChangesCountry] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.COUNTRY,'changes','change',0);
[barmatrix_changes_model, DChangesModel] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.MODEL,'changes','change',0);
[barmatrix_changes_date, DChangesDate] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.DEPLOYMENTDATE,'changes','change',0);


%%%%%%%%%%%% pie chart per float %%%%%%%%%%%%
% all different values of config variable
notchanged_data = Floats.(calcul_param{1}).data(CONFIG_changed == 1);
notchanged_country = Floats_list.COUNTRY(CONFIG_changed == 1,:);
notchanged_model = Floats_list.MODEL(CONFIG_changed == 1,:);
notchanged_ddate = Floats_list.DEPLOYMENTDATE(CONFIG_changed == 1,:);

[barmatrix_notchanged_country, DNotChangesCountry] = get_matrix_barplot(notchanged_data,notchanged_country, units, 'values', 0);
[barmatrix_notchanged_model, DNotChangesModel] = get_matrix_barplot(notchanged_data,notchanged_model, units, 'values', 0);
[barmatrix_notchanged_date, DNotChangesDate] = get_matrix_barplot(notchanged_data,notchanged_ddate, units, 'values', 0);

% sample size (only not changing floats)
n_floats_notchang = sum(sum(barmatrix_notchanged_country));
n_cycles_notchang = sum(cellfun(@(x)length(x),Floats.(calcul_param{1}).data(CONFIG_changed == 1)));


% %%%%%%%%%%%% pie chart per cycle %%%%%%%%%%%%
% % with all floats
[barmatrix_cycles_country, DCyclesCountry] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.COUNTRY, units, 'cycles', 1);
[barmatrix_cycles_model, DCyclesModel] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.MODEL, units, 'cycles', 1);
[barmatrix_cycles_date, DCyclesDate] = get_matrix_barplot(Floats.(calcul_param{1}).data,Floats_list.DEPLOYMENTDATE, units, 'cycles', 1);

n_cycles_after = sum(sum(barmatrix_cycles_country));

   
%% Output folder

working_date = datestr(now,'yyyymmdd'); 
short_name = strsplit(calcul_param{1},'_');
short_name = short_name{end -1};

% check if output folder exits
if ~exist([output_folder '/' calcul_param{1}], 'dir')
    mkdir([output_folder '/' calcul_param{1}])
end


%% make 3 plots: One if the configuration changed, one for the floats that did not change config and one with y-axis as number of cycles
disp(' ')
disp('Plotting...')


%% config parameter has changed?
figure(1)

super_title = [calcul_param{1} ' value changed?'];
h_title = suptitle(super_title);
h_title.Interpreter = 'none';
h_title.Units = 'normalized';
h_title.Position = [0.5, -0.0200, 0];


subplot(2,2,1)
conv_country_labels = str2double(DChangesCountry.labels); %% Converted labels from a cell array filled with character vectors to double digital array.
bar(conv_country_labels,DChangesCountry.floats/n_floats*100)
text(0.70,0.9,['Sample: ' num2str(n_floats) ' floats' newline ...
         '              ' num2str(sum(n_cycles)) ' cycles'], 'Units', 'normalized')
title('Per float' , 'Interpreter', 'none','FontSize',12)

% percentage labels
for ibar = 1: length(DChangesCountry.floats)
    txt = text(conv_country_labels(ibar), DChangesCountry.floats(ibar)/n_floats*100 + 2,[num2str(round(DChangesCountry.floats(ibar)/n_floats*100,1)) ' %'], ...
        'HorizontalAlignment', 'center','FontSize', 9);
    set(txt, 'Rotation', 0);
    
end

ylabel('% of floats','FontSize',10)
ylim([0,max(DChangesCountry.floats/n_floats*100)+5])
% xlim([min(conv_country_labels)-1 max(conv_country_labels)+1])
xlim([min(conv_country_labels)-1 15])
xlabel('Number of changes','FontSize',10)
%xaxis = (min(conv_country_labels):1:max(conv_country_labels))
xticks(conv_country_labels)
xtickangle(45)


subplot(2,2,2)
if length(DChangesCountry.labels) == 1
    bar([1,NaN],[barmatrix_changes_country'; nan(1,length(barmatrix_changes_country))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(DChangesCountry.labels))
    xlim([0.25 1.75]);
else
    conv_country_labels= str2double(DChangesCountry.labels);
    bar(conv_country_labels, barmatrix_changes_country', 'stacked')
end
legend(DChangesCountry.all_meta,'Location','eastoutside', 'FontSize', 8);
% xlim([min(conv_country_labels)-1 max(conv_country_labels)+1])
xlim([min(conv_country_labels)-1 15])
%xaxis = (min(conv_country_labels):1:max(conv_country_labels))
xticks(conv_country_labels)
xlabel('Number of changes','FontSize',10)
ylabel('Number of floats','FontSize',10)
title('Per country' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',9)

subplot(2,2,3)
if length(DChangesCountry.labels) == 1
    bar([1,NaN],[barmatrix_changes_model'; nan(1,length(barmatrix_changes_model))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(DChangesCountry.labels))
    xlim([0.25 1.75]);
else
    conv_model_labels= str2double(DChangesModel.labels);
    bar(conv_model_labels, barmatrix_changes_model', 'stacked')
end
legend(DChangesModel.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none');
% xlim([min(conv_model_labels)-1 max(conv_model_labels)+1])
xlim([min(conv_country_labels)-1 15])
xlabel('Number of changes', 'FontSize',10)
%xaxis = (min(conv_model_labels):1:max(conv_model_labels))
xticks(conv_model_labels)
ylabel('Number of floats','FontSize',10)
title('Per model' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',9)

subplot(2,2,4)
if length(DChangesCountry.labels) == 1
    bar([1,NaN],[barmatrix_changes_date'; nan(1,length(barmatrix_changes_date))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(DChangesCountry.labels))
    xlim([0.25 1.75]);
else
    conv_date_labels = str2double(DChangesDate.labels);
    bar(conv_date_labels, barmatrix_changes_date', 'stacked')
end
legend(DChangesDate.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none');
% xlim([min(conv_date_labels)-1 max(conv_date_labels)+1])
xlim([min(conv_country_labels)-1 15])
%xaxis = (min(conv_date_labels):1:max(conv_date_labels))
xticks(conv_date_labels)
xlabel('Number of changes', 'FontSize',10)
ylabel('number of floats','FontSize',10)
title('Per deployment year' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'jet')
xtickangle(45)
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'fontsize',9)

%format
%full screen
set(gcf, 'Position', get(0, 'Screensize'));

% figure name
set(gcf,'Name',[calcul_param{1} ' changes'])
% background color
set(gcf,'color','w');

% Export of the figure
fprintf('Exporting figure as png, comment if not needed ... \n')

fig1_path = [output_folder '/' calcul_param{1} '/' short_name '_valuechanges_' working_date '.png'];
print('-f1',fig1_path,'-dpng','-r300');

%% per float, not changed floats 
figure(2)

super_title = [calcul_param{1} ' for not changed floats'];
h_title = suptitle(super_title);
h_title.Interpreter = 'none';
h_title.Units = 'normalized';
h_title.Position = [0.5, -0.0200, 0];


subplot(2,2,1)
x_labels = categorical(DNotChangesCountry.labels);
x_labels = reordercats(x_labels, DNotChangesCountry.labels);
bar(x_labels, DNotChangesCountry.floats/n_floats_notchang*100, 0.5)
ylabel('% of floats','FontSize',10)
text(0.05,0.9,['Sample: ' num2str(n_floats_notchang) ' floats' newline ...
         '              ' num2str(sum(n_cycles_notchang)) ' cycles'], 'Units', 'normalized')
% text(0.1,0.9,['Sample: ' num2str(n_floats_notchang) ' floats' newline ...
%           '              ' num2str(sum(n_cycles_notchang)) ' cycles'], 'Units', 'normalized', 'Color', 'w')
% text(0.64,0.9,['Sample: ' num2str(n_floats_notchang) ' floats' newline ...
%          '              ' num2str(sum(n_cycles_notchang)) ' cycles'], 'Units', 'normalized')
title('Per float' , 'Interpreter', 'none','FontSize',12)
       
% percentage labels
for ibar = 1: length(DNotChangesCountry.floats)
    txt = text(ibar, DNotChangesCountry.floats(ibar)/n_floats_notchang*100 + 2.5,[num2str(round(DNotChangesCountry.floats(ibar)/n_floats_notchang*100,1)) ' %'], ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    set(txt, 'Rotation', 0);
end
ylim([0,max(DNotChangesCountry.floats/n_floats_notchang*100)+5])
%xlim([0.25 1.75])
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)

subplot(2,2,2)
if length(x_labels) == 1
    bar([1,NaN],[barmatrix_notchanged_country'; nan(1,length(barmatrix_notchanged_country))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(x_labels))
    xlim([0.25 1.75]);
else
    bar(x_labels, barmatrix_notchanged_country', 0.5, 'stacked')
end
legend(DNotChangesCountry.all_meta,'Location','eastoutside', 'FontSize', 8)
ylabel('number of floats','FontSize',10)
title('Per country' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)

subplot(2,2,3)
if length(x_labels) == 1
    bar([1,NaN],[barmatrix_notchanged_model'; nan(1,length(barmatrix_notchanged_model))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(x_labels))
    xlim([0.25 1.75]);
else
    bar(x_labels, barmatrix_notchanged_model', 0.5,  'stacked')
end
legend(DNotChangesModel.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none')
ylabel('number of floats','FontSize',10)
title('Per model' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)


subplot(2,2,4)
if length(x_labels) == 1
    bar([1,NaN],[barmatrix_notchanged_date'; nan(1,length(barmatrix_notchanged_date))], 'stacked')
    set(gca,'xtick',1)
    set(gca,'xticklabel',char(x_labels))
    xlim([0.25 1.75]);
else
    bar(x_labels, barmatrix_notchanged_date', 0.5, 'stacked')
end
legend(DNotChangesDate.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none')
ylabel('number of floats','FontSize',10)
title('Per deployment year' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'jet')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)

%format
% full screen
% set(gcf, 'Position', get(0, 'Screensize'));
% figure name
set(gcf,'Name',[calcul_param{1} ' per float'])
% background color
set(gcf,'color','w');

% Export of the figure
fprintf('Exporting figure as png, comment if not needed ... \n')

fig2_path = [output_folder '/' calcul_param{1} '/' short_name '_notchangedfloats_' working_date '.png'];
print('-f2',fig2_path,'-dpng','-r300');

%% per cycle 

figure(3)

super_title = [calcul_param{1} ' per cycle'];
h_title = suptitle(super_title);
h_title.Interpreter = 'none';
h_title.Units = 'normalized';
h_title.Position = [0.5, -0.0200, 0];


subplot(2,2,1)
x_labels = categorical(DCyclesCountry.labels);
x_labels = reordercats(x_labels, DCyclesCountry.labels);
bar(x_labels, DCyclesCountry.cycles/n_cycles_after*100)
ylabel('% of cycles','FontSize',16)
text(0.05,0.9,['Sample: ' num2str(n_floats) ' floats' newline ...
           '              ' num2str(sum(n_cycles_after)) ' cycles'], 'Units', 'normalized')
% text(0.64,0.9,['Sample: ' num2str(n_floats) ' floats' newline ...
%            '              ' num2str(sum(n_cycles_after)) ' cycles'], 'Units', 'normalized')
title('Per cycle' , 'Interpreter', 'none')
       
% percentage labels
for ibar = 1: length(DCyclesCountry.cycles)
    txt = text(ibar, DCyclesCountry.cycles(ibar)/n_cycles_after*100 + 2,[num2str(round(DCyclesCountry.cycles(ibar)/n_cycles_after*100,1)) '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
    set(txt, 'Rotation', 0);
end
ylim([0,max(DCyclesCountry.cycles/n_cycles_after*100)+5])
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize', 9)

subplot(2,2,2)
bar(x_labels, barmatrix_cycles_country', 'stacked')
legend(DCyclesCountry.all_meta,'Location','eastoutside', 'FontSize', 8)
ylabel('Number of cycles','FontSize',10)
title('Per country' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize', 9)

subplot(2,2,3)
bar(x_labels, barmatrix_cycles_model', 'stacked')
legend(DCyclesModel.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none')
ylabel('Number of cycles','FontSize',10)
title('Per model' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'colorcube')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)

subplot(2,2,4)
bar(x_labels, barmatrix_cycles_date', 'stacked')
legend(DCyclesDate.all_meta,'Location','eastoutside', 'FontSize', 8, 'Interpreter', 'none')
ylabel('Number of cycles','FontSize',10)
title('    Per deployment year' , 'Interpreter', 'none','FontSize',12)
colormap(gca,'jet')
xtickangle(45)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',9)

%format
% full screen
set(gcf, 'Position', get(0, 'Screensize'));
% figure name
set(gcf,'Name',[calcul_param{1} ' per cycle'])
% background color
set(gcf,'color','w');



% Export of the figure
fprintf('Exporting figure as png, comment if not needed ... \n')

fig3_path = [output_folder '/' calcul_param{1} '/' short_name '_percycle_' working_date '.png'];
print('-f3',fig3_path,'-dpng','-r300');




