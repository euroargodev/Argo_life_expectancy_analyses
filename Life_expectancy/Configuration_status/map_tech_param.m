% NOT FINISHED - IN DEVELOPMENT

%map_tech_param
% plot a scatter of a tech, traj or config parameter in a map (lon, lat) 
% using a threshold
%
% Input (below)
% - floats_list: csv (or .txt) file, with in header WMO, RT (dac in charge 
%        of real time processing) or DAC (dac in charge with coriolis 
%        database format). Exemples: MOCCA_rt_dm.csv and ARVOR_5years.csv. 
%        Path to this file should be modified below
% - delimiter: floats_list delimiter (ex : ',',';')
% - gdac_path: path to gdac
% - param_name: (cell variable) name of parameter to be plot. It should be 
%        a numerical non discrete parameter. 
% - where_file: (cell variable) where to find the variable in param_name. 
%        It can be 'tech', 'aux' , 'traj' or 'config'. 
% - threshold: chosen threshold. It can be 'mean', 'median' or a number
% - colormap_thres: if 'above', colormap is use for points above threshold.
%        If 'under', colormap is use for points under threshold.
% - subtitle: subtitle in figure
% - grounded: if 'yes', it plots grounded cycles, if 'no', it does not plot
%        grounded cycles. It only works for ARVOR floats (NOT WORKING AT THE MOMENT).
% - working_folder: folder where all outputs will be saved
% - QGIS_csv: 'yes' or 'no'. If 'yes', it creates a .csv file with data 
%        easy to use in QGIS.
%
% Output
% - Map saved in (working_folder)/Figures folder (in .fig and .png formats)
% - data_(paramname)_(threshold)_(date).txt file : contains a list of floats 
%        above or under threshold, cycles above or under threshold and 
%        minimun and maximun technical values for each float (saved in 
%        (working_folder)/Data folder)
% - qgis_(paramname)_(threshold)_(date).csv file : (if QGIS_csv = 'yes') 
%        .csv with data easy to use in qgis for reproducing the same figure
%        (saved in (working_folder)/Data folder)
% - Screen output saved in screenout_paramname_threshold_date.txt, including
%        not plotted floats and grounded floats (saved in (today_date)/Data folder)
%
% Auxiliary functions needed:
%    read_csv
%    get_floats_data_gdac
%    M_MAP matlab package
%    export_fig package
%    get_csv_QGIS
%
% NOTE 
% (1) Depending on the number of floats, the script will take some minutes
% (2) If there is only one float in input list, trajectory lines are plotted.
% (3) Be carreful with RAM memory. If you use so many floats the program can
% run out of memory.
% (4) For logical or discrete parameters use map_tech_param_LOGICAL
%
% Modified on 19/09/2019

% add paths (packages and auxiliary functions)
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/seawater
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/m_map
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/MyTools % all functions I developed
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/addaxis % more than 2 y axis in a figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/legendflex % more than 1 column in legend
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/setgetpos_V1.2
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/ezyfit/ezyfit % fit data to a curve
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/aux_functions

close all
clear variables
% 
% % add paths (packages and auxiliary functions)
% aux_functions_path = [pwd '/aux_functions'];
% addpath(genpath(aux_functions_path))


% screen output in to a file
diary screenoutput.txt

%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input paths
list_dir = '/home1/datahome/co_arg/larduini/Lists/Comparison_EA_Intern/EA_ARVOR_I_ALL.csv'
delimiter = ',';
gdac_path = '/home/ref-argo/gdac/'
% parameters
param_name = {'NUMBER_RepositionsDuringPark_COUNT'} % NUMBER_PumpActionsDuringPark_COUNT ; NUMBER_RepositionsDuringPark_COUNT ; TIME_IridiumGPSFix_seconds ; VOLTAGE_BatteryPumpStartProfile_volts; CONFIG_ParkPressure_dbar; CONFIG_PumpActionTimeBuoyancyAcquisitionForInAirMeasCycle_csec
where_file = {'tech'} % tech, aux, traj or config
threshold = 0 % 'mean', 'median' or number
colormap_thres = 'above' % 'above' or 'under' (defaut above)
subtitle = 'Floats deployed in the Baltic Sea';  %'ARVOR floats last 5 years'
grounded = 'yes' % 'yes' or 'no' (default yes)
% output paths
working_folder = '/home1/datahome/co_arg/larduini/MAP_technical_param/'
QGIS_csv = 'no' ; % 'yes or 'no'

dac_dir = gdac_path;
% mc_code = {'all'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: grounded for APEX
% TODO: only one format for the input float list file?
% TODO: colormap in legend
% TODO: MC code for traj variables
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

%% get floats data
% get data from gdac
disp('Getting data ...')

% variables names and where to find them
traj_variables = {'LONGITUDE','LATITUDE'};
variables_list = [traj_variables, param_name];
where_file = [{'traj', 'traj'}, where_file];
mc_code={'all','all'};

[Floats] = get_floats_data_gdac_v3_FINAL(Floats_list, variables_list, where_file, dac_dir, mc_code);

% number of floats and cycles
n_floats = size(Floats.WMO.data,1);
n_cycles = cell2mat(cellfun(@(x)length(x),Floats.LATITUDE.data,'un',0));

%memory_in_use = monitor_memory_whos
%whos

%% Plot in map

disp(' ')
disp('Plotting ...')

close all

figure(1)
% % full screen
% set(gcf, 'Position', get(0, 'Screensize'));

% map and tech parameter limits
Limits.lon_min = min(cat(1,Floats.LONGITUDE.data{:}));
Limits.lon_max = max(cat(1,Floats.LONGITUDE.data{:}));
Limits.lat_min = min(cat(1,Floats.LATITUDE.data{:}));
Limits.lat_max = max(cat(1,Floats.LATITUDE.data{:}));
Limits.tech.max = max(cat(1,Floats.(param_name{1}).data{:}));
Limits.tech.min = min(cat(1,Floats.(param_name{1}).data{:}));
Limits.tech.mean = mean(cat(1,Floats.(param_name{1}).data{:}),'omitnan');
Limits.tech.median = median(cat(1,Floats.(param_name{1}).data{:}),'omitnan');

% if strange values in latitude
if Limits.lat_max > 90
    Limits.lat_max = 90;
end

% background map
m_proj('mercator','long',[Limits.lon_min-1  Limits.lon_max+1],'lat',[Limits.lat_min-1  Limits.lat_max+1]);
[CS,CH]=m_etopo2('contourf',[-10000:250:250 200 150 100 50 25 0],'edgecolor','none');
colormap(m_colmap('blues'));
m_coast('patch',[.7 .6 .4],'edgecolor',[.7 .6 .4]);
m_grid;

% m_elev('contour',[-1000 -1000],'color',[0.6 0.6 0.6]);
% hold on
% m_elev('contour',[-2000 -2000],'color',[0.7 0.7 0.7]);
% m_gshhs_i('patch',[0.9290 0.6940 0.1250]);
% m_grid('xtick',3,'ytick',3,'fontsize',7,'color',[0.5 0.5 0.5]);


hold on


%%%%%% tech param plot %%%%%

% data in blue regarding threshold
if contains(colormap_thres,'under')
    colorblue = 'above';
else
    colorblue = 'under';
end

% initialisations
Counters.floats.(colorblue) = cell(0);
Counters.cycles.(colorblue) = cell(0);
Counters.max_tech_value.(colorblue) = cell(0);
Counters.min_tech_value.(colorblue) = cell(0);
Counters.plotted_obs = 0;

for ifloat =1:n_floats % plot in blue loop 
    disp(ifloat)
    if contains(colormap_thres,'under')
         % threshold (big values in blue)
        if ischar(threshold)
            t_blue=double(Floats.(param_name{1}).data{ifloat} > Limits.tech.(threshold));
        else
            t_blue=double(Floats.(param_name{1}).data{ifloat} > threshold);
        end
    else
        % threshold (small values in blue)
        if ischar(threshold)
            t_blue=double(Floats.(param_name{1}).data{ifloat} <= Limits.tech.(threshold));
        else
            t_blue=double(Floats.(param_name{1}).data{ifloat} <= threshold);
        end
    end
    
    % do not plot grounded cycles
    if contains(grounded,'no')
        t_ground = double(~contains(string(Floats.GROUNDED.data{ifloat}),{'Y','B'})); % = 0 if gounded so cycle will not be plotted
        %t_ground = double(~contains(string(TrajData.GROUNDED.data{ifloat}),{'U'}));
        if ~all(t_ground == 1) % if grounded
            t_blue = t_blue.*t_ground; % counter in second loop
        end
    end
    t_blue(t_blue==0)= NaN;
    
    if n_floats == 1 % if only one float plot trajectory line
        m_plot(Floats.LONGITUDE.data{ifloat},Floats.LATITUDE.data{ifloat})
    end
     
    % plot in blue
    tech_blue = Floats.(param_name{1}).data{ifloat}.*t_blue;
    h1 = m_plot(Floats.LONGITUDE.data{ifloat}.*t_blue,Floats.LATITUDE.data{ifloat}.*t_blue,'Marker','o',...
        'MarkerEdgeColor','b','markersize',2);

    % list of floats and observations under theshold
    if ~all(isnan(t_blue)) % if plotted
       Counters.floats.(colorblue){end + 1} = Floats_list.WMO(ifloat,:);
       Counters.cycles.(colorblue){end + 1} = sum(t_blue,'omitnan');
       Counters.max_tech_value.(colorblue){end + 1} = max(tech_blue);
       Counters.min_tech_value.(colorblue){end + 1} = min(tech_blue);
    end
    
    % number of observations taken into account
    Counters.plotted_obs = Counters.plotted_obs + sum(t_blue,'omitnan');
    
end % under threshold loop


%initialisations
Counters.floats.(colormap_thres) = cell(0);
Counters.cycles.(colormap_thres) = cell(0);
Counters.max_tech_value.(colormap_thres) = cell(0);
Counters.min_tech_value.(colormap_thres) = cell(0);
Counters.grounded_floats = cell(0);
Counters.grounded_cycles = cell(0);

for ifloat =1:n_floats % colormap loop
            
    
    if contains(colormap_thres,'under')
        % threshold (small values in colormap)
        if ischar(threshold)
            t_color=double(Floats.(param_name{1}).data{ifloat} <= Limits.tech.(threshold));
        else
            t_color=double(Floats.(param_name{1}).data{ifloat} <= threshold);
        end
    else
        % threshold (big values in colormap)
        if ischar(threshold)
            t_color=double(Floats.(param_name{1}).data{ifloat} > Limits.tech.(threshold));
        else
            t_color=double(Floats.(param_name{1}).data{ifloat} > threshold);
        end
    end
    
    % do not plot grounded cycles
    if contains(grounded,'no')
        t_ground = double(~contains(string(Floats.GROUNDED.data{ifloat}),{'Y','B'})); % = 0 if gounded so cycle will not be plotted
        %t_ground = double(~contains(string(TrajData.GROUNDED.data{ifloat}),{'U'})); 
        if ~all(t_ground == 1) % if grounded
            t_color = t_color.*t_ground;
            Counters.grounded_floats{end+1} = Floats_list.WMO(ifloat,:);
            Counters.grounded_cycles{end+1} = sum(~t_ground);
            % screen output
            disp(' ')
            disp(Floats_list.WMO(ifloat,:))
            %TrajData.GROUNDED.data{ifloat}'
            disp('GROUNDED')
            fprintf('Number of cycles grounded: %i\n', Counters.grounded_cycles{end})
        end
    end
    
    t_color(t_color==0)= NaN;
        
        
    % plot above threshold
    tech_color = Floats.(param_name{1}).data{ifloat}.*t_color;
    % use colormap
    h2 = m_scatter(Floats.LONGITUDE.data{ifloat}.*t_color,Floats.LATITUDE.data{ifloat}.*t_color,...
            30,tech_color,'filled');
        
    
    % list of floats and observations above theshold
    if ~all(isnan(t_color)) % if plotted
       Counters.floats.(colormap_thres){end + 1} = Floats_list.WMO(ifloat,:);
       Counters.cycles.(colormap_thres){end + 1} = sum(t_color,'omitnan');
       Counters.max_tech_value.(colormap_thres){end + 1} = max(tech_color);
       Counters.min_tech_value.(colormap_thres){end + 1} = min(tech_color);
    end
    
    % number of observations taken into account
    Counters.plotted_obs = Counters.plotted_obs + sum(t_color,'omitnan');
         
end % above threshold loop

% figure format
if contains(grounded,'no')
    subtitle = [subtitle,' (without GROUNDED cycles)'];
end
title([sprintf('%s (updated %s)', param_name{1}, datestr(now(),'yyyy-mm-dd')) newline ...
     subtitle],'FontSize', 20, 'Interpreter', 'none')
% colorbar
if contains(colormap_thres,'above') 
    cmap = colormap(flipud(autumn));
else
    cmap = colormap(autumn);
end
h = colorbar;
%    units
units = strsplit(param_name{1},'_');
units = units(end);
ylabel(h, units,'FontSize', 16)
if ~isempty(min([Counters.min_tech_value.(colormap_thres){:}]))
    caxis([min([Counters.min_tech_value.(colormap_thres){:}]) max([Counters.max_tech_value.(colormap_thres){:}])])
end
% legend
if ischar(threshold)
    lh = legend([h1,h2],[colorblue ' threshold (<= ', num2str(Limits.tech.(threshold)),' [', threshold, '])'],...
        [colormap_thres ' threshold (> ', num2str(Limits.tech.(threshold)),' [', threshold, '])'],...
        'Location','southoutside');
else
    lh = legend([h1,h2],[colorblue ' threshold (<= ', num2str(threshold),')'],...
        [colormap_thres ' threshold (> ', num2str(threshold),')'],...
        'Location','southoutside');
end
% colormap in legend
% hold off
% %set(gcf,'units','pix')
% % get the position of the legend, and calculate the place for the colormaps:
% legpos = get(lh, 'Position');
% pos = [legpos(1)+0.008 legpos(2)-0.01 legpos(3)*0.1 legpos(4)*0.35];
% % Create a 'picture' of what you want to appear in the legend:
% legax = axes('Position',pos); % place the new picture above the legend
% imagesc(legax,1:max(size(autumn))) % Create the picture
% colormap(cmap) % appy custom colormap
% axis off % remove all axes details

% background color
set(gcf,'color','w');




%% Screen results and annotation in figure

%disp(' ')
fprintf('\nList of floats %s theshold : \n', colormap_thres)
for i = 1:length(Counters.floats.(colormap_thres))
   fprintf('  %s\n',Counters.floats.(colormap_thres){i})
end
per_floats = length(Counters.floats.(colormap_thres))/(n_floats)*100;
per_obs = sum([Counters.cycles.(colormap_thres){:}])/Counters.plotted_obs*100;
fprintf('\nNumber of floats in list: %i \n',size(Floats_list.WMO,1))
fprintf('Number of floats plotted: %i \n',n_floats)
fprintf('Number of floats %s threshold: %i (%f %%)\n', colormap_thres, length(Counters.floats.(colormap_thres)),per_floats)
fprintf('Number of observations plotted: %i \n', Counters.plotted_obs)
fprintf('Number of observations %s threshold: %i (%f %%)\n', colormap_thres, sum([Counters.cycles.(colormap_thres){:}]),per_obs)
% grounded cycles information
if contains(grounded,'no')
    per_floats_grounded = length(Counters.grounded_floats)/(n_floats)*100;
    per_obs_grounded = sum([Counters.grounded_cycles{:}])/(Counters.plotted_obs + sum([Counters.grounded_cycles{:}]))*100;
    fprintf('\nNumber of grounded floats: %i (%f %%)\n',length(Counters.grounded_floats),per_floats_grounded)
    fprintf('Number of grounded cycles: %i (%f %%)\n',sum([Counters.grounded_cycles{:}]),per_obs_grounded)
end

% annotation in figure
set(gcf,'units','pix')
% position relative to the legend
legpos = get(lh, 'Position');
annx = legpos(1) + legpos(3) + 0.02;
anny = legpos(2) - legpos(4) - 0.04;
annotation('textbox', [annx, anny, .1, .1], 'string', ['Number of floats in list: ',num2str(size(Floats_list.WMO,1)), newline,...
    'Number of floats plotted: ', num2str(n_floats), newline,...
    'Number of floats ', colormap_thres,' threshold: ', num2str(length(Counters.floats.(colormap_thres))),' (',num2str(round(per_floats,3)),'%)', newline,...
    'Number of observations plotted:  ', num2str(Counters.plotted_obs),newline,...
    'Number of observations ', colormap_thres,' threshold: ', num2str(sum([Counters.cycles.(colormap_thres){:}])),' (', num2str(round(per_obs,3)),'%)'],...
    'FitBoxToText','on','FontWeight','bold','FontSize',9,'LineStyle', 'none')


%% save figure

% get smaller param name for file name
file_param = strsplit(param_name{1},'_');
if isnumeric(threshold)
    threshold = num2str(threshold);
end

working_date = datestr(now,'yyyymmdd');
% check if output folder exits
if ~exist([working_folder '/Figures'], 'dir')
    mkdir([working_folder '/Figures'])
end
if ~exist([working_folder '/Data'], 'dir')
    mkdir([working_folder '/Data'])
end
% out files names
if contains(grounded,'no')
    fig_file = [working_folder '/Figures/map_' file_param{2} '_thr' threshold '_notGROUNDED_' working_date '.fig'];
    png_file = [working_folder '/Figures/map_' file_param{2} '_thr' threshold '_notGROUNDED_' working_date '.png'];
    txt_file = [working_folder '/Data/data_' file_param{2} '_thr' threshold '_notGROUNDED_' working_date '.txt'];
    QGIS_file = [working_folder '/Data/qgis_' file_param{2} '_thr' threshold '_notGROUNDED_' working_date '.csv'];
    screen_file = [working_folder '/Data/screenout_' file_param{2} '_thr' threshold '_notGROUNDED_' working_date '.txt'];
else
    fig_file = [working_folder '/Figures/map_' file_param{2} '_thr' threshold '_' working_date '.fig'];
    png_file = [working_folder '/Figures/map_' file_param{2} '_thr' threshold '_' working_date '.png'];
    txt_file = [working_folder '/Data/data_' file_param{2} '_thr' threshold '_' working_date '.txt'];
    QGIS_file = [working_folder '/Data/qgis_' file_param{2} '_thr' threshold '_' working_date '.csv'];
    screen_file = [working_folder '/Data/screenout_' file_param{2} '_thr' threshold '_' working_date '.txt'];
end

fprintf('\nSaving figure in %s ...\n',fig_file)
export_fig(fig_file)
fprintf('Saving figure in %s ...\n',png_file)
export_fig(png_file)
fprintf('Figure saved\n')

%% text file
fprintf('\nSaving results in %s ...\n',txt_file)
fid=fopen(txt_file,'w');

% title and information
fprintf(fid, '# \n');
fprintf(fid, '# Technical parameter: %s \n', param_name{1});
fprintf(fid, '# \n');
fprintf(fid, '# Floats list path : %s\n', floats_list);
fprintf(fid, '# %s\n', subtitle);
fprintf(fid, '# Update date : %s\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid, '# Threshold : %s\n', threshold);
fprintf(fid, '# Grounded cycles plotted : %s\n', grounded);
fprintf(fid, '# Colormap used for values %s threshold\n', colormap_thres);
fprintf(fid, '# \n');
fprintf(fid, '# \n');

% Table: Floats above or under threshold
fprintf(fid, '# Floats %s threshold\n', colormap_thres);
fprintf(fid, '# \n');
fprintf(fid, '# WMO: float reference\n');
fprintf(fid, '# n_cycles_%s_threshold: number of cycles %s threshold\n', colormap_thres, colormap_thres);
fprintf(fid, '# min_tech_param_value: minimun technical parameter value of cycles %s threshold\n',colormap_thres);
fprintf(fid, '# max_tech_param_value: maximun technical parameter value of cycles %s threshold\n',colormap_thres);
fprintf(fid, '# \n');
% header
header = ['WMO,' sprintf('n_cycles_%s_threshold,',colormap_thres) 'min_tech_param_value,' 'max_tech_param_value'];
fprintf(fid, '%s \n', header);
% data
table = [Counters.floats.(colormap_thres); Counters.cycles.(colormap_thres);...
         Counters.min_tech_value.(colormap_thres);Counters.max_tech_value.(colormap_thres)];
if ~isempty(table) % 0 plots under or above threshold 
    table = sortrows(table',2,'descend')'; % sort by number of cycles
    fprintf(fid, '%s,%d,%f,%f\n',table{:});
end

% sum up
fprintf(fid, '# \n');
fprintf(fid, '# Number of floats in list: %i \n',size(Floats_list.WMO,1));
fprintf(fid, '# Number of floats plotted: %i \n',n_floats);
fprintf(fid, '# Number of floats %s threshold: %i (%f %%)\n',colormap_thres, length(Counters.floats.(colormap_thres)),per_floats);
fprintf(fid, '# Number of observations plotted: %i \n',Counters.plotted_obs);
fprintf(fid, '# Number of observations %s threshold: %i (%f %%)\n',colormap_thres, sum([Counters.cycles.(colormap_thres){:}]),per_obs);

% if grounded
if contains(grounded,'no')
    fprintf(fid, '# \n');
    fprintf(fid, '# Number of grounded floats: %i (%f %%)\n',length(Counters.grounded_floats),per_floats_grounded);
    fprintf(fid, '# Number of grounded cycles: %i (%f %%)\n',sum([Counters.grounded_cycles{:}]),per_obs_grounded);
end
fclose(fid);
fprintf('Results saved\n')


%% save data for QGIS use

if strcmp(QGIS_csv,'yes')
    get_csv_QGIS(Floats, QGIS_file)
end


%% save screen output
% stop screen output file
diary off
movefile('screenoutput.txt',screen_file)
fprintf('\nScreen output saved in %s\n',screen_file)

