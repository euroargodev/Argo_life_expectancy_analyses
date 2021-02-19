%% Description of the script 

% This script aims to plot on a map the different types of groundings flags per cycle according to a WMO list furnished in the input.
% Groudings flags are stored in the traj files of a float. It can take either: 
% Y= Cycle DID ground; N= Cycle did NOT ground; U= Unknown (corresponds to floats that do not store this information. Apex floats before APF11 Version)
% You can either put a WMOs list as an input or just a single WMO number (by commenting the list part)
%
% AUTHOR: Luca Arduini Plaisant, Euro-Argo ERIC
%         (luca.arduini.plaisant@euro-argo.eu)
%
% Modified on 2020/12/14

%% Add additional functions
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/M_map
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/MyTools % all functions I developed
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/seawater
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/addaxis % more than 2 y axis in a figure
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/legendflex % more than 1 column in legend
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/setgetpos_V1.2
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/ezyfit/ezyfit % fit data to a curve
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/aux_functions

clear variables
close all
 

%% WMO input

% Raw export from JCOMMOPS with the following columns: COUNTRY,DEPLOYMENT DATE,MODEL,REF,STATUS

% List of floats 
list_dir = '/home1/datahome/co_arg/larduini/Lists/Marginal_seas/Marginalseas_EA_Iridium_floats_since_2008.csv'

% Baltic floats list
%list_dir = '/home1/datahome/co_arg/larduini/Lists/Export_from_jcommops_baltic_floats.csv'

% Black_sea floats list
% list_dir_black_sea = '/home1/datahome/co_arg/larduini/Lists/Marginal_seas/Export_from_jcommops_black_sea_floats.csv';

% Or single WMO number
% WMO_number = 6903704

dac_dir= '/home/ref-argo/gdac/';
dac_path = [dac_dir 'dac/'];
export_output_folder = '/home1/datahome/co_arg/larduini/Exports/Groundings/';
study_region = 'European Marginal Seas'; % Location of the sample, for the figures titles
study_region_file = 'European Marginal Seas'; %Location of the sample, with underscore in the middle to name file properly


% Single file path
%file_prof = '/home/ref-argo/gdac/dac/coriolis/7900466/7900466_prof.nc';
% file_traj = '/home/ref-argo/gdac/dac/coriolis/7900466/7900466_Rtraj.nc';

% Get information of the NetCdf file specified
%ncinfo(file_prof);
% ncinfo(file_traj)
% ncdisp(file_traj)
% 
% ncread(file_traj,'CYCLE_NUMBER_ACTUAL')

%% read list of floats

[Floats_list] = read_csv(list_dir,',');
Floats_list.WMO = cellstr(Floats_list.REF);

% [Floats_list_black_sea] = read_csv(list_dir_black_sea,',');
% Floats_list_black_sea.WMO = cellstr(Floats_list_black_sea.REF);

[Floats_list] = rmfield(Floats_list,'REF');
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


%% Determine the parameters wanted and get the data from the GDAC

var_params = {'LATITUDE', 'LONGITUDE', 'POSITION_QC', 'CYCLE_NUMBER', 'CYCLE_NUMBER_INDEX',  'CYCLE_NUMBER_ACTUAL','GROUNDED'};
where_file = {'traj', 'traj', 'traj', 'traj',  'traj', 'traj', 'traj'};
mc_code = {'all', 'all', 'all', 'all', 'all', 'all', 'all'};


[floats,WMO2delete_index] = get_floats_data_gdac_v3_FINAL(Floats_list, var_params, where_file, dac_dir, mc_code);
n_floats = length(floats.WMO.data);

% % WMO2delete_index is the index of the WMOs position in the array which has no traj file...
% 
% % In this loop, I want to update the index list to delete with the WMO of the black sea, as I am focusing here on the mediterranean basin (to comment if sea region ~= mediterranean)
% for j=1:n_floats
%     for jj=1:length(Floats_list_black_sea.WMO)
%         if strcmp(floats.WMO.data(j), Floats_list_black_sea.WMO(jj)) == 1
%             disp(j)
%             idx_black_sea(j) = j;
%         else
%             disp('no correspondance')
%         end
%     end
% end
% 
% idx_black_sea = idx_black_sea(idx_black_sea >0);
% 
% % Update "floats" struct array by deleting WMOs without traj file and WMOs corresponding to the black sea basin
% 
% WMO2delete_index = WMO2delete_index(~isnan(WMO2delete_index));
% index_list2delete = sort([WMO2delete_index idx_black_sea]);
% floats2 = floats;
% fields_list = fieldnames(floats);
% 
% for i= 1:length(fields_list)
%     disp(i)
%     field = fields_list(i);
%     floats2.(field{1}).data(index_list2delete) = [];
% end




% Update Floats list array by deleting WMOs without traj file

% Floats_list2 = Floats_list;
% fields_names = fieldnames(Floats_list);
% 
% for ii = 1:length(fields_names)
%     disp(ii)
%     name = fields_names(ii);
%     Floats_list2.(name{1})(index_list2delete,:) = [];
% end


floats2 = floats;
Floats_list2 = Floats_list
n_floats = length(floats2.WMO.data);




%% Looping through floats

fprintf('Starting calculations : \n')

for ifloats = 1:n_floats
    
    disp(ifloats)
    fprintf('Calculation progress: %.2f %%\n',(ifloats/max(n_floats))*100)
    
    
    % Définition des variables

    % On prend la totalité des values, pour tous les cycles, et tous les
    % flotteurs
    longitude = floats2.LONGITUDE.data(ifloats);
    % On transpose le cell array d'une ligne/n_floats columns à n_floats ligne/
    % 1 column
    longitude = longitude';
    % On convertit le cell array en matrice
    
    % Arbritary limit on latitude to limit to floats still in the mediterranean basin
    longitude = cell2mat(longitude);

    if any(longitude <= -5) == 1
        continue
    end
    
    
    latitude = floats2.LATITUDE.data(ifloats);
    latitude = latitude';
    latitude = cell2mat(latitude);
    
    % Arbritary limit on latitude to limit to floats still in the mediterranean basin
    if any(latitude >= 45) == 1
        continue
    end
    
    position_qc = floats2.POSITION_QC.data(ifloats);
    position_qc = position_qc';
    position_qc = cell2mat(position_qc);
    
    % Condition if position QC is flagged 4 (Bad data, with positions tests failed...)
    if any(position_qc == '4') == 1
        continue
    end
    
    cycle_nb = floats2.CYCLE_NUMBER.data(ifloats);
    cycle_nb = cycle_nb';
    cycle_nb = cell2mat(cycle_nb);

    if cellfun(@(x) isempty(x), floats2.CYCLE_NUMBER_INDEX.data(ifloats))
        cycle_nb_index = floats2.CYCLE_NUMBER_ACTUAL.data(ifloats);
    else    
        cycle_nb_index = floats2.CYCLE_NUMBER_INDEX.data(ifloats);
    end
    cycle_nb_index = cycle_nb_index';
    cycle_nb_index = cell2mat(cycle_nb_index);

    if max(cycle_nb_index) ~= max(cycle_nb)
        continue
    elseif max(cycle_nb) <= 0
        continue
    end
      
    grounded = floats2.GROUNDED.data(ifloats);
    grounded = grounded';
    grounded = cell2mat(grounded);

    % Arrangement en table des variables de meme dimension (cycle_nb, lat, lon)
    % de dimension N_measurements :
    table1 = table(cycle_nb, latitude, longitude, position_qc);
    
    % Suppression des cycles -1 et 0
    table1(table1.cycle_nb <= 0 ,:) = [];
    
    
    % Arrangement en table des variables de la dimension N_Cycle (grounded et
    % cyle_number_index)
    table4 = table(cycle_nb_index,grounded);
  
    % On supprime la ligne correspondant au cycle 0
    table4(table4.cycle_nb_index<=0 ,:) = [];

    
    % Création d'une condition si il n'y a pas de lat/lon pour certains
    % cycles (juste nans). Car une fois les nans enlevés, les dimensions
    % entre cycle_nb et cycle_nb_index seront différentes, rendant
    % impossible la concaténation des matrices
    
    empty_cycles = [];
    for i = 1:max(cycle_nb)
        %disp(i)
        data1 = table1(table1.cycle_nb == i,:);
        if isempty(data1(~isnan(data1.latitude),:))
            %disp('EMPTY CYCLE WITH NO LAT/LON')
            empty_cycles(i,:) = i;
        end
    end
    empty_cycles = empty_cycles(empty_cycles>0);

    for i=empty_cycles'
        table4(table4.cycle_nb_index == i,:) = [];
    end

    
    % Suppression des lignes avec des nans en lat/lon
    table2 = table1(~isnan(table1.latitude),:);
 
    % On supprime les lignes pour le même cycle number (dérive à la surface et
    % acquisition de plusieurs points GPS)
    %C = unique(table2);
    [C,index_table2] = unique(table2.cycle_nb);
    table3 = table2(index_table2,:);

   


    % On ajoute le WMO number à la table
    WMO_str = char(floats2.WMO.data(ifloats));
    WMO = cell(length(table4.grounded),1);
    WMO(:) = {(WMO_str)};

    % On ajoute le type de flotteur
    MODEL_str = char(Floats_list2.MODEL(ifloats,1:10));
    MODEL = cell(length(table4.grounded),1);
    MODEL(:) = {MODEL_str};

    table3 = [WMO MODEL table3];
    table3.Properties.VariableNames([1 2]) = {'WMO','MODEL'};
    
    % On merge les deux tables initialement de dimensions différentes:
    all_tables{ifloats} = [table3  table4];
   
    
end

%% Définition des tableaux finaux et des cycles grounded ou non, ou unknown :

table_final = vertcat(all_tables{:});
grounded_table = table_final(table_final.grounded == 'Y',:);
non_grounded_table = table_final(table_final.grounded == 'N',:);
unknown_grounded_table = table_final(table_final.grounded == 'U',:);

fprintf('End of calculations ! \n')
floats_number = length(unique(table_final.WMO));
fprintf('Number of remaining floats taken into account: %d \n', floats_number)

%% PLOTS
% map and tech parameter limits
Limits.lon_min = min(grounded_table.longitude);
Limits.lon_max = max(grounded_table.longitude);
Limits.lat_min = min(grounded_table.latitude);
Limits.lat_max = max(grounded_table.latitude);
% Limits.tech.max = max(cat(1,Floats.(param_name{1}).data{:}));
% Limits.tech.min = min(cat(1,Floats.(param_name{1}).data{:}));
% Limits.tech.mean = mean(cat(1,Floats.(param_name{1}).data{:}),'omitnan');
% Limits.tech.median = median(cat(1,Floats.(param_name{1}).data{:}),'omitnan');

% if strange values in latitude
if Limits.lat_max > 90
    Limits.lat_max = 90;
end



close all
%% 1ere figure: Statistical repartition of the grounding's flags

fprintf('Plotting First output \n')
figure(1)
fig_title = '_Pie_chart_grounding_flags_repartition';
% % full screen
% set(gcf, 'Position', get(0, 'Screensize'));

total_cycles = length(table_final.cycle_nb);
grounded_cycles = length(grounded_table.cycle_nb);
grounded_cycles_ratio = (grounded_cycles/total_cycles)*100;

non_grounded_cycles = length(non_grounded_table.cycle_nb);
non_grounded_cycles_ratio = (non_grounded_cycles/total_cycles)*100;

unknown_grounded_cycles = length(unknown_grounded_table.cycle_nb);
unknown_grounded_cycles_ratio = (unknown_grounded_cycles/total_cycles)*100;

ratio = [grounded_cycles_ratio non_grounded_cycles_ratio unknown_grounded_cycles_ratio];
% explode = [1 0];
labels = {'Grounded cycles', 'Non-grounded cycles', 'Unknown grounding'};

P = pie(ratio);
patch1 = P(1);
patch1.FaceColor = 'r';
text1 = P(2);
text1.Position = 0.5 * text1.Position;
text1.FontSize = 20;
patch3 = P(3);
patch3.FaceColor= 'g';
text2 = P(4);
text2.Position = 0.5 * text2.Position;
text2.FontSize = 20;
patch5 = P(5);
patch5.FaceColor= 'b';
text3 = P(6);
text3.Position = 0.5 * text3.Position;
text3.FontSize = 20;

str = ['Total cycle number : ' num2str(length(table_final.cycle_nb))];
text(1.25,0.2,str,'FontSize',14);
% annotation('textbox',[.2 .3],'String',str,'FitBoxToText','on');
title_str = ['Repartition of groundings flags for the EA marginal seas floats']; %%newline 'in the ' study_region];
title(title_str,'FontSize', 20);
legend(labels, 'Location','eastoutside', 'FontSize', 14);

% Export of the figure
% fprintf('Exporting figure as png, comment if not needed ... \n')
% print('-f1',[export_output_folder study_region_file fig_title],'-dpng','-r300');


%% 2eme figure : Geographical repartition of the groudning'Répartition géographique des grounding's flags

fprintf('Plotting Second output \n')
figure(2)
fig_title = '_Grounding_flags_geographic_repartition';
% full screen
% set(gcf, 'Position', get(0, 'Screensize'));

% x= [Limits.lon_min-2 Limits.lon_max+2 Limits.lon_max+2 Limits.lon_min-2]
% y = [Limits.lat_max+1 Limits.lat_max+1 Limits.lat_min-1 Limits.lat_min-1]
% m_patch(x,y,'red')
m_proj('mercator', 'long',[Limits.lon_min-2 Limits.lon_max+2],'lat',[Limits.lat_min-2 Limits.lat_max+1]);
[CS,CH]=m_etopo2('contourf',[-6000:250:250 200 150 100 50 25 0],'edgecolor','none');
% [bathy,lat,lon] = mygrid_sand2([Limits.lon_min-2 Limits.lon_max+2 Limits.lat_min-1 Limits.lat_max+1]);
% m_contourf(lon,lat,bathy);
colormap(m_colmap('blues'));
%m_contfbar(-500 0 'blue');
m_gshhs_i('patch',[0.9290 0.6940 0.1250]);
scatter3 = m_line(unknown_grounded_table.longitude,unknown_grounded_table.latitude,'Marker','o','color','b','LineStyle','none','markersize',4,'markerfacecolor','w');
scatter2 = m_line(non_grounded_table.longitude,non_grounded_table.latitude,'Marker','o','color','g','LineStyle','none','markersize',4,'markerfacecolor','w');
scatter1 = m_line(grounded_table.longitude,grounded_table.latitude,'Marker','o','color','r','LineStyle','none','markersize',4,'markerfacecolor','w');
m_grid;
title_str = ['Map of the Groundings flags per cycle for EA RISE floats' newline 'in the ' study_region];
title(title_str,'fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
labels = {'Grounded','Non grounded', 'Unknown'};
lgd = legend([scatter1 scatter2 scatter3], labels, 'Location','NorthWest');
lgd.FontSize = 15;

% Export of the figure
% fprintf('Exporting figure as png, comment if not needed ... \n')
% print('-f2',[export_output_folder study_region_file fig_title],'-dpng','-r300');




%% 3eme figure: Répartition géographique des grounded cycles avec une
% coloration des marqueur par WMO
fprintf('Plotting Third output \n')
figure(3)
fig_title = '_Grounded_cycles_repartition_per_WMOs';
% full screen
% set(gcf, 'Position', get(0, 'Screensize'));


grounded_WMOs = unique(grounded_table.WMO);
bins = length(grounded_WMOs);
ticks = [1:bins];
labels_WMO = cell2mat(grounded_WMOs);
A = str2double(grounded_table.WMO);
B = unique(A);
T = table(A,grounded_table.longitude, grounded_table.latitude);
T.Properties.VariableNames([2 3]) = {'longitude','latitude'};
CM = flipud(lines(bins));
CM = CM(randperm(length(CM)),:);
m_proj('mercator', 'long',[Limits.lon_min-2 Limits.lon_max+2],'lat',[Limits.lat_min-2 Limits.lat_max+1]);
[CS,CH]=m_etopo2('contourf',[-6000:250:250 200 150 100 50 25 0],'edgecolor','none');
% [lon,lat,bathy]=m_elev([Limits.lon_min-2 Limits.lon_max+2],'lat',[Limits.lat_min-1 Limits.lat_max+1]])
% [bathy,lat,lon] = mygrid_sand2([Limits.lon_min-2 Limits.lon_max+2 Limits.lat_min-1 Limits.lat_max+1]);
% m_contourf(lon,lat,bathy);
colormap(m_colmap('blues'));
m_gshhs_i('patch',[0.9290 0.6940 0.1250]);
scatter_x = zeros(bins,1);
for ii = 1:bins
    scatter_x(ii) = m_line(T.longitude(T.A==B(ii)), T.latitude(T.A==B(ii)),'Marker','o','color',CM(ii,:),'LineStyle','none','markersize',4,'markerfacecolor',CM(ii,:), 'MarkerEdgeColor','k');
end
m_grid;
title_str = ['Repartition of the grounded cycles according to their WMOs number' newline 'for EA RISE floats in the ' study_region];
title(title_str,'fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
% columnlegend(3,labels_WMO,'Location','southwest')
lgd = legend(scatter_x,labels_WMO,'Location','northwestoutside');
lgd.FontSize=15;

% Export of the figure
% fprintf('Exporting figure as png, comment if not needed ... \n')
% print('-f3',[export_output_folder study_region_file fig_title],'-dpng','-r300');


%%

fprintf('Plotting Fourth output \n')
figure(4)
fig_title = '_Unknown_Grounded_cycles_repartition_per_WMOs';
% full screen
set(gcf, 'Position', get(0, 'Screensize'));


unknown_grounded_WMOs = unique(unknown_grounded_table.WMO);
bins = length(unknown_grounded_WMOs);
ticks = [1:bins];
labels_WMO = cell2mat(unknown_grounded_WMOs);
A = str2double(unknown_grounded_table.WMO);
B = unique(A);
T = table(A,unknown_grounded_table.longitude, unknown_grounded_table.latitude);
T.Properties.VariableNames([2 3]) = {'longitude','latitude'};
CM = flipud(lines(bins));
CM = CM(randperm(length(CM)),:);
m_proj('mercator', 'long',[Limits.lon_min-2 Limits.lon_max+2],'lat',[Limits.lat_min-1 Limits.lat_max+1]);
[bathy,lat,lon] = mygrid_sand2([Limits.lon_min-2 Limits.lon_max+2 Limits.lat_min-1 Limits.lat_max+1]);
m_contourf(lon,lat,bathy);
colormap(m_colmap('blues',80));
m_gshhs_i('patch',[0.9290 0.6940 0.1250]);
scatter_x = zeros(bins,1);
for ii = 1:bins
    scatter_x(ii) = m_line(T.longitude(T.A==B(ii)), T.latitude(T.A==B(ii)),'Marker','o','color',CM(ii,:),'LineStyle','none','markersize',4,'markerfacecolor',CM(ii,:), 'MarkerEdgeColor','k');
end
m_grid;
title_str = ['Repartition of the unknow grounded cycles according to their WMOs number in the ' study_region];
title(title_str,'fontsize',10);
xlabel('Longitude');
ylabel('Latitude');
legend(scatter_x, labels_WMO,'Location','northwest');

% % Export of the figure
% fprintf('Exporting figure as png, comment if not needed ... \n')
% print('-f3',[export_output_folder study_region_file fig_title],'-dpng','-r300');


figure(5)
m_proj('mercator', 'long',[Limits.lon_min-2 Limits.lon_max+2],'lat',[Limits.lat_min-1 Limits.lat_max+1]);
[bathy,lat,lon] = mygrid_sand2([Limits.lon_min-2 Limits.lon_max+2 Limits.lat_min-1 Limits.lat_max+1]);
m_contourf(lon,lat,bathy);
colormap(m_colmap('blues',80));
% m_contfbar(-3000:250:0 ,'blue');
m_gshhs_i('patch',[0.9290 0.6940 0.1250]);
m_grid;


