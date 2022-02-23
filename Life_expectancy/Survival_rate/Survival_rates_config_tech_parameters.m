% CONFIG_survival_rate
% plots survival rates depending on config values for a given config 
% parameter (given in number of cycles, vertical_km and float age)
%
% Input
% - list_dir: path to floats list. Floats list should contains at least WMO
%       and STATUS (list from jcommops web site)
% - dac_dir: path to gdac
% - param: Name of configuration parameter to be analysed
% - sample_size_limit: only config parameter values used by a number of
%       floats > sample_size_limit will be plot
% - output_folder: path to output folder
%
% Output
% - Figures with survival rates depending on number of cycles, vertical km 
%       and float age 
%
% Auxiliary functions:
%    read_csv
%    get_floats_filespath
%    get_floats_data_gdac
%    calculate_CyclePeriod
%    calculate_CTDPoints
%    get_vertical_km_multiprof
%    format_data_for_plotting
%
% 
%
% Modified on 2021/01/05 by Luca Arduini


clear variables
close all


% add paths (packages and auxiliary functions)
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/MyTools
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list from JCOMMOPS with fields name = REF, STATUS, COUNTRY, DEPLOYMENT DATE, MODEL
% list_dir = '/home1/datahome/co_arg/larduini/Lists/Comparison_models/dead_Arvor_I.csv'
list_dir = '/home1/datahome/co_arg/larduini/Lists/Battery_dead_ARVOR_platform_type_floats.csv'
dac_dir = '/home/ref-argo/gdac/'

% Parameter to investigate in terms of impact on the survival rate opf a group of floats
param = {'NUMBER_PumpActionsDuringAscentToSurface_COUNT'} % If the number of CTD points is investigated, type 'CTDPoints' as a config parameter.
                     % If the groundings is investigated, type {'GROUNDED'} as a parameter, {'traj'} in the where file 
                     % and modify the mc_code on line 155 from {'all'} to {'all';'all'}   

% File abbreviation where the parameter to investigate is stored. Either: {'config'}, {'traj'}, {'tech'}, {'meta'}
where_file = {'tech'};

% Is the parameter to 6901181investigate a boolean or not ? The script will produce an output per values taken by the parameter if the param is not a boolean
% (like for example a config parameter such as: CONFIG_Park_pressure_dbar).
% If the param is a boolean, then, the scripts will plot survival rates in function of the percentages of cycles where the value of this parameter is 1
%(like for example : 'NUMBER_RepositionsDuringPark_COUNT').
is_boolean = 0; % Either 1 or 0 (Yes or No)

% Threshold for the interval parameter values: Used if the input parameter is not a boolean. The values taken by the parameter are going to be
% stacked in different intervals based on the Thresh_1 defined below in order to regroup multiple floats with parameter values in the same range

%Thresh_1 = 150; % based on the parameter values units

% Threshold based on the percentage of cycles affected by the parameter : Thresh_2 below will only be used if the input parameter is a boolean.
% This will allow to stack floats in groups where xx% of the total cycles have the input parameter as TRUE/1/ON,etc...

%Thresh_2 = 25; % percentage(%) of cycles

% Only consider sample size superior to xx floats(default is 10: the survival rate from a smaller sample would be less reliable).
% However, if all the outputs are wanted, put the sample_size_limit to 0.
sample_size_limit = 0; % floats
sample_size_limit2 = 0;

plot_group_model =1; % if 1: plots one output per model, with different parameter values curves. If 0, plots as many graphs as there are different param
                     % values and each curves represents a model.

% Output folder:
output_folder = '/home1/datahome/co_arg/larduini/Exports/Survival_rates/config_tech_param';

is_export = 0; % Either 0 or 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% get floats from list
[Floats_list] = read_csv(list_dir,',');

% number of floats
Floats_list.WMO = cellstr(Floats_list.REF);
% Floats_list.COUNTRY = cellstr(Floats_list.RT);

Floats_list_ini = Floats_list;

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
index_dac = ~ismember(cellstr(Floats_list.REF), Floats_paths.WMO);
for inames = 1 : length(Floats_list_field)
    Floats_list.(Floats_list_field{inames})(index_dac,:) = [];
end


% sort by wmo before merge structs
[Floats_list.WMO, index_sort] = sort(Floats_list.WMO);
for inames = 1 : length(Floats_list_field)
    Floats_list.(Floats_list_field{inames}) = Floats_list.(Floats_list_field{inames})(index_sort,:);
end
[Floats_paths.WMO, index_sort] = sort(Floats_paths.WMO);
Floats_list.DAC = Floats_paths.DAC(index_sort)';

Floats_list.WMO = char(Floats_paths.WMO); 


%% get data from gdac
disp(' ')
disp('Getting data from gdac...')


% if CycleTime param
if contains(param,'CycleTime')
    calcul_param = param;
    param = {'CONFIG_CycleTime_hours','CONFIG_CycleTime_days', 'CONFIG_CycleTime_minutes'};
    where_file = {'config','config','config'};
elseif contains(param,'CTDPoints')
    % use also CONFIG_CycleTime_days
    calcul_param = param;
    param = {'CONFIG_CycleTime_hours', 'CONFIG_CycleTime_days', 'CONFIG_CycleTime_minutes'...
                    'CONFIG_ParkSamplingPeriod_hours', 'CONFIG_ProfilePressure_dbar',...
                    'CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar',...
                    'CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar',...
                    'CONFIG_ProfileSurfaceSlicesThickness_dbar',...
                    'CONFIG_ProfileIntermediateSlicesThickness_dbar',...
                    'CONFIG_ProfileBottomSlicesThickness_dbar'};
    where_file = {'config','config','config','config','config','config',...
                  'config','config','config', 'config'};
else
    calcul_param = param;
end


% date from traj file
all_params = [param {'JULD'}];
where_file = [where_file {'traj'}];

if contains(param, 'GROUNDED')
	mc_code = {'all','all'};
else
	mc_code = {'all'};
end

% mc_code = {'all','all'};

% if strcmp(where_file{1}, 'aux') ==1
%     dac_dir = '/home/ref-argo/gdac/aux';
% end


[Floats,WMO2delete_index, notused_floats] = get_floats_data_gdac_v3_FINAL(Floats_list, all_params, where_file, dac_dir, mc_code);

%%% Warnings %%%
% If a float is     
idx_empty = cellfun(@isempty, Floats.(param{1}).data);
warning_string_1 = ['Number of floats whithout the input parameter requested: ' num2str(sum(idx_empty))];
disp(warning_string_1)
empty_models= Floats_list.MODEL(idx_empty,:);
empty_WMO= Floats_list.WMO(idx_empty,:);
empty_DAC = Floats_list.DAC(idx_empty,:);
disp('WMO and models of the floats without the input parameter requested :');
disp(strcat(empty_WMO,' : ',empty_models))


%% calculations
disp(' ')
disp('Calculations...')

% if contains(calcul_param,'CycleTime')
%     [Floats2] = calculate_CyclePeriod(Floats,'hours');

if contains(calcul_param,'CTDPoints')
%     [Floats2] = calculate_CyclePeriod(Floats,'days');
%     param(contains(param,'CONFIG_CycleTime_days')) = [];
    [Floats2] = calculate_CTDPoints(Floats,param,0);
    
    Floats = Floats2;
    calcul_param = {'CONFIG_CTDPoints_NUMBER'};
end


%% get vertical km and float age

% vertical km
[Analysis.vertical_km, Analysis.vertical_km_mean] = get_vertical_km_multiprof(Floats, dac_dir);


% float age
index = cellfun(@(x) isempty(x), Floats.JULD.data);
Floats.JULD.data(index) = {NaN};

Analysis.float_age = cellfun(@(x) max(x) - min(x), Floats.JULD.data)/365;

Floats = rmfield(Floats,'JULD');
Floats = rmfield(Floats,'REFERENCE_DATE_TIME');

if contains(calcul_param,'CTDPoints')
    Floats = rmfield(Floats, param);
end

%% make format
disp(' ')
disp('Formatting...')

remove_cycles = [0,1];

if strcmp(param{1},'GROUNDED') ==1
    for ifloat = 1:length(Floats.GROUNDED.data)
        for icycle = 1:length(Floats.GROUNDED.data{ifloat})
            if strcmp(Floats.GROUNDED.data{ifloat}(icycle), 'N') == 1
                Floats.GROUNDED.data{ifloat}(icycle) = 0;
            elseif strcmp(Floats.GROUNDED.data{ifloat}(icycle), 'U') == 1
                Floats.GROUNDED.data{ifloat}(icycle) = 0;
            elseif strcmp(Floats.GROUNDED.data{ifloat}(icycle), 'Y') == 1
                Floats.GROUNDED.data{ifloat}(icycle) = 1;
            end    
        end
    end
    
% elseif contains(calcul_param,'CTDPoints')
%     disp('no specific formatting')  
else
    [Floats, notused_floats] = format_data_for_plotting(Floats, remove_cycles);
end

    

n_floats = length(Floats.WMO.data);
n_cycles = cellfun(@(x) length(x), Floats.(calcul_param{1}).data);

% same floats for country, model and date
index_data = ismember(cellstr(Floats_list.WMO), notused_floats.WMO);
for inames = 1 : length(Floats_list_field)
    Floats_list.(Floats_list_field{inames})(index_data,:) = [];
end
Analysis.vertical_km(index_data) = [];
Analysis.vertical_km_mean(index_data) = [];
Analysis.float_age(index_data) = [];



n_floats = length(Floats.(calcul_param{1}).data);


%% get survival rate

if sum(contains(where_file, 'config')) ~= 0
    
    % first config value and number of cycles
    Analysis.first_configvalue = cellfun(@(x) x(1), Floats.(calcul_param{1}).data);
    Analysis.last_cycle = cellfun(@(x) x(end), Floats.(calcul_param{1}).cycle);

    % config changed or not
    Analysis.CONFIG_changed = NaN(n_floats,1);
    for ifloat = 1: n_floats
        Analysis.CONFIG_changed(ifloat) = length(unique(Floats.(calcul_param{1}).data{ifloat}));
    end
    n_not_changed = sum(Analysis.CONFIG_changed == 1);

    % use not changing floats
    notchan_configvalue = Analysis.first_configvalue(Analysis.CONFIG_changed == 1);
    notchan_last_cycle = Analysis.last_cycle(Analysis.CONFIG_changed == 1);
    notchan_verticalkm = Analysis.vertical_km(Analysis.CONFIG_changed == 1);
    notchan_float_age = Analysis.float_age(Analysis.CONFIG_changed == 1);
    notchan_WMO = Floats_list.WMO(Analysis.CONFIG_changed == 1,:);
    notchan_MODEL = Floats_list.MODEL(Analysis.CONFIG_changed == 1,:);
    notchan_status = Floats_list.STATUS(Analysis.CONFIG_changed == 1,:);

    % get all different config values
    unique_values = unique(notchan_configvalue);
    unique_values(isnan(unique_values)) = [];

    % get biggest cycles number
    max_cycle = max(notchan_last_cycle);
    max_verticalkm = max(notchan_verticalkm);
    max_float_age = max(notchan_float_age);

    plot_label = cell(1,length(unique_values));
    sample_size = NaN(1,length(unique_values));
end

if is_boolean == 0
    Analysis.last_cycle = cellfun(@(x) x(end), Floats.(calcul_param{1}).cycle);
    
    Analysis.min_value_param = cellfun(@min, Floats.(calcul_param{1}).data);
    Analysis.max_value_param = cellfun(@max, Floats.(calcul_param{1}).data);
    Analysis.mean_value_param = cellfun(@mean, Floats.(calcul_param{1}).data);
    
    Analysis.sum_value_param = cellfun(@sum, Floats.(calcul_param{1}).data);
    
    max_cycle =  Analysis.last_cycle;
    max_verticalkm = max(Analysis.vertical_km);
    max_float_age = max(Analysis.float_age);
end

if is_boolean == 1
    Analysis.last_cycle = cellfun(@(x) x(end), Floats.(calcul_param{1}).cycle);

    Analysis.sum_value_param = cellfun(@sum, Floats.(calcul_param{1}).data); % The total sum of the parameter value for all the cycles
    
    non_nul_cycle = cellfun(@find, Floats.(calcul_param{1}).data, 'UniformOutput', false); % numero of cycles having a parameter value different than 0
    Analysis.numb_cycles_affected = cellfun(@length, non_nul_cycle); % Number of cycles affected by the parameter (counting only one per cycle)
                                                            % This variable will be used for the percentage of cycles affected calculation
    
    percen = NaN(length(Analysis.last_cycle),1);
    for i = 1:length(Analysis.last_cycle)
        percen(i) = (Analysis.numb_cycles_affected(i) / Analysis.last_cycle(i))*100;
        Analysis.percentage_total_cycles(i) = round(percen(i),2);
    end
    
    max_cycle =  Analysis.last_cycle;
    max_verticalkm = max(Analysis.vertical_km);
    max_float_age = max(Analysis.float_age);         
end

%% Writting an information output text box %%%

msg1 = ['- You started with an input of ' num2str(length(Floats_list_ini.WMO)) ' floats'];
msg2 = ['- After connecting to the DAC, only ' num2str(length(Floats.WMO.data)) ' floats were found'];
warning_string_1 = ['- Number of floats whithout the input parameter requested: ' num2str(sum(idx_empty))];
warning_string_2= ['- WMO, models and DAC of the floats without the input param requested:'];
warning_string_3 = strcat(empty_WMO,',',empty_models,',',empty_DAC);
warning_string_3 = string(warning_string_3(:,:));
skip_line_str = ['--------------'];
max_cycle_string = ['The maximum number of cycle is : ' num2str(max(Analysis.last_cycle)) ' cycles'];
max_vkm_string = ['The maximum vertical km traveled is : ' num2str(round(max(Analysis.vertical_km)),2) ' kms'];
max_age_string= ['The maximum age reached is : ' num2str(round(max(Analysis.float_age(find(Analysis.float_age < 10)))),2) ' years'];
bottom_msg_string = ['Click the OK button below to continue'];

all_str = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s',msg1,msg2,warning_string_1,warning_string_2,warning_string_3,skip_line_str,...
    max_cycle_string,max_vkm_string,max_age_string,bottom_msg_string);

% info_msg = msgbox({msg1;msg2;warning_string_1;warning_string_2;warning_string_3;skip_line_str;max_cycle_string;max_vkm_string;...
%     max_age_string;skip_line_str;bottom_msg_string},'Main informations on your Data', 'modal');

info_msg = msgbox(all_str,'Main informations on your Data', 'modal');
uiwait(info_msg)    % Make the script pause until the user click ok on the display window

%% Creating the output directory

formatOut = 'ddmmyyyy';
working_date = datestr(now,formatOut);

newdir = [output_folder '/' calcul_param{1} '/' working_date];
if ~exist(newdir, 'dir')
  mkdir(newdir);
end

newdir2 = [output_folder '/' calcul_param{1} '/' working_date '/WMO_masks'];
if ~exist(newdir2, 'dir')
  mkdir(newdir2);
end

inputs_str = sprintf('%s%s\n%s%s\n%s%f\n%s%f\n%s%f\n%s\n','parameter: ',calcul_param{1},'parameter file type: ',where_file{1},'Parameter is boolean: ',is_boolean,...
                     'Sample size limit: ',sample_size_limit,'Plot grouping per models: ',plot_group_model,skip_line_str);
export_str1 =[inputs_str,all_str];


if is_export == 1
    filename = [output_folder '/' calcul_param{1} '/' working_date '/' 'information_message.txt'];
    fid = fopen(filename, 'wt');
    fprintf(fid, export_str1);
    fclose(fid);
end

disp('Information display over')


%% Test Luca %%
disp('computing...')
% close all

%%% Chose the value of the config param wanted to compare survival rate %%%
index_values = {};
cycles_vector = {};
verticalkm_vector = {};
floatage_vector = {};
floatage_vector2 = {};
max_float_age2 = {};
WMO_vector = {};
MODEL_vector = {};
status_vector = {};
sample_size = {};
plot_label = {};
death_floats = {};
alive_index = {};
MODEL_vector2 = {};
max_cycle={};
vect_cycle = {};
alive_idx = {};
sample_s_test = {};
floats_survived = {};
plot_data_cycle = {};
plot_data_table_cycle= {};
cycle_numbers= {};
table_cycle_numbers= {};
vect_vkm={};
alive_idx_vkm = {};
sample_s_vkm={};
plot_data_vkm={};
max_verticalkm= {};
plot_data_table_vkm= {};
vkm_traveled= {};
table_vkm= {};
floatage_intervals = {};
vect_age= {};
alive_idx_age= {};
sample_s_age= {};
plot_data_age= {};
plot_data_table_age= {};
age_reached= {};
table_age= {};
unique_model= {};
unique_model_nb= {};
model_idx={}; 
str_model={}; 
table_model_idx={};
sorted_param_values = {};



%% Loop used if the input parameter is not a CONFIG one and is not a boolean

if sum(contains(where_file, 'config')) == 0 && is_boolean == 0
    
    msg1 = ['According to your inputs, the parameter you chose is neither a config one, nor a boolean'];
    msg2 = ['You can either compute survival rates according to an interval defined on:'];
    choice1 = [' - The sum of this parameter value for all the floats cycles [enter 1 below]'];
    choice2 = [' - The mean value taken by this parameter over the cycles [enter 2 below]'];
    choice3 = [' - The parameter occured during cycle (groundings for example) [enter 3 below]'];
    skip_line_str = ['----------------------------------------------'];
    user_choice = ['User choice: either 1, 2 or 3'];
    
    all_str = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s',msg1,msg2,choice1,choice2,choice3,skip_line_str,user_choice);
    
    input_user = inputdlg(all_str,'User inputs required');
   
    choice_values_grouping = str2num(input_user{1});
    
    export_str2 = sprintf('%s\n%s\n%s\n%s\n%s\n%s%s%f\n',skip_line_str,choice1,choice2,choice3,skip_line_str,user_choice,': ',choice_values_grouping);
%     fid = fopen(filename, 'wt');
%     fprintf(fid, export_str2);
%     fclose(fid);
    
    if choice_values_grouping == 1
        max_value = max(Analysis.sum_value_param);
        
        msg21= ['You chose to define an interval based on the sum of the parameter values'];
        msg22= ['You are now going to define a threshold in order to form intervals depending on this sum'];
        msg23= ['The finer your threshold will be, the longer the computations will take'];
        skip_line_str = ['----------------------------------------------'];
        msg24= ['The maximum sum of the param value reached is : ' num2str(max_value)];
        user_choice = ['Threshold : '];

        all_str = sprintf('%s\n%s\n%s\n%s\n%s\n%s',msg21,msg22,msg23,msg24,skip_line_str,user_choice);
        input_threshold = inputdlg(all_str,'Interval threshold selection');
        
        Thresh_1 = str2num(input_threshold{1});
        
    elseif choice_values_grouping ==2
        max_mean = max(Analysis.mean_value_param);
        
        msg21= ['You chose to define an interval based on the mean of the parameter value per cycle'];
        msg22= ['You are now going to define a threshold in order to form intervals depending on this average value'];
        msg23= ['The finer your threshold will be, the longer the computations will take'];
        skip_line_str = ['----------------------------------------------'];
        msg24= ['The maximum average reached is: ' num2str(max_mean)];
        user_choice = ['Threshold : '];

        all_str = sprintf('%s\n%s\n%s\n%s\n%s\n%s',msg21,msg22,msg23,msg24,skip_line_str,user_choice);
        input_threshold = inputdlg(all_str,'Interval threshold selection');
        
        Thresh_1 = str2num(input_threshold{1});
       
   %%% For choice 3 : just count if the parameter occured during cycle (0 or 1), no threshold needed since the index will only take form on a 0,1 matrix
    
    end
    
    
    if choice_values_grouping == 3
        
    else
        
        export_str3 = sprintf('%s\n%s\n%s\n%s%f\n',skip_line_str,msg24,skip_line_str,user_choice,Thresh_1);
        export_str_all = [export_str1, export_str2, export_str3];
        
        if is_export == 1
            fid = fopen(filename, 'wt');
            fprintf(fid, export_str_all);
            fclose(fid);
        end
    end
    
    
    if choice_values_grouping == 1
        %%% Definition of the intervals to compute the survival rate if the input param is not a configuration one and is not a boolean :
        interval = 0:Thresh_1:max(Analysis.sum_value_param);
        for ivalue = 1:length(interval)

            if ivalue<length(interval)
                index_values{ivalue} = find(Analysis.sum_value_param > interval(ivalue) & Analysis.sum_value_param < interval(ivalue+1));
            else ivalue == length(interval)
                index_values{ivalue} = find(Analysis.sum_value_param == interval(ivalue));
            end

        end
    end
    
    if choice_values_grouping == 2
        %%% Definition of the intervals to compute the survival rate if the input param is not a configuration one and is not a boolean :
        interval = 0:Thresh_1:max(Analysis.mean_value_param);
        for ivalue = 1:length(interval)

            if ivalue<length(interval)
                index_values{ivalue} = find(Analysis.mean_value_param > interval(ivalue) & Analysis.mean_value_param < interval(ivalue+1));
            else ivalue == length(interval)
                index_values{ivalue} = find(Analysis.mean_value_param == interval(ivalue));
            end

        end
    end
    
    if choice_values_grouping == 3
        index_values{1} = find(Analysis.sum_value_param  == 0); %% Finding floats were the parameter did not occured during life of the float
        index_values{2} = find(Analysis.sum_value_param > 0); %% Finding floats were the parameter did occur during the life of the float
    end
    
    
    if choice_values_grouping ~= 3   
        %%% Dropping index values that are empty (no param values in certain intervals)
        interval = interval(~cellfun('isempty',index_values));
        index_values = index_values(~cellfun('isempty',index_values));
    
    
        %%% Dropping index_values that have a total float number under sample size limit threshold defined as input
        interval = interval(cellfun(@(x)length(x)>sample_size_limit, index_values));
        index_values = index_values(cellfun(@(x)length(x)>sample_size_limit, index_values));
    
    else
        
    end

%%
    for ivalue = 1:length(index_values)

        cycles_vector{ivalue} = Analysis.last_cycle(index_values{ivalue});
        verticalkm_vector{ivalue} = Analysis.vertical_km(index_values{ivalue});
        floatage_vector{ivalue} = Analysis.float_age(index_values{ivalue});

        max_cycle{ivalue} = max(cycles_vector{ivalue});
        max_verticalkm{ivalue}= max(verticalkm_vector{ivalue});
        
        %Mise en place d'un âge limite au-delà duquel on considère les valeurs d'âge erronnée (ici 10ans max)
        floatage_vector2{ivalue} = floatage_vector{ivalue};    
        for ii = 1:length(floatage_vector2{ivalue})
            if floatage_vector2{ivalue}(ii)>= 10
                floatage_vector2{ivalue}(ii) = NaN;
            end
        end
        
        
        max_float_age2{ivalue} = max(floatage_vector2{ivalue});
        
        WMO_vector{ivalue} = Floats_list.WMO(index_values{ivalue},:); 
        MODEL_vector{ivalue} = Floats_list.MODEL(index_values{ivalue},:);
        status_vector{ivalue} = Floats_list.STATUS(index_values{ivalue},:);
        
        sample_size{ivalue} = length(index_values{ivalue});
        
        if choice_values_grouping == 3
            plot_label{1} = {"Parameter did not occur during float lifetime"};
            plot_label{2} = {"Parameter did occur at least once during float lifetime"};
        else
            if ivalue < length(interval)
                plot_label{ivalue} = {[num2str(interval(ivalue)) ' <= PARAM <= ' num2str(interval(ivalue +1)) ' (' num2str(sample_size{ivalue}) ' floats)']};
            elseif ivalue == length(interval)
                plot_label{ivalue} = {[ 'PARAM <= ' num2str(interval(ivalue)) ' (' num2str(sample_size{ivalue}) ' floats)']};
            end
        end

        death_floats{ivalue} = sum(contains(cellstr(status_vector{ivalue}), 'INACTIVE')) + sum(contains(cellstr(status_vector{ivalue}), 'CLOSED'));
        alive_index{ivalue} = ~(contains(cellstr(status_vector{ivalue}), 'INACTIVE') | contains(cellstr(status_vector{ivalue}), 'CLOSED'));
    
        
        %%% Create a table indexing which float is alive or dead according to the model and WMO
        WMO_vector2{ivalue} = cellstr(WMO_vector{ivalue});

        if choice_values_grouping ~= 3
            if ivalue<length(interval)
                str_interval = ['interval' num2str(interval(ivalue)) '_' num2str(interval(ivalue +1))];
            elseif ivalue == length(interval)
                str_interval = ['interval' num2str(interval(ivalue))];
            end
        end
        
        MODEL_vector2{ivalue} = cellstr(MODEL_vector{ivalue});
        unique_model{ivalue} = unique(MODEL_vector2{ivalue});
        unique_model_nb{ivalue}= length(unique_model{ivalue});
        
        
       
    end
    
    
    
    for ivalue = 1:length(index_values)
        
        for imodel = 1:unique_model_nb{ivalue}
            
            model_idx{ivalue}{imodel} = (strcmp(cellstr(MODEL_vector2{ivalue}), unique_model{ivalue}(imodel)));
            str_model{ivalue}(imodel) = [string(unique_model{ivalue}(imodel))];
            table_model_idx{ivalue} = table(model_idx{ivalue}{:});
            
            
            WMO_vector3{ivalue}{imodel} = WMO_vector2{ivalue};
            table_WMO_idx{ivalue} = table(WMO_vector3{ivalue}{:});
        end

        cell_string{ivalue} = cellstr(str_model{ivalue});
        for imodel = 1:unique_model_nb{ivalue}
            disp(cell_string{ivalue}{imodel})
            table_WMO_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};
            table_model_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};            
        end
        
        % Duplicate the table of model index in order to update it later without names of models that contains less than sample size limit
        table_model_idx2{ivalue} = table_model_idx{ivalue};
        table_WMO_idx2{ivalue} = table_WMO_idx{ivalue};
        
        for imodel = 1:unique_model_nb{ivalue}
            
            % Creation of a variable containing the names of the models to delete from the table of model because they have less than the sample size limit
            if sum(table_model_idx{ivalue}.(str_model{ivalue}{imodel})) < sample_size_limit
                vars2delete{ivalue}{imodel} = str_model{ivalue}{imodel};
            else
                vars2delete{ivalue}{imodel}= {};
            end
            vars2delete{ivalue} = vars2delete{ivalue}(cellfun(@(x) ~isempty(x), vars2delete{ivalue}));
        end
        % Updating the table of models index by deleting fields (models)
        for j= 1:length(vars2delete{ivalue})
            table_WMO_idx2{ivalue} = removevars(table_WMO_idx2{ivalue}, vars2delete{ivalue}{j});
            table_model_idx2{ivalue} = removevars(table_model_idx2{ivalue}, vars2delete{ivalue}{j});
        end
        
        
        % updating the cell string of models by deleting floats models not meeting the sample size requirement
        cell_string{ivalue} = table_model_idx2{ivalue}.Properties.VariableNames;
        
    end
    
end

%% Loop used if the input parameter is not a CONFIG one but is a boolean

if sum(contains(where_file, 'config')) == 0 && is_boolean == 1
    
    %%% Definition of the intervals to compute the survival rate if the input param is not a configuration one and is a boolean (intervals on % of cycles affected) :
%     interval2 = 0:Thresh_2:max(Analysis.percentage_total_cycles);

    msg1 = ['According to your inputs, the parameter you chose is a boolean'];
    msg2 = ['You will compute survival rates according to an interval defined on:'];
    choice1 = ['the percentage of cycles affected by this value (value were the parameter = 1)'];
    skip_line_str = ['----------------------------------------------'];
    msg22= ['You are now going to define a threshold in order to form intervals'];
    msg23= ['The finer your threshold will be, the longer the computations will take'];
    msg24= ['The interval is defined on a percentage of the cycles affected, maximum is : 100%'];
    user_choice = ['Threshold : '];

    all_str = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s',msg1,msg2,choice1,skip_line_str,msg22,msg23,msg24,user_choice);
    
    input_threshold = inputdlg(all_str,'Interval threshold selection');

    Thresh_2 = str2num(input_threshold{1});
    

    interval2 = 0:Thresh_2:100;
    for ivalue = 1:length(interval2)
        
        if ivalue<length(interval2)
            index_values{ivalue} = find(Analysis.percentage_total_cycles >= interval2(ivalue) & Analysis.percentage_total_cycles <= interval2(ivalue+1));
        elseif ivalue == length(interval2)
            continue
%             index_values{ivalue} = find(Analysis.percentage_total_cycles > interval2(ivalue-1) & Analysis.percentage_total_cycles < interval2(ivalue));
        end

    end
    
    %%% Droping the first value of the interval2 variable (0) in order to have the same dimension as the index_values one
    %%% However, we keep the interval2 variable intact as it will be used for the naming of the plots
    interval2b = interval2(2:end);
    
    %%% Dropping index values that are empty (no param values in certain intervals)
    index_values = index_values(~cellfun('isempty',index_values));
    interval2b = interval2b(~cellfun('isempty',index_values));
    
    %%% Dropping index_values that have a total float number under sample size limit threshold defined as input
    index_values = index_values(cellfun(@(x)length(x)>sample_size_limit, index_values));
    interval2b = interval2b(cellfun(@(x)length(x)>sample_size_limit, index_values));
    
    for ivalue = 1:length(index_values)

        cycles_vector{ivalue} = Analysis.last_cycle(index_values{ivalue});
        verticalkm_vector{ivalue} = Analysis.vertical_km(index_values{ivalue});
        floatage_vector{ivalue} = Analysis.float_age(index_values{ivalue});

        max_cycle{ivalue} = max(cycles_vector{ivalue});
        max_verticalkm{ivalue}= max(verticalkm_vector{ivalue});
        
        %Mise en place d'un âge limite au-delà duquel on considère les valeurs d'âge erronnée (ici 10ans max)
        floatage_vector2{ivalue} = floatage_vector{ivalue};    
        for ii = 1:length(floatage_vector2{ivalue})
            if floatage_vector2{ivalue}(ii)>= 10
                floatage_vector2{ivalue}(ii) = NaN;
            end
        end
        max_float_age2{ivalue} = max(floatage_vector2{ivalue});
        
        WMO_vector{ivalue} = Floats_list.WMO(index_values{ivalue},:); 
        MODEL_vector{ivalue} = Floats_list.MODEL(index_values{ivalue},:);
        status_vector{ivalue} = Floats_list.STATUS(index_values{ivalue},:);
        
        sample_size{ivalue} = length(index_values{ivalue});
        if ivalue < length(interval2)
            plot_label{ivalue} = {[num2str(interval2(ivalue)) ' <= PARAM <= ' num2str(interval2(ivalue +1)) ' (' num2str(sample_size{ivalue}) ' floats)']};
        elseif ivalue == length(interval2)
            continue
%             plot_label{ivalue} = {[ 'PARAM = ' num2str(interval2(ivalue)) ' (' num2str(sample_size{ivalue}) ' floats)']};
        end

        death_floats{ivalue} = sum(contains(cellstr(status_vector{ivalue}), 'INACTIVE')) + sum(contains(cellstr(status_vector{ivalue}), 'CLOSED'));
        alive_index{ivalue} = ~(contains(cellstr(status_vector{ivalue}), 'INACTIVE') | contains(cellstr(status_vector{ivalue}), 'CLOSED'));
    
        
    %%% Create a table indexing which float is alive or dead according to the model
        
        MODEL_vector2{ivalue} = cellstr(MODEL_vector{ivalue});
        unique_model{ivalue} = unique(MODEL_vector2{ivalue});
        unique_model_nb{ivalue}= length(unique_model{ivalue});
        
        WMO_vector2{ivalue} = cellstr(WMO_vector{ivalue});

        for imodel = 1:unique_model_nb{ivalue}
            model_idx{ivalue}{imodel} = (contains(cellstr(MODEL_vector2{ivalue}), unique_model{ivalue}(imodel)));
            str_model{ivalue}(imodel) = [string(unique_model{ivalue}(imodel))];
            table_model_idx{ivalue} = table(model_idx{ivalue}{:});
            
            WMO_vector3{ivalue}{imodel} = WMO_vector2{ivalue};
            table_WMO_idx{ivalue} = table(WMO_vector3{ivalue}{:});
        end

        cell_string{ivalue} = cellstr(str_model{ivalue});
        for imodel = 1:unique_model_nb{ivalue}
            table_model_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};
            table_WMO_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};
        end
        
        % Duplicate the table of model index in order to update it later without names of models that contains less than sample size limit
        table_model_idx2{ivalue} = table_model_idx{ivalue};
        table_WMO_idx2{ivalue} = table_WMO_idx{ivalue};
        
        for imodel = 1:unique_model_nb{ivalue}
            
            % Creation of a variable containing the names of the models to delete from the table of model because they have less than the sample size limit
            if sum(table_model_idx{ivalue}.(str_model{ivalue}{imodel})) < sample_size_limit
                vars2delete_idx{ivalue}{imodel} = str_model{ivalue}{imodel};
            else
                vars2delete_idx{ivalue}{imodel} = {};
            end
            var2delete_idx{ivalue} = cellfun(@(x) ~isempty(x), vars2delete_idx{ivalue});
            interval2delete = cellfun(@sum, var2delete_idx);
            all_intervals = cellfun(@width,table_model_idx);
            
            if interval2delete(ivalue) - all_intervals(ivalue) == 0
                empty_interval(ivalue) = 1;
            else
                empty_interval(ivalue) = 0;
            end
            
            
            vars2delete{ivalue} = vars2delete_idx{ivalue}(cellfun(@(x) ~isempty(x), vars2delete_idx{ivalue}));
        end
        
        % Updating the table of models index by deleting fields (models)
        for j= 1:length(vars2delete{ivalue})
%             if length(vars2delete{ivalue}) == size(table_model_idx2{ivalue},2)
%                 disp('The number of floats for all the models within this interval is below the sample limit: Deleting table entry')
%                 table_model_idx2{ivalue} = [];
%             else    
            table_model_idx2{ivalue} = removevars(table_model_idx2{ivalue}, vars2delete{ivalue}{j});
            table_WMO_idx2{ivalue} = removevars(table_WMO_idx2{ivalue}, vars2delete{ivalue}{j});
%             end    
        end
        
        % updating the cell string of models by deleting floats models not meeting the sample size requirement
        cell_string{ivalue} = table_model_idx2{ivalue}.Properties.VariableNames;
        
%         if isempty(table_model_idx2{ivalue}) ==1
%             table_model_idx2(ivalue) = []; %deleting interval from table if none of the floats model is above sample limit
%         end
        

        
        % updating intervals string
        cell_interval2 = interval2b(~empty_interval);
        cell_interval2 = strsplit(num2str(cell_interval2));
      
    end
end

%% Loop used if the input parameter is a CONFIG one

if sum(contains(where_file, 'config')) ~= 0
    for ivalue = 1:length(unique_values)
%          index_values{ivalue} = [notchan_configvalue == unique_values(ivalue)];
        index_values{ivalue} = find(notchan_configvalue == unique_values(ivalue));
    end
    
    
    %%% Dropping index_values that have a total float number under sample size limit threshold defined as input
    index_values2 = index_values(cellfun(@(x) length(x)> sample_size_limit, index_values));
    unique_values2 = unique_values(cellfun(@(x) length(x)> sample_size_limit, index_values));
    
    for ivalue= 1:length(index_values2)
%         index_values{ivalue} = [notchan_configvalue == unique_values(ivalue)];
        cycles_vector{ivalue} = notchan_last_cycle(index_values2{ivalue});
        verticalkm_vector{ivalue} = notchan_verticalkm(index_values2{ivalue});
        floatage_vector{ivalue} = notchan_float_age(index_values2{ivalue});
        
        

        
        max_cycle{ivalue} = max(cycles_vector{ivalue});
        max_verticalkm{ivalue}= max(verticalkm_vector{ivalue});

        %Mise en place d'un âge limite au-delà duquel on considère les valeurs d'âge erronnée (ici 10ans max)
        floatage_vector2{ivalue} = floatage_vector{ivalue};    
        for ii = 1:length(floatage_vector2{ivalue})
            if floatage_vector2{ivalue}(ii)>= 10;
                floatage_vector2{ivalue}(ii) = NaN;
            end
        end
        max_float_age2{ivalue} = max(floatage_vector2{ivalue});

        WMO_vector{ivalue} = notchan_WMO(index_values2{ivalue},:); 
        MODEL_vector{ivalue} = notchan_MODEL(index_values2{ivalue},:);
        status_vector{ivalue} = notchan_status(index_values2{ivalue},:);

        sample_size{ivalue} = sum(index_values2{ivalue});
        plot_label{ivalue} = {['CONFIG = ' num2str(unique_values2(ivalue)) ' (' num2str(sample_size{ivalue}) ' floats)']};

        death_floats{ivalue} = sum(contains(cellstr(status_vector{ivalue}), 'INACTIVE')) + sum(contains(cellstr(status_vector{ivalue}), 'CLOSED'));
        alive_index{ivalue} = ~(contains(cellstr(status_vector{ivalue}), 'INACTIVE') | contains(cellstr(status_vector{ivalue}), 'CLOSED'));

    %%% Create a table indexing which float is alive or dead according to the model
        MODEL_vector2{ivalue} = cellstr(MODEL_vector{ivalue});
        WMO_vector2{ivalue} = cellstr(WMO_vector{ivalue});
        unique_model{ivalue} = unique(MODEL_vector2{ivalue});
        unique_model_nb{ivalue}= length(unique_model{ivalue});

        for imodel = 1:unique_model_nb{ivalue}
            model_idx{ivalue}{imodel} = (contains(cellstr(MODEL_vector2{ivalue}), unique_model{ivalue}(imodel)));
            str_model{ivalue}(imodel) = [string(unique_model{ivalue}(imodel))];
            table_model_idx{ivalue} = table(model_idx{ivalue}{:});
            
            WMO_vector3{ivalue}{imodel} = WMO_vector2{ivalue};
            table_WMO_idx{ivalue} = table(WMO_vector3{ivalue}{:});
        end

        cell_string{ivalue} = cellstr(str_model{ivalue});
        for imodel = 1:unique_model_nb{ivalue}
            table_model_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};
            table_WMO_idx{ivalue}.Properties.VariableNames{imodel} = cell_string{ivalue}{imodel};
        end
        
        % Duplicate the table of model index in order to update it later without names of models that contains less than sample size limit
        table_model_idx2{ivalue} = table_model_idx{ivalue};
        table_WMO_idx2{ivalue} = table_WMO_idx{ivalue};
        
        for imodel = 1:unique_model_nb{ivalue}
            
            % Creation of a variable containing the names of the models to delete from the table of model because they have less than the sample size limit
            if sum(table_model_idx{ivalue}.(str_model{ivalue}{imodel})) < sample_size_limit
                vars2delete{ivalue}{imodel} = str_model{ivalue}{imodel};
            else
                vars2delete{ivalue}{imodel} = {};
            end
            vars2delete{ivalue} = vars2delete{ivalue}(cellfun(@(x) ~isempty(x), vars2delete{ivalue}));
        end
        
        % Updating the table of models index by deleting fields (models)
        for j= 1:length(vars2delete{ivalue})
            table_model_idx2{ivalue} = removevars(table_model_idx2{ivalue}, vars2delete{ivalue}{j});
            table_WMO_idx2{ivalue} = removevars(table_WMO_idx2{ivalue}, vars2delete{ivalue}{j});
        end
    end
    
    
    % Deleting intervals that became empty after deleting the models under the sample limit size in them.
%     table_model_idx2 = table_model_idx2(cellfun(@(x) ~isempty(x), table_model_idx2));
%     table_WMO_idx2 = table_WMO_idx2(cellfun(@(x) ~isempty(x), table_WMO_idx2));


%   for ivalue = 1:length(table_model_idx2)
%       if size(table_model_idx2{ivalue},2) == 0
%           table_model_idx2{ivalue} = [];
%       end
%   end
%   table_model_idx2 = table_model_idx2(~cellfun('isempty',table_model_idx2));
  
  for ivalue = 1:length(table_model_idx2)  
      % updating the cell string of models by deleting floats models not meeting the sample size requirement
      cell_string{ivalue} = table_model_idx2{ivalue}.Properties.VariableNames;   
  end  

end


%% Creation of plots arrays 
for ivalue = 1:length(table_model_idx2)
%     if isempty(table_model_idx2{ivalue}) ==1
%         disp('isempty')
%         continue
%     end
    disp(ivalue)

%     for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
    for imodel = 1:size(table_model_idx2{ivalue},2)    
%         if isempty(imodel) == 1
%             continue
%         end
        
        for icycle = 1:max_cycle{ivalue}
            vect_cycle{ivalue}{imodel}= cycles_vector{ivalue}(table_model_idx2{ivalue}{:,imodel});
            alive_idx_cycle{ivalue}{imodel} = alive_index{ivalue}(table_model_idx2{ivalue}{:,imodel});
            sample_s_cycle{ivalue}{icycle,imodel} =  sum(vect_cycle{ivalue}{imodel}(~alive_idx_cycle{ivalue}{imodel}) <= icycle) + sum(vect_cycle{ivalue}{imodel} >= icycle);
            floats_survived_cycle{ivalue}{icycle,imodel} = sum(vect_cycle{ivalue}{imodel} >= icycle);
            plot_data_cycle{ivalue}{icycle,imodel} = sum(vect_cycle{ivalue}{imodel} >= icycle)/sample_s_cycle{ivalue}{icycle,imodel}*100;
        end
    end
   


    plot_data_table_cycle{ivalue} = cell2table(plot_data_cycle{ivalue}, 'VariableNames', cell_string{ivalue});
    cycle_numbers{ivalue} = (1:max_cycle{ivalue})';
    table_cycle_numbers{ivalue} = table(cycle_numbers{ivalue}, 'VariableNames',{'Cycle_Number'});
    %%% Final table for cycles below:
    plot_data_table_cycle{ivalue} = [table_cycle_numbers{ivalue},plot_data_table_cycle{ivalue}]


    %%% loop vertical km
    % vkm_intervals = 0:10:max_verticalkm; %%%% Possibility to create some verticalkm intervals to simplify lisibility (every 10 Kms or something like that.
    for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
        for ivkm = 1:max_verticalkm{ivalue}
            vect_vkm{ivalue}{imodel} = verticalkm_vector{ivalue}(table_model_idx2{ivalue}{:,imodel});
            alive_idx_vkm{ivalue}{imodel} = alive_index{ivalue}(table_model_idx2{ivalue}{:,imodel});
            sample_s_vkm{ivalue}{ivkm,imodel} =  sum(vect_vkm{ivalue}{imodel}(~alive_idx_vkm{ivalue}{imodel}) <= ivkm) + sum(vect_vkm{ivalue}{imodel} >= ivkm);
            plot_data_vkm{ivalue}{ivkm,imodel} = sum(vect_vkm{ivalue}{imodel} > ivkm)/sample_s_vkm{ivalue}{ivkm,imodel}*100;
        end
    end

    plot_data_table_vkm{ivalue} = cell2table(plot_data_vkm{ivalue}, 'VariableNames', cell_string{ivalue});
    vkm_traveled{ivalue} = (1:max_verticalkm{ivalue})';
    table_vkm{ivalue} = table(vkm_traveled{ivalue}, 'VariableNames', {'Vertical_km_traveled'});
    %%% Final table for vertical distance below:
    plot_data_vkm{ivalue} = [table_vkm{ivalue},plot_data_table_vkm{ivalue}]


    %%% loop float age
    floatage_intervals{ivalue} = 0:0.1:max_float_age2{ivalue};
    for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
        for iage = 1:length(floatage_intervals{ivalue})
            vect_age{ivalue}{imodel} = floatage_vector{ivalue}(table_model_idx2{ivalue}{:,imodel});
            alive_idx_age{ivalue}{imodel} = alive_index{ivalue}(table_model_idx2{ivalue}{:,imodel});
            sample_s_age{ivalue}{iage,imodel} =  sum(vect_age{ivalue}{imodel}(~alive_idx_age{ivalue}{imodel}) <= floatage_intervals{ivalue}(iage)) + ... 
            sum(vect_age{ivalue}{imodel} >= floatage_intervals{ivalue}(iage));
            plot_data_age{ivalue}{iage,imodel} = sum(vect_age{ivalue}{imodel} > floatage_intervals{ivalue}(iage))/sample_s_age{ivalue}{iage,imodel}*100;
        end   
    end

    plot_data_table_age{ivalue} = cell2table(plot_data_age{ivalue}, 'VariableNames', cell_string{ivalue});
    age_reached{ivalue} = (floatage_intervals{ivalue})';
    table_age{ivalue} = table(age_reached{ivalue}, 'VariableNames', {'Age_reached'});
    %%% Final table for age below:
    plot_data_age{ivalue} = [table_age{ivalue},plot_data_table_age{ivalue}]
end

%% Re-arranging the plot arrays in order to make one plot per model and display the survival rates in function
%%% of the parameter value taken

if plot_group_model == 1
    disp('grouping survival rates per models')

    all_models = unique(cellstr(Floats_list.MODEL));

    A_cycles= {};
    A_age_reached = {};
    A_plot_table_age= {};
    A_table_age={};
    A_plot_table_age={};
    A_age= {};
    A_vkm={};

    for ivalue= 1:length(table_model_idx2)
        %%% le cell string ne doit pas changer d'un ivalue à l'autre. On doit
        %%% prendre sans cesse le meme modèle, donc cell_string = {6x1} :
        %%% {Arvor} {APEX} , etc.. et si ivalue out of range de cell_string, skip
    %     cell_string{ivalue} = cellstr(str_model{ivalue});
        disp('ivalue =')
        disp(ivalue)
        for imodel = 1:length(all_models)
            if any(strcmp(all_models{imodel},plot_data_table_cycle{ivalue}.Properties.VariableNames)) == 0
                disp('No data for this model in this interval:')
                disp(imodel)
                A_WMO{imodel}{ivalue} = {};
                B_WMO{ivalue}{imodel} = {};
                
                A_cycles{imodel}{ivalue}= {};
                A_age{imodel}{ivalue}= {};
                A_vkm{imodel}{ivalue}= {};
            else
                disp('Data found for this model, within this range:')
                disp(imodel)

                A_WMO{imodel}{ivalue} = table_WMO_idx2{ivalue}.(all_models{imodel})(table_model_idx2{ivalue}.(all_models{imodel}));
                B_WMO{ivalue}{imodel} = table_WMO_idx2{ivalue}.(all_models{imodel})(table_model_idx2{ivalue}.(all_models{imodel}));
                
                A_cycles{imodel}{ivalue} = plot_data_table_cycle{ivalue}.(all_models{imodel});
                A_age_reached{imodel}{ivalue} = (floatage_intervals{ivalue})';
%                 A_table_age_reached{imodel} = cell2table(A_age_reached{imodel});
                A_age{imodel}{ivalue} = plot_data_age{ivalue}.(all_models{imodel});
%                 A_table_age{imodel} = cell2table(A_age{imodel});
%                 A_plot_table_age{imodel}{ivalue} = [A_age_reached{imodel}{ivalue},A_table_age{imodel}];
                
                A_vkm{imodel}{ivalue} = plot_data_vkm{ivalue}.(all_models{imodel});
            end

            A_WMO2{imodel} = cat(1,A_WMO{imodel}{:});
            table_WMO_per_model{imodel} = cell2table(A_WMO2{imodel});
            
            B_WMO2{ivalue} = cat(1,B_WMO{ivalue}{:});
            table_WMO_per_interval{ivalue} = cell2table(B_WMO2{ivalue});
            
            if isempty(table_WMO_per_model{imodel}) == 0
                table_WMO_per_model{imodel}.Properties.VariableNames = {'REF'};
            end
            
            if isempty(table_WMO_per_interval{ivalue}) == 0
                table_WMO_per_interval{ivalue}.Properties.VariableNames = {'REF'};
            end
            
            filename = [output_folder '/' calcul_param{1} '/' working_date '/WMO_masks/' 'WMO_mask_for_model_' all_models{imodel} '.csv'];
            writetable(table_WMO_per_model{imodel}, filename, 'Delimiter', ';')
            
            if sum(contains(where_file, 'config')) == 0
                if choice_values_grouping == 3
                    filename2 = [output_folder '/' calcul_param{1} '/' working_date '/WMO_masks/' 'WMO_mask_for_groundings' num2str(ivalue) '.csv'];
                    writetable(table_WMO_per_interval{ivalue}, filename2, 'Delimiter', ';')
                else
                    filename2 = [output_folder '/' calcul_param{1} '/' working_date '/WMO_masks/' 'WMO_mask_for_interval_' num2str(ivalue) '.csv'];
                    writetable(table_WMO_per_interval{ivalue}, filename2, 'Delimiter', ';')
                end
            end
            
        end     
    end
    


    
end


%%
%%% PLOTS

%% Specific plots if the parameter is not a CONFIG one and not a boolean 
%%% These plots are specifically designed with titles and legends to plot survival rates grouped by interval of values taken by the input parameter

for ivalue = 1:length(table_model_idx2)
    if sum(contains(where_file, 'config')) == 0 && is_boolean == 0 && plot_group_model==0
        %%% Plot
        figure(ivalue)
        
        % full screen
        set(gcf, 'Position', get(0, 'Screensize'));
        %     clear title
        lgd_str={};
        super_title=[];
        %     title_str= {};

        
        %     set(gcf, 'Position', get(0, 'Screensize'));
        total_floats = size(WMO_vector{ivalue},1);
        
        % Verify if the total floats for the parameter value interval is above threshold defined in input (default is 10)
        if total_floats < sample_size_limit
            disp('The number of floats within the parameter value interval is below the sample size limit threshold')
            continue
        end
        
        if ivalue < length(index_values)
            super_title= [interval(ivalue) + " <= " + cellstr(calcul_param) + " <= " + interval(ivalue+1) + " [" + total_floats + " floats total]"];
        elseif ivalue == length(index_values)
            super_title= [cellstr(calcul_param) + " > " + interval(ivalue) + " [" + total_floats + " floats total]"];
        end
        
        annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 18, 'FitBoxToText', 'on', 'Interpreter', 'none')

        subplot(2,2,1)
        %Survival rates per cycles
        %     plot_data_table_cycle2 = plot_data_table_cycle;
        %     plot_data_table_cycle2{5}.PROVOR_II = [];
        %     plot_data_table_cycle2{5}.PROVOR_III = [];
        %     plot_data_table_cycle2{10}.PROVOR_II = [];
        % % %     
        %     unique_model_nb2 = unique_model_nb
        %     unique_model_nb2{5} = unique_model_nb{5} -2
        %     unique_model_nb2{10} = unique_model_nb{10} -1

        for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_table_cycle{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Number of cycles reached')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per cycle per float type')

        subplot(2,2,2)
        %Survival rates per age
        %     plot_data_age2 = plot_data_age;
        %     plot_data_age2{5}.PROVOR_II = [];
        %     plot_data_age2{5}.PROVOR_III = [];
        %     plot_data_age2{10}.PROVOR_II = [];

        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_age{ivalue}{:,1}, plot_data_age{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        xlim([0 max_float_age2{ivalue}])
        ylim([0 100])
        xlabel('Age reached (year)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per age reached per float type')

        subplot(2,2,3:4)
        %Survival rate per vertical distance traveled
        %     plot_data_vkm2 = plot_data_vkm;
        %     plot_data_vkm2{5}.PROVOR_II = [];
        %     plot_data_vkm2{5}.PROVOR_III = [];
        %     plot_data_vkm2{10}.PROVOR_II = [];
        %     
        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_vkm{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Vertical distance traveled (Km)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per vertical Kms traveled per float type')
        
       %%% export
        if is_export == 1
            for ivalue = 1:length(index_values)
                if ivalue < length(index_values)
                    filename = [cellstr(calcul_param) + "_interval_" + interval(ivalue) + interval(ivalue+1) + ".jpg"];
                elseif ivalue == length(index_values)
                    filename = [cellstr(calcul_param) + "_interval_equal_to_" + interval(ivalue) + ".jpg"];
                end
            end
        end
    
    end
end
    
if sum(contains(where_file, 'config')) == 0 && is_boolean == 0 && plot_group_model == 1
    for imodel = 1:length(all_models)
        figure(imodel)
        
        %%% Plot
        % full screen
        set(gcf, 'Position', get(0, 'Screensize'));
        % White background
        set(gcf,'color','white')
%         set(0, 'DefaultFigureRenderer', 'opengl');
        %     clear title
        lgd_str={};
        super_title=[];
        %     title_str= {};

        



        A = cellstr(Floats_list.MODEL);
        for ii = 1:length(A)
            if strcmp(all_models{imodel},A{ii}) == 1
                sum_model{imodel}(ii) = 1;
            else
                sum_model{imodel}(ii) = 0;
            end
        end    
        total_floats_per_model{imodel} = sum(sum_model{imodel});
       
        super_title= ["Survival rates for " + all_models{imodel} +  " floats, per " + ...
            calcul_param{1} + " parameter values" + " [" + total_floats_per_model{imodel} + " floats total]" ];
        annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 16, 'FitBoxToText', 'on', 'Interpreter', 'none')

%%%Survival rates per cycles
        subplot(2,2,1)
        
        for ivalue = 1:length(table_model_idx)
            if choice_values_grouping == 3
                plot(A_cycles{imodel}{ivalue}, 'LineWidth', 2)
                lgd_str{ivalue} = [plot_label{ivalue}{:} + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)"];
                total_floats_considered{imodel}{ivalue} = sum(table_model_idx2{ivalue}.(all_models{imodel}));
                hold on
            else
                if ivalue < length(table_model_idx)
                    if isempty(A_cycles{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
                        plot(A_cycles{imodel}{ivalue}, 'LineWidth', 2)
                        lgd_str{ivalue} = "Parameter value within [" + interval(ivalue) + ";" + interval(ivalue+1) + "] (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
                        total_floats_considered{imodel}{ivalue} = sum(table_model_idx2{ivalue}.(all_models{imodel}));

                        hold on
                    else
                        plot(0,0,'x')
                        hold on
%                         continue
    %                     lgd_str{ivalue} = "No floats in the parameter value interval [" + interval(ivalue) + ";" + interval(ivalue+1) + "]";
    %                     total_floats_considered{imodel}{ivalue} = 0;
    %                     hold on
                    end

                elseif ivalue == length(table_model_idx)

    %                 plot(0,0,'rx')
    %                 lgd_str{ivalue} = "No floats within the [" + interval(ivalue) + ";" + interval(ivalue+1) + "]% interval";
                    if isempty(A_cycles{imodel}{ivalue}) == 0
                        total_floats_considered{imodel}{ivalue} = sum(table_model_idx2{ivalue}.(all_models{imodel}));
                    else
                        total_floats_considered{imodel}{ivalue} = 0;
                    end
                    hold on
                    continue
                end
            end     
        end
        
        total_floats_considered2{imodel} = sum(cellfun(@sum, total_floats_considered{imodel}));
        
        xlim([0 500])
        ylim([0 100])
        xlabel('Number of cycles reached')
        ylabel('Survival rate (%)')
%         legend(lgd_str{:})
        title('Survival rate per cycles made')

%%% Survival rates per age
        subplot(2,2,2)

        for ivalue = 1:length(table_model_idx)
            if choice_values_grouping == 3
                plot(A_age_reached{imodel}{ivalue}, A_age{imodel}{ivalue},'LineWidth', 2)
                lgd_str{ivalue} = [plot_label{ivalue}{:} + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)"];
                hold on
            else   
                if ivalue < length(table_model_idx)
                    if isempty(A_age{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
                        plot(A_age_reached{imodel}{ivalue},A_age{imodel}{ivalue}, 'LineWidth', 2)
                        lgd_str{ivalue} = "Parameter value within [" + interval(ivalue) + ";" + interval(ivalue+1) + "] (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";

                        hold on
                    else
                        plot(0,0,'x')
                        hold on
%                         continue
    %                     lgd_str{ivalue} = "No floats in the parameter value interval [" + interval(ivalue) + ";" + interval(ivalue+1) + "]";
    %                     hold on
                    end

                elseif ivalue == length(table_model_idx)
    %                 plot(0,0,'rx')
    %                 lgd_str{ivalue} = "No (or not enough) floats within the [" + interval(ivalue) + ";" + interval(ivalue+1) + "]% interval";
                    hold on
                    continue
                end
            end
        end
        xlim([0 (max(cellfun(@max, max_float_age2)))])
        ylim([0 100])
        xlabel('Age reached (year)')
        ylabel('Survival rate (%)')
%         legend(lgd_str{:})
        title('Survival rate per age reached')

%%% Survival rate per vertical distance traveled
        subplot(2,2,3:4)
  
        for ivalue = 1:length(table_model_idx)
            if choice_values_grouping == 3
                plot(A_vkm{imodel}{ivalue}, 'LineWidth', 2)
                lgd_str{ivalue} = [plot_label{ivalue}{:} + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)"];
                hold on
                
            else
                if ivalue < length(table_model_idx)
                    if isempty(A_vkm{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
                        plot(A_vkm{imodel}{ivalue}, 'LineWidth', 2)
                        lgd_str{ivalue} = "Parameter value within [" + interval(ivalue) + ";" + interval(ivalue+1) + "] (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";

                        hold on
                    else
                        
                        plot(0,0,'x')
                        lgd_str{ivalue} = ["No floats in the parameter value interval [" + interval(ivalue) + ";" + interval(ivalue+1) + "]"];
                        hold on
%                         continue
    %                     hold on
                    end

                elseif ivalue == length(table_model_idx)
    %                 plot(0,0,'rx')
    %                 lgd_str{ivalue} = "No floats within the [" + interval(ivalue) + ";" + interval(ivalue+1) + "]% interval";
                    hold on
                    continue
                end 
            end
        end
        ylim([0 100])
        xlabel('Vertical distance traveled (Km)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:}, 'Location','EastOutside')
        title('Survival rate per vertical Kms traveled')
        
%         super_title= ["Survival rates for " + all_models{imodel} +  " floats, per " + ...
%             param{1} + " parameter values" + " [" + total_floats_considered2{imodel} + "/" + total_floats_per_model{imodel} + "considered/input floats]" ];
%         annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 16, 'FitBoxToText', 'on', 'Interpreter', 'none')
        
        %%% export
        if is_export==1
            disp('exporting...')        
            filename = [all_models{imodel} + "_" + cellstr(calcul_param) + "_surv_rates_groundings_intervals_50.jpg"];
            fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/" + filename];        
            export_fig(fig_path)
        end
    end 
end


%% Specific plots if the input parameter is not a CONFIG one but is a boolean 
%%% These plots are specifically designed with titles and legends to plot survival rates grouped by percentage of cycles affected by the input parameter
for ivalue = 1:length(table_model_idx2)
    if sum(contains(where_file, 'config')) == 0 && is_boolean == 1 && plot_group_model == 0
        %%% Plot
        figure(ivalue)
        % full screen
        set(gcf, 'Position', get(0, 'Screensize'));
        %     clear title
        lgd_str={};
        super_title=[];
        %     title_str= {};

        
        %     set(gcf, 'Position', get(0, 'Screensize'));
        for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            total_floats_per_model{ivalue}{imodel} = sum(table_model_idx2{ivalue}{:,imodel});
        end
        total_floats{ivalue} = sum(cellfun(@sum, total_floats_per_model{ivalue}));

%         % Verify if the total floats for the parameter value interval is above threshold defined in input (default is 10)
%         if total_floats < sample_size_limit
%             disp('The number of floats within the parameter value interval is below the sample size limit threshold')
%             continue
%         end
        
%         if ivalue < length(index_values)
        super_title= [interval2(ivalue) + "%" + " <= " + " cycles where " + cellstr(calcul_param) + " is TRUE" + " <= " + interval2(ivalue+1) + "%" + " [" + total_floats{ivalue} + " floats total]"];
%         elseif ivalue == length(index_values)
%             super_title= ["Cycles where " + cellstr(calcul_param) + " is TRUE " + " = " + interval2(ivalue) + "%" + " [" + total_floats{ivalue} + " floats total]"];
%         end
        
        annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 18, 'FitBoxToText', 'on', 'Interpreter', 'none')

        subplot(2,2,1)
        %Survival rates per cycles
        %     plot_data_table_cycle2 = plot_data_table_cycle;
        %     plot_data_table_cycle2{5}.PROVOR_II = [];
        %     plot_data_table_cycle2{5}.PROVOR_III = [];
        %     plot_data_table_cycle2{10}.PROVOR_II = [];
        % % %     
        %     unique_model_nb2 = unique_model_nb
        %     unique_model_nb2{5} = unique_model_nb{5} -2
        %     unique_model_nb2{10} = unique_model_nb{10} -1

        for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_table_cycle{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Number of cycles reached')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per cycle per float type')

        subplot(2,2,2)
        %Survival rates per age
        %     plot_data_age2 = plot_data_age;
        %     plot_data_age2{5}.PROVOR_II = [];
        %     plot_data_age2{5}.PROVOR_III = [];
        %     plot_data_age2{10}.PROVOR_II = [];

        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_age{ivalue}{:,1}, plot_data_age{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        xlim([0 max_float_age2{ivalue}])
        ylim([0 100])
        xlabel('Age reached (year)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per age reached per float type')

        subplot(2,2,3:4)
        %Survival rate per vertical distance traveled
        %     plot_data_vkm2 = plot_data_vkm;
        %     plot_data_vkm2{5}.PROVOR_II = [];
        %     plot_data_vkm2{5}.PROVOR_III = [];
        %     plot_data_vkm2{10}.PROVOR_II = [];
        %     
        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_vkm{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Vertical distance traveled (Km)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per vertical Kms traveled per float type')
        
        %%% export
        if is_export==1
            for ivalue = 1:length(index_values)
                if ivalue < length(index_values)
                    filename = [string(calcul_param) + "_within_interval_" + interval2(ivalue) + "_" + interval2(ivalue+1) + ".jpg"];
                elseif ivalue == length(index_values)
                    filename = [string(calcul_param) + "_interval_equal_to_" + interval2(ivalue) + ".jpg"];
                end
            end
        end
    end
end

%%% If group model == 1
if sum(contains(where_file, 'config')) == 0 && is_boolean == 1 && plot_group_model == 1
    for imodel = 1:length(all_models)
        
        %%% Plot
        figure(imodel)
        % full screen
        set(gcf, 'Position', get(0, 'Screensize'));
        %     clear title
        lgd_str={};
        super_title=[];
        %     title_str= {};

        



        A = cellstr(Floats_list.MODEL);
        for ii = 1:length(A)
            if strcmp(all_models{imodel},A{ii}) == 1
                sum_model{imodel}(ii) = 1;
            else
                sum_model{imodel}(ii) = 0;
            end
        end    
        total_floats_per_model{imodel} = sum(sum_model{imodel});
       
        super_title= ["Survival rates for " + all_models{imodel} +  " floats, per " + ...
            calcul_param{1} + " parameter values" + " [" + total_floats_per_model{imodel} + "floats total]" ];
        annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 16, 'FitBoxToText', 'on', 'Interpreter', 'none')

%%%Survival rates per cycles
        subplot(2,2,1)
        
        for ivalue = 1:length(table_model_idx)
            if ivalue < length(table_model_idx)
                if isempty(A_cycles{imodel}{ivalue}) == 0
                    plot(A_cycles{imodel}{ivalue}, 'LineWidth', 2)
                    lgd_str{ivalue} = "Floats within [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% total cycles affected by param (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
                    hold on
                else
                    plot(0,0,'rx')
                    lgd_str{ivalue} = "No floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                    hold on
                end
                
            elseif ivalue == length(table_model_idx)
                plot(0,0,'rx')
                lgd_str{ivalue} = "No floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                hold on
            end
        end
        xlim([0 500])
        ylim([0 100])
        xlabel('Number of cycles reached')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per cycles made')

%%% Survival rates per age
        subplot(2,2,2)

        for ivalue = 1:length(table_model_idx)
            if ivalue < length(table_model_idx)
                if isempty(A_age{imodel}{ivalue}) == 0
                    plot(A_age_reached{imodel}{ivalue},A_age{imodel}{ivalue}, 'LineWidth', 2)
                    lgd_str{ivalue} = "Floats within [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% total cycles affected by param (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
                    hold on
                else
                    plot(0,0,'rx')
                    lgd_str{ivalue} = "No floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                    hold on
                end
                
            elseif ivalue == length(table_model_idx)
                plot(0,0,'rx')
                lgd_str{ivalue} = "No (or not enough) floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                hold on
            end
        end
        xlim([0 (max(cellfun(@max, max_float_age2)))])
        ylim([0 100])
        xlabel('Age reached (year)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per age reached')

%%% Survival rate per vertical distance traveled
        subplot(2,2,3:4)
  
        for ivalue = 1:length(table_model_idx)
            if ivalue < length(table_model_idx)
                if isempty(A_vkm{imodel}{ivalue}) == 0
                    plot(A_vkm{imodel}{ivalue}, 'LineWidth', 2)
                    lgd_str{ivalue} = "Floats within [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% total cycles affected by param (" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
                    hold on
                else
                    plot(0,0,'rx')
                    lgd_str{ivalue} = "No floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                    hold on
                end
                
            elseif ivalue == length(table_model_idx)
                plot(0,0,'rx')
                lgd_str{ivalue} = "No floats within the [" + interval2(ivalue) + ";" + interval2(ivalue+1) + "]% interval";
                hold on
            end
        end
        ylim([0 100])
        xlabel('Vertical distance traveled (Km)')
        ylabel('Survival rate (%)')
        legend(lgd_str{:})
        title('Survival rate per vertical Kms traveled')
        
%         
                
        %%% export
        if is_export==1
            disp('exporting...')        
            filename = [all_models{imodel} + "_" + cellstr(calcul_param) + "_surv_rates.jpg"];
            fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/" + filename];        
            export_fig(fig_path)
        end
    end
end    
    
%% Specific plots if the parameter is a CONFIG one
%%% These plots are specifically designed with titles and legends to plot survival rates grouped by values taken by the CONFIG parameter
for ivalue = 1:length(table_model_idx2)
    if sum(contains(where_file, 'config')) ~= 0 && plot_group_model==0
        %%% Plot
        figure(ivalue)
        % full screen
        set(gcf, 'Position', get(0, 'Screensize'));
        % White background
        set(gcf,'color','white')
        %     clear title
        lgd_str={};
        super_title=[];
        %     title_str= {};

        
        %     set(gcf, 'Position', get(0, 'Screensize'));
        total_floats = size(WMO_vector{ivalue},1);
        
        % Verify if the total floats for the parameter value is above threshold defined in input (default is 10)
        if total_floats < sample_size_limit
            disp('The number of floats for the config parameter value is below the sample size limit threshold')
            continue
        end
        
        super_title= [cellstr(calcul_param) + " = " + unique_values2(ivalue) + " [" + total_floats + " floats total]"];
        annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 18, 'FitBoxToText', 'on', 'Interpreter', 'none')

        subplot(2,2,1)
        %Survival rates per cycles
        %     plot_data_table_cycle2 = plot_data_table_cycle;
        %     plot_data_table_cycle2{5}.PROVOR_II = [];
        %     plot_data_table_cycle2{5}.PROVOR_III = [];
        %     plot_data_table_cycle2{10}.PROVOR_II = [];
        % % %     
        %     unique_model_nb2 = unique_model_nb
        %     unique_model_nb2{5} = unique_model_nb{5} -2
        %     unique_model_nb2{10} = unique_model_nb{10} -1

        for imodel = 1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_table_cycle{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Number of cycles reached')
        ylabel('Survival rate (%)')
%         legend(lgd_str{:}, 'Interpreter', 'none')
        title('Survival rate per cycle per float type')

        subplot(2,2,2)
        %Survival rates per age
        %     plot_data_age2 = plot_data_age;
        %     plot_data_age2{5}.PROVOR_II = [];
        %     plot_data_age2{5}.PROVOR_III = [];
        %     plot_data_age2{10}.PROVOR_II = [];

        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_age{ivalue}{:,1}, plot_data_age{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        xlim([0 max_float_age2{ivalue}])
        ylim([0 100])
        xlabel('Age reached (year)')
        ylabel('Survival rate (%)')
%         legend(lgd_str{:}, 'Interpreter', 'none')
        title('Survival rate per age reached per float type')

        subplot(2,2,3:4)
        %Survival rate per vertical distance traveled
        %     plot_data_vkm2 = plot_data_vkm;
        %     plot_data_vkm2{5}.PROVOR_II = [];
        %     plot_data_vkm2{5}.PROVOR_III = [];
        %     plot_data_vkm2{10}.PROVOR_II = [];
        %     
        for imodel=1:(unique_model_nb{ivalue}-length(vars2delete{ivalue}))
            plot(plot_data_vkm{ivalue}{:,imodel+1}, 'LineWidth', 2)
            lgd_str{imodel} = cell_string{ivalue}{imodel} + " (" + sum(table_model_idx2{ivalue}{:,imodel}) + " floats)";
            hold on
        end
        ylim([0 100])
        xlabel('Vertical distance traveled (Km)')
        ylabel('Survival rate (%)')
        lgd = legend(lgd_str{:}, 'Interpreter', 'none', 'Location', 'eastoutside');
        lgd.FontSize = 11;
        title('Survival rate per vertical Kms traveled per float type')
        
         %%% export
        if is_export==1
            disp('exporting...')        
            filename = [calcul_param{1} + "_" + num2str(unique_values2(ivalue)) + "_surv_rates.jpg"];
            fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/" + "per_value_" + filename];        
            export_fig(fig_path)
        end
    end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if sum(contains(where_file, 'config')) ~= 0 && plot_group_model==1
   for imodel = 1:length(all_models)

    %%% Plot
    figure(imodel)
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % White background
    set(gcf,'color','white')
    %     clear title
    lgd_str={};
    super_title=[];
    %     title_str= {};





    A = cellstr(Floats_list.MODEL);
    for ii = 1:length(A)
        if strcmp(all_models{imodel},A{ii}) == 1
            sum_model{imodel}(ii) = 1;
        else
            sum_model{imodel}(ii) = 0;
        end
    end    
    total_floats_per_model{imodel} = sum(sum_model{imodel});

    super_title= ["Survival rates for " + all_models{imodel} +  " floats, per " + ...
        calcul_param{1} + " parameter values" + " [" + total_floats_per_model{imodel} + " floats total]" ];
    annotation('textbox', [0.5, 1, 0, 0], 'string', super_title, 'HorizontalAlignment', 'center', 'Fontsize', 16, 'FitBoxToText', 'on', 'Interpreter', 'none')

%%%Survival rates per cycles
    subplot(2,2,1)

    for ivalue = 1:length(table_model_idx2)
        if isempty(A_cycles{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
            plot(A_cycles{imodel}{ivalue}, 'LineWidth', 2)
%               lgd_str{ivalue} = "Configuration parameter = " + num2str(unique_values2(ivalue)) + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
            hold on
        else
            continue
% %             plot(0,0,'x')
% %             lgd_str{ivalue} = "No floats with a configuration parameter = " + num2str(unique_values2(ivalue));
%             lgd_str{ivalue} = "";
            hold on
        end
    end
    xlim([0 300])
    ylim([0 100])
    xlabel('Number of cycles reached')
    ylabel('Survival rate (%)')
%     lgd_str = lgd_str(~cellfun('isempty',lgd_str));
%     legend(lgd_str{:})
    title('Survival rate per cycles made')

%%% Survival rates per age
    subplot(2,2,2)

    for ivalue = 1:length(table_model_idx)
        if isempty(A_age{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
            plot(A_age_reached{imodel}{ivalue},A_age{imodel}{ivalue}, 'LineWidth', 2)
%             lgd_str{ivalue} = "Configuration parameter = " + num2str(unique_values2(ivalue)) + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
            hold on
        else
            continue
%             plot(0,0,'x')
%             lgd_str{ivalue} = "No floats with a configuration parameter = " + num2str(unique_values2(ivalue));
            hold on
        end
    end
    xlim([0 (max(cellfun(@max, max_float_age2)))])
    ylim([0 100])
    xlabel('Age reached (year)')
    ylabel('Survival rate (%)')
%     lgd_str = lgd_str(~cellfun('isempty',lgd_str));
%     legend(lgd_str{:})
    title('Survival rate per age reached')

%%% Survival rate per vertical distance traveled
    subplot(2,2,3:4)

    for ivalue = 1:length(table_model_idx)
        if isempty(A_vkm{imodel}{ivalue}) == 0 && size(A_WMO{imodel}{ivalue},1) >= sample_size_limit2 == 1
            plot(A_vkm{imodel}{ivalue}, 'LineWidth', 2)
            lgd_str{ivalue} = "Configuration parameter = " + num2str(unique_values2(ivalue)) + "(" + sum(table_model_idx2{ivalue}.(all_models{imodel})) + " floats)";
            hold on
        else
            continue
%             plot(0,0,'x')
%             lgd_str{ivalue} = "No floats with a configuration parameter = " + num2str(unique_values2(ivalue));
            hold on
        end
    end
    
    ylim([0 100])
    xlabel('Vertical distance traveled (Km)')
    ylabel('Survival rate (%)')
    lgd_str = lgd_str(~cellfun('isempty',lgd_str));
    lgd = legend(lgd_str{:}, 'Interpreter', 'none', 'Location', 'eastoutside');
    lgd.FontSize = 11;
    title('Survival rate per vertical Kms traveled')

%         

    %%% export
    if is_export==1
        disp('exporting...')        
        filename = [all_models{imodel} + "_" + calcul_param{1} + "_surv_rates.jpg"];
        fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/" + "per_model_" + filename];        
        export_fig(fig_path)
    end
   end

end
% end     





%% Histograms plots and statistical repartition

mask_model = cellstr(Floats_list.MODEL);
ARVOR_A_idx = strcmp(mask_model, 'ARVOR_A');
ARVOR_I_idx = strcmp(mask_model, 'ARVOR_I');
ARVOR_L_idx = strcmp(mask_model, 'ARVOR_L');

for ifloats = 1:length(Floats.WMO.data)
    sum_param{ifloats} = sum( Floats.(calcul_param{1}).data{ifloats});
    average_param_per_cycle{ifloats} = mean(Floats.(calcul_param{1}).data{ifloats});
    total_cycles{ifloats} = size(Floats.(calcul_param{1}).data{ifloats}, 1);
end

%%% Figure 4
figure(4)
% full screen
set(gcf, 'Position', get(0, 'Screensize'));
% White background
set(gcf,'color','white')

sum_param2 = cell2mat(sum_param);

%%% Histograms from Arvor-A models
% sum_param2 = sum_param2(ARVOR_A_idx)';

h = histogram(sum_param2');
l = xline(mean(sum_param2),'--','Average');

% xlim([0 7000])
xlabel(['Total of ' calcul_param{1} ' done'], 'Interpreter', 'none')
ylabel('Number of floats')
legend([h,l], ['Sum ' calcul_param{1}], ['Average = ' num2str(round(mean(sum_param2)))], 'Interpreter', 'none')
title(['Statistical repartition of the ' calcul_param{1} ' throughout the whole life of a float'], 'Interpreter', 'none')

%%% export
if is_export==1
    disp('exporting...')        
    filename = ["statistical_repartition_sum_" + calcul_param{1} + "_whole_life.jpg"];
    fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/"  + filename];        
    export_fig(fig_path)
end
    
%%% Figure 5
figure(5)
% full screen
set(gcf, 'Position', get(0, 'Screensize'));
% White background

average_param_per_cycle2 = cell2mat(average_param_per_cycle);

%%% Histograms from Arvor-A models
% average_param_per_cycle2 = average_param_per_cycle2(ARVOR_A_idx)';

h = histogram(average_param_per_cycle2');
l = xline(mean(average_param_per_cycle2),'--','Average','DisplayName',['Average: ' num2str(round(mean(average_param_per_cycle2)))]);

% xlim([0 30])
xlabel(['Average of ' calcul_param{1} ' per cycles'],'Interpreter', 'none')
ylabel('Number of floats')
legend([h,l], ['Sum ' calcul_param{1}], ['Average = ' num2str(mean(average_param_per_cycle2))], 'Interpreter', 'none')
title(['Statistical repartition of the ' calcul_param{1} ' per cycles'], 'Interpreter', 'none')

%%% export
if is_export==1
    disp('exporting...')        
    filename = ["statistical_repartition_average_" + calcul_param{1} + "_per_cycle.jpg"];
    fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/"  + filename];        
    export_fig(fig_path)
end
    
%%% Figure 6
figure(6)
% full screen
set(gcf, 'Position', get(0, 'Screensize'));
% White background
set(gcf,'color','white')

total_cycles2 = cell2mat(total_cycles);

%%% Histograms from Arvor-A models
% total_cycles2 = total_cycles2(ARVOR_A_idx)';

h = histogram(total_cycles2');
l = xline(mean(total_cycles2),'--','Average','DisplayName',['Average: ' num2str(round(mean(total_cycles2))) ' cycles']);

xlim([0 500])
xlabel('Number of cycles reached')
ylabel('Number of floats')
legend([h,l], ['Sum ' calcul_param{1}], ['Average = ' num2str(round(mean(total_cycles2)))], 'Interpreter', 'none')
title('Statistical repartition of the number of cycles reached', 'Interpreter', 'none')

%%% export
if is_export==1
    disp('exporting...')        
    filename = ["statistical_repartition_number_of_cycles.jpg"];
    fig_path = [output_folder + "/" + calcul_param{1} + "/" + working_date + "/"  + filename];        
    export_fig(fig_path)
end






