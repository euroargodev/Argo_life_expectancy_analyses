# EA_RISE_WP2.1
Repository for the scripts created during the Work Package 2.1 of the Euro-Argo RISE project
This repository will be spearated in 3 different axes:
- DMQC
- Life expectancy
- Toolbox


      A.	Life expectancy related scripts
1.	Plot fleet configuration values (CONFIG_fleet_status.m)
Script description:
% plots the status of a list of floats regarding a configuration parameter and splitting the results depending on the country, deployment date and float model
Output:
Figure config param changed? / Figure config param value for floats which did not change configuration / Figure config param value per cycle. 

 
2.	Mapping of the technical parameters
Description of the script:
Plot a scatter of a tech, traj or config parameter in a map (lon, lat) using a threshold.
•	colormap_thres: if 'above', colormap is use for points above threshold. If 'under', colormap is use for points under threshold.
•	grounded: if 'yes', it plots grounded cycles, if 'no', it does not plot grounded cycles. It only works for ARVOR floats (since the grounding flags on Apex floats are only recorded since the APF11).
NOTES:
(1) Depending on the number of floats, the script will take some minutes
(2) If there is only one float in input list, trajectory lines are plotted.
(3) Be careful with RAM memory. If you use so many floats the program can run out of memory.
 
 

3.	Plot configuration parameters survival rates (survival_config.m)
Pb ligne 394 sur la séparation par modèles !! A regarder de plus près
Script description:
Plots survival rates depending on config values for a given config parameter (given in number of cycles, vertical_km and float age). This script takes as an input a WMOs list.
NOTES: 
(1) Using as input config_param = {'CONFIG_CTDPoints_NUMBER'} the script calculates the number of theoretical CTD points per profile using the function calculate_CTDPoints
(2) Multiple red messages by come up in the prompt windows saying that the variable was not found. It is normal, the variable in question is a cycle number index which might be named differently for another float type than ARVOR. 
(3) Possibility to fix a sample size limit (min recommended=10 floats). Computing a survival rate for a smaller sample than that would not be very reliable as one float would could have a huge impact on the resulting survival rate curve and not represent a general trend of the sample.
Output: 
Figures with survival rates depending on number of cycles, vertical km and float age 


4.	Plot survival rates (Plot_survival_rates.m)
Description of the script:
This script permits to plot different survival rates computations stored as variables (.mat). It permits comparison between different floats samples, models, etc… The survival rates of the floats samples, differentiation between models, etc. is done in an auxiliary script (to rename): /home1/datahome/co_arg/larduini/Andrea/Life_expectancy/plot_survival_rate.m
The survival rates are computed according to different x-axis: number of cycles, float age and vertical distance traveled (in km)
You will find below two examples of the plots produced by this script: one according to the area of deployment (1) and the other, every float with one precise config parameter and differentiated by floats models (2):

 
5.	Groundings repartition (Map_groundings.m)
The script was developed by Luca Arduini Plaisant in order to provides some discussing materials for the WP6 (Marginal Seas) of the EA RISE project. It was then used in the WP2.1 report to represent the impact of the groundings on the life expectancy of a float sample and to help quantify it.
Description of the script:
% This script aims to plot on a map the different types of groundings flags per cycle according to a WMO list furnished in the input.
% Groundings flags are stored in the traj files of a float. It can take either: 
% Y= Cycle DID ground; N= Cycle did NOT ground; U= Unknown (corresponds to floats that do not store this information. Apex floats before APF11 Version)
% You can either put a WMOs list as an input or just a single WMO number (by commenting the list part)
NOTES:
(1)	One of the “difficulties” of the script was to convert a variable from the N-cycle dimension to the N-measurements one. The cycle number and latitude/longitude variables are following the N-measurements dimension and the grounding flags are following the N-cycle dimension. The variable “Cycle_nb_INDEX” stored in the traj file permits to link the two dimensions, following a few steps presented throughout the script.

This script produces different types of outputs:
-	A statistical repartition of the grounding’s flag for the WMOs list considered
-	A geographical repartition of these flags (see below)
-	A geographical repartition of the flags per WMOs to permit a specific float tracking for example.
Example of a geographical repartition output from this script:

      B. DMQC related scripts
1.	DMQC statistics (get_DMQC_stats.m)
Description of the script:
Gets DMQC statistics for given floats.
Notes:
(1)	Reading index file may take 2 or 3 minutes
(2)	It is better to use only ascent profiles (Dprof = 0) for coherence with figures from “get_DMQC_adjustment.m” script
    
   
 
2.	Salinity adjustments for a list of floats (get_DMQC_adjustment.m)
Description of the script:
Shows psal adjustment for a list of floats
Notes: 
- (1) Getting data process may take some minutes (depending on number of floats) because we are checking all profile files
- (2) Descent profiles are not included
- (3) Format options (number of floats per figure and yaxis ticks size) can be modified
 
 
      C. TOOLBOX
    
Name	Description	Auxiliary functions
Scripts
map_tech_param:	plot a scatter of a tech, traj or config parameter in a map (lon, lat) using a threshold	read_csv
get_floats_data_gdac
M_MAP: matlab package
export_fig: package
get_csv_QGIS
get_phases_times:	calculates duration and speed of floats in each cycle phase (descent to park, parking, descent to profile, profile drift and ascent to surface)	suptitle
get_floats_data_gdac
get_floats_filespath
export_fig
get_last_config_values:	get list of last configuration parameters values	read_csv
get_floats_data_gdac

Important functions
read_csv:	read a csv file and generates an struct with file variables	
get_floats_files_paths:	gets files paths from ar_index_global_meta.txt file for given floats and optionally creates a .txt file	
get_data_from_index:	gets chosen variables from index file for given floats	
get_floats_data_gdac:	gets data from tech, traj, aux, meta or index files	get_traj_param
get_tech_param
get_configparam_meta
get_data_from_index
get_floats_data_gadc_format:	gets data from tech, traj, aux or meta file and formats it to have the same number of cycles for the same float	get_traj_param
get_tech_param
get_configparam_meta
get_data_from_index
get_configparam_meta:	gets all cycles values of a list of given configuration parameters from a meta file (one float)	
get_tech_param:	gets all cycles values of a list of given technical parameters from a tech file (one float)	
get_traj_param:	gets all values of a list of given parameters from a trajectory file	
get_traj_param_AUX		

