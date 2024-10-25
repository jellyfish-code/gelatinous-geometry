%{
======================================================================
    Example usage of the code for an offset graft and a butterfly graft.
======================================================================
%}
clear; 

%% Declare Parameters
contraction_rate = 20;      % Contraction rate of jellyfish (contrations per minute). Parameter measured from experiments, see figure S6 in Supplementary Material.
muscle_strain = 0.25;       % Muscle strain of jellyfish (dimensionless). Parameter measurement from experiments, see figure S6 in Supplementary Materials.      
offset = 4;                 % Offset of jellyfish grafts

% length_scale = 0.001;       % length scale conversion from mm to m. 
% area_conversion_factor = length_scale^2; 

bulk_modulus = 0.05*1e3;    % Bulk modulus of jellyfish in Pascals. 
area0 = 1.001;              % Initial area of jellyfish in mm^2 (?).
elast0 = 0.05*0.2*1e3;      % Elasticity of spring in Pascals.
elast1 = 0.05*0.8*1e3;      % Elasticity of spring in Pascals.
vis = 500*1e3;              % Viscosity of dashpot in Maxwell model. Units in Pascal*seconds.
damping_coefficient = 500;  % Damping coefficient of jellyfish. Units in Newton*seconds/meter.
max_dR = 1.55;              % Maximum change in radius of jellyfish during contraction.
dR_rate = 0.15;             % Increase in radius change with distance from anchored end (Figure S6c in paper).


%% Parameters to predict whether shape will reorganise, see Chapter 4.1 in M. Gong's thesis (2022). 
characteristic_ratio = muscle_strain*(elast0 + elast1)/vis;     % Should be betweeen 0.4e-6 and 0.3e-5 (1/s) for full reorganisation. 
relaxation_time = (vis/elast1)/3600;                            % Relaxation time greater than 11.6 hours did not fully reorganise.
elasticity_ratio = elast0/elast1;                               % 1:4 ratio resulted in full reorganisation, whereas 1:1 ratio did not.

% Assert correct ranges 
if characteristic_ratio < 4e-7 || characteristic_ratio > 3e-6
    warning("Characteristic ratio not in range for full reorganisation."); 
elseif relaxation_time > 11.6
    warning("Relaxation time greater than 11.6 hours does not lead to reorganise."); 
elseif elasticity_ratio == 1
    warning("Elasticity ratio of 1:1 does not lead to reorganisation.")
end


%% Uncomment to simulate an offset graft
% % Specify datapath of directory to save data in
% datapath = pwd; % Set current directory as datapath
% 
% % Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. Creates folder example_usage_data if not already created.
% folder_save = ['example_usage_data/', 'offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate), '_using_visco_offset_SLM_newmus_remote_version_updated_files_corrected_conversion_factors_nonessential_files_deleted'];
% 
% visco_offset_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, offset, folder_save, datapath)

%% Uncomment to simulate a butterfly graft
% Specify datapath of directory to save data in
datapath = pwd; % Set current directory as datapath

% Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
% Creates folder Data if not already created.
folder_save = ['example_usage_data/', 'butterfly_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate), '_using_visco_offset_SLM_newmus_remote_version_updated_files_corrected_conversion_factors_nonessential_files_deleted'];

visco_butterfly_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, folder_save, datapath)

%% (Optional) Stitch images into a .avi video
dataDir = fullfile([pwd, '/' ,folder_save]); % Specify location of images
video_name = ['animation_offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_damping_coefficient_', num2str(damping_coefficient), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)]; 
create_animation(dataDir, video_name); 