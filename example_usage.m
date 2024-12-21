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

bulk_modulus = 0.05*1e3;    % Bulk modulus of jellyfish in Pascals. 
area0 = 1.001;              % Relaxed area as a percentage of jellyfish area. Greater than 1.
elast0 = 0.05*0.2*1e3;      % Elasticity of spring in Pascals.
elast1 = 0.05*0.8*1e3;      % Elasticity of spring in Pascals.
vis = 500*1e3;              % Viscosity of dashpot in Maxwell model. Units in Pascal*seconds.
damping_coefficient = 500;  % Damping coefficient of jellyfish. Units in Newton*seconds/meter.
max_dR = 1.55;              % Maximum change in radius of jellyfish during contraction.
dR_rate = 0.15;             % Increase in radius change with distance from anchored end (Figure S6c in paper).

tau_SLM = vis*(elast0 + elast1)/(elast0*elast1); 
tau_maxwell = vis/elast1; 
tau = min(tau_SLM, tau_maxwell); 
timestep_fraction_of_tau = tau/(30*60); % current timestep is 30 mins

%% Uncomment to simulate an offset graft
% % Specify datapath of directory to save data in
% datapath = pwd; % Set current directory as datapath
% graft_type = '_offset_graft'; 
% % Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. Creates folder example_usage_data if not already created.
% date = string(datetime("today")); 
% folder_save = ['example_usage_data/', date, graft_type, '_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)];
% 
% visco_offset_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, offset, folder_save, datapath)

%% Uncomment to simulate a butterfly graft
% % Specify datapath of directory to save data in
% datapath = pwd; % Set current directory as datapath
% graft_type = '_butterfly_graft'; 
% % Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
% % Creates folder Data if not already created.
% date = string(datetime("today")); 
% folder_save = ['example_usage_data/', date, graft_type, '_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)];
% 
% visco_butterfly_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, folder_save, datapath)

%% (Optional) Stitch images into a .avi video
% dataDir = fullfile([pwd, '/' ,folder_save]); % Specify location of images
% video_name = ['animation', graft_type, 'corrected_numerical_algorithm_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_damping_coefficient_', num2str(damping_coefficient), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)]; 
% create_animation(dataDir, video_name); 

%% Uncomment to start a parallel process
elasticity_range = [0.05*1e3, 1e3]; 

parfor i = 1:2 
    elast0 = elasticity_range(i)*0.2;      % Elasticity of spring in Pascals.
    elast1 = elasticity_range(i)*0.8;      % Elasticity of spring in Pascals.

    % Specify datapath of directory to save data in
    datapath = pwd; % Set current directory as datapath
    graft_type = '_butterfly_graft'; 
    
    % Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
    % Creates folder Data if not already created.
    date = sprintf('%s', datetime("today"));
    folder_save = ['example_usage_data/', date, graft_type, '_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)];

    visco_offset_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, offset, folder_save, datapath)
end