%{
======================================================================
    Example usage of the code for an offset graft.
    TO DO: Add units.
======================================================================
%}

%% Declare Parameters
contraction_rate = 40;  % Contraction rate of jellyfish (contrations per minute). Parameter measured from experiments, see figure S6 in Supplementary Material.
muscle_strain = 0.25;    % Muscle strain of jellyfish (dimensionless). Parameter measurement from experiments, see figure S6 in Supplementary Materials.      
offset = 3;            % Offset of jellyfish grafts

bulk_modulus = 0.05;    % (!) Units need to be updated.
area = 1.001;           % Initial area of jellyfish.
elast0 = 0.01;            % (Pa) Elasticity of spring. (!) Units need to be updated.
elast1 = 0.04;            % (Pa) Elasticity of spring in Maxwell model. (!) Units need to be updated.
vis = 500;             % (Pa.s) Dashpot in Maxwell model. (!) Units need to be updated.

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
%% Simulation 
% Specify datapath of directory to save data in
datapath = pwd; % Set current directory as datapath

% Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
% Creates folder Data if not already created.
folder_save = ['example_usage_data/', 'offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)];

% Simulate jellyfish
visco_offset_SLM_newmus(elast0, ...
                        elast1, ...
                        vis, ...
                        bulk_modulus, ...
                        area, ...
                        muscle_strain, ...
                        contraction_rate, ...
                        offset, ...
                        folder_save, ...
                        datapath);

%% (Optional) Stitch images into a .avi video
dataDir = fullfile([pwd, '/' ,folder_save]); % Specify location of images
video_name = ['animation_offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)]; 
create_animation(dataDir, video_name); 
