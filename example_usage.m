%{
======================================================================
    Example usage of the code for an offset graft.
    TO DO: Add units.
======================================================================
%}

%% Declare Parameters
contraction_rate = 20;  % Contraction rate of jellyfish (contrations per minute). Parameter measured from experiments, see figure S6 in Supplementary Material.
muscle_strain = 0.21;    % Muscle strain of jellyfish (dimensionless). Parameter measurement from experiments, see figure S6 in Supplementary Materials.      
offset = 2;            % Offset of jellyfish grafts

bulk_modulus = 2e6;      
area = 1.001;           % Initial area of jellyfish.
elast0 = 5;            % (Pa) Elasticity of spring. 
elast1 = 15;            % (Pa) Elasticity of spring in Maxwell model. 
vis = 1e5;             % (Pa.s) Dashpot in Maxwell model. 

%% Parameters to predict whether shape will reorganise, see Chapter 4.1 in M. Gong's thesis (2022). 
characteristic_ratio = muscle_strain*(elast0 + elast1)/vis;     % Should be betweeen 0.4e-6 and 0.3e-5 (1/s) for full reorganisation. 
relaxation_time = (vis/elast1)/3600;                            % Relaxation time greater than 11.6 hours did not fully reorganise.
elasticity_ratio = elast0/elast1;                               % 1:4 ratio resulted in full reorganisation, whereas 1:1 ratio did not.
%% Simulation 

% Specify datapath of directory to save data in
datapath = pwd; % Set current directory as datapath

% Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
% Creates folder Data if not already created.
folder_save = ['example_usage_data/', 'offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus)];

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
create_animation(dataDir); 
