%{
======================================================================
    Example usage of the code for an offset graft.
    TO DO: Add units.
======================================================================
%}

%% Declare Parameters
contraction_rate = 20;  % Contraction rate of jellyfish (contrations per minute). Parameter measured from experiments, see figure S6 in Supplementary Material.
muscle_strain = 0.2;    % Muscle strain of jellyfish (dimensionless). Parameter measurement from experiments, see figure S6 in Supplementary Materials.      
offset = 2;            % Offset of jellyfish grafts

bulk_modulus = 0.1;
area = 1.001;           % Initial area of jellyfish.
elast0 = 0.05;          % Elasticity of spring. 
elast1 = 0.05;          % Elasticity of spring in Maxwell model. 
vis = 1000;              % Dashpot in Maxwell model. 

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
