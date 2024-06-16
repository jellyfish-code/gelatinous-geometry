%{
======================================================================
    Example usage of the code for an offset graft.
    TO DO: Add units.
======================================================================
%}

%% Declare Parameters
contraction_rate = 20;  % Contraction rate of jellyfish
offset = 2;             % Offset of jellyfish grafts
muscle = 0.2;          

bulk_modulus = 0.1;
area = 1.001;           % Initial area of jellyfish.
elast0 = 0.05;          % Elasticity of spring. 
elast1 = 0.05;          % Elasticity of spring in Maxwell model. 
vis = 1000;             % Dashpot in Maxwell model. 

%% Simulation 

% Specify datapath of directory to save data in
datapath = pwd; % Set current directory as datapath

% Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
% Creates folder Data if not already created.
folder_save = ['example_usage_data/', 'offset_graft_', 'elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_vis_', num2str(vis), '_bulk_mod_', num2str(bulk_modulus)];

% Simulate jellyfish
visco_offset_SLM_newmus(elast0, ...
                        elast1, ...
                        vis, ...
                        bulk, ...
                        area_range, ...
                        muscle_strain, ...
                        contraction_rate, ...
                        offset, ...
                        folder_save, ...
                        datapath);