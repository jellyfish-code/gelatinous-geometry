%{
======================================================================
    Runs parallelized parameter sweep. 
======================================================================
%}

%% Declare Non-Sweep Parameters
muscle_strain = 0.25;       % Muscle strain of jellyfish (dimensionless). Parameter measurement from experiments, see figure S6 in Supplementary Materials.      

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

%% Declare Sweep Parameters
contraction_rate_sweep = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]; 
graft_diameter = 10;
offset_to_diameter_ratio_sweep = [0.1, 0.2, 0.3, 0.4, 0.5];
offset_sweep = offset_to_diameter_ratio_sweep*10;
parameter_sweep_table = combinations(contraction_rate_sweep, offset_sweep); 

%% Start a parallel process

for i = 1:height(parameter_sweep_table) 
    contraction_rate = parameter_sweep_table(i,:).contraction_rate_sweep; % Contraction rate of jellyfish (contrations per minute.
    offset = parameter_sweep_table(i,:).offset_sweep;                 % Offset of graft in Pascals.

    % Specify datapath of directory to save data in
    datapath = pwd; % Set current directory as datapath
    graft_type = '_offset_graft'; 

    % Specify subfolder to save data in. Data is saved in datapath (here, the current directory) inside folder Data. 
    % Creates folder Data if not already created.
    date = sprintf('%s', datetime("today"));
    folder_save = ['HPC_runs/parameter_sweep_figure5e/', date, graft_type, '_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)];

    % visco_offset_SLM_newmus(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, offset, folder_save, datapath)
    
    %% Post Processing
    
    % % List number of images
    % 
    % % Get the list of all files and folders in the directory
    % fileList = dir([folder_save, '/graft_reorganization']);
    % 
    % % Exclude directories ('.' and '..' and any other subdirectories)
    % fileList = fileList(~[fileList.isdir]);
    % 
    % % Count the number of files
    % numFiles = numel(fileList);
    % 
    % % Display the result
    % fprintf('Number of files in the directory: %d\n', numFiles);

    % Generate animation
    dataDir = fullfile([pwd, '/' ,folder_save, '/graft_reorganization']); % Specify location of images
    video_name = [date, graft_type, '_elast0_', num2str(elast0), '_elast1_', num2str(elast1), '_viscosity_', num2str(vis), '_bulk_modulus_', num2str(bulk_modulus), '_offset_', num2str(offset), '_contraction_rate_', num2str(contraction_rate)]; 
    create_animation(dataDir, video_name); 
end