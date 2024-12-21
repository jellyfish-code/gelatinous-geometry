% To test slurm job set-up

disp("Successful Run.")

%% Writing images
path0 = pwd();
folder_save = 'test_slurm_script'; 
Dr = dir([path0 '/' folder_save]);
disp("directory"); 
disp(Dr);
S = mkdir(path0,folder_save);  
