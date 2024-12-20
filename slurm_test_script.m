% To test slurm job set-up

print("Successful Run.")

%% Writing images
path0 = pwd();
folder_save = 'test_slurm_script'; 
Dr = dir([path0 '/' folder_save]);
S = mkdir(path0,folder_save);  