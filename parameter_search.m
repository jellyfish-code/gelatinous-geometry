%{
======================================================================
    Master file that runs model with different parameter sets to find viscous/elastic moduli that lead to various solutions
======================================================================
%}

%% Declare Parameters
% A0 going from 10^3 to 10^6 (collagen gels from 5x10^2 to 2x10^4)
%Viscosity from 10^4 - 10^7
%Bulk modulus rubber is around 2x10^6 Pa, water at 2x10^9 Pa

%contraction_range = [1, 20, 40, 60];
contraction_rate = 20;
offset = 2;
muscle_range = [0.05 0.1 0.15 0.2 0.25 0.3 0.35];

bulk_range = [0.05 0.1 0.5];
area_range = 1.001; %, 1.005, 1.01, 1.05, 1.1];
elast_range = [0.025 0.05 0.075 0.1 0.15 0.2 0.25];
vis_range = [100 500 1000 5000 10000];

params = [0.025 1; 0.025 1.5; 0.05 0.75; 0.05 1; 0.075 0.5; 0.075 0.75; 0.1 0.5; 0.1 0.75; 0.15 0.5; 0.2 0.5; 0.25 0.25; 0.25 0.5];

% Specify datapath of directory to save data in
datapath = './Data'; % Saves data in current directory inside folder Data. Creates folder Data if not already created.

%% Simulate Processes

% Execute for-loop iterations in parallel processes
parfor a = 1:length(elast_range) %bleh = 1:length(params_left) %1:length(params) %a = 
    %a = params_left(bleh);
    %a = params_left(x);
    %elast = params(a, 1);
    elast = elast_range(a);
    elast0 = elast*0.5;
    elast1 = elast*0.5;
    %muscle_strain = params(a, 2);
    for b = 4:length(vis_range)
        %vis = params(a,2);
        vis = vis_range(b);
        for c = 1:length(bulk_range) 
            %bulk = params(a,3);
            bulk = bulk_range(c);
            %bulk = vis/1000;
            for d = 1:length(muscle_range)
                muscle_strain = muscle_range(d);
                %folder_save = ['realtime_021221', num2str(a), '_', num2str(b), num2str(c), num2str(d)]; %, 
                %realtime_contraction_SLM(elast0, elast1, vis, bulk, muscle_strain, folder_save);

                % Specify subfolder to save data in
                folder_save1 = ['SLM_off_031921', num2str(a), num2str(b), num2str(c), num2str(d)]; %
                
                % Simulate jellyfish
                visco_offset_SLM_newmus(elast0, ...
                                        elast1, ...
                                        vis, ...
                                        bulk, ...
                                        area_range, ...
                                        muscle_strain, ...
                                        contraction_rate, ...
                                        offset, ...
                                        folder_save1, ...
                                        datapath);
                %folder_save2 = ['SLM_but_031921', num2str(a), num2str(b), num2str(c), num2str(d)]; %
                %visco_butterfly_SLM_newmus(elast0, elast1, vis, bulk, area_range, muscle_strain, contraction_rate, folder_save2);
            end
        end
    end
end
%          %% A0 going from 10^3 to 10^6 (collagen gels from 5x10^2 to 2x10^4)
% %Viscosity from 10^4 - 10^7
% %Bulk modulus rubber is around 2x10^6 Pa, water at 2x10^9 Pa
% 
% contraction_range = [1, 20, 40, 60];
% %contraction_rate = 20;
% offset_range = [1, 2, 3];
% muscle_range = [0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8]; %[2.0 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6];
% 
% bulk_range = [8];
% area_range = 1.001; %, 1.005, 1.01, 1.05, 1.1];
% elast_range = [0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2]; %[0.02 0.04 0.06 0.08 0.1 0.15 0.2 0.25 0.3];
% vis_range = [7500];
% 
% %params = [0.15, 8100, 9.5, 1.35; 0.2 8600  10 1.35; 0.1 8500 10 1.4; 0.2 7900 9 1.5; 0.225 7800 9 1.5];
% 
% parfor a = 1:length(elast_range) %a = 1:length(params) %a = 1:length(elast_range) %
%     %a = params_left(x);
%     %elast = params(a, 1);
%     elast = elast_range(a);
%     elast0 = elast*0.5;
%     elast1 = elast*0.5;
%     for b = 1:length(vis_range)
%         %vis = params(a, 2);
%         vis = vis_range(b);
%         for c = 1:length(bulk_range) 
%             %bulk = params(a,3);
%             bulk = bulk_range(c);
%             for d = 1:length(muscle_range)
%                 %muscle_strain = params(a,4);
%                 %for b = 1:length(offset_range)
%                 %offset = offset_range(b);
%                 %for c = 1:length(contraction_range)
%                 %contraction_rate = contraction_range(c);
%                 muscle_strain = muscle_range(d);
%                 folder_save = ['realtime_011821', '1', num2str(a), '_', num2str(d)]; %, num2str(b) num2str(c), 
%                 realtime_contraction_SLM(elast0, elast1, vis, bulk, muscle_strain, folder_save);
%                 %folder_save1 = ['SLM_off_111920', 'params', num2str(a), 'off', num2str(c), 'rate', num2str(b)]; %num2str(a), num2str(b), num2str(c), num2str(d)]; %
%                 %visco_offset_SLM(elast0, elast1, vis, bulk, area_range, muscle_strain, contraction_rate, offset, folder_save1);
%                 %folder_save2 = ['SLM_but_111520_run3', num2str(a), num2str(b), num2str(c), num2str(d)]; %'params', num2str(a)]; %
%                 %visco_butterfly_SLM(elast0, elast1, vis, bulk, area_range, muscle_strain, contraction_rate, folder_save2);
%             end
%         end
%     end
% end
% % end