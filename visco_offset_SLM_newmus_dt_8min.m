%{
======================================================================
    Main function evolving jellyfish offset mesh.  
======================================================================

INPUT:
        elast0 (double):                 Elasticity of spring (in Pascals).
        elast1 (double):                 Elasticity of spring (in Pascals).
        vis (double):                    Viscosity of dashpot (in Pascal*seconds).
        damping_coefficient (double):    Damping coefficient of jellyfish (n Newton*seconds/meter).
        bulk_modulus (double):           Bulk modulus of jellyfish (in Pascals).
        area0 (double):                  Relaxed area as a percentage of jellyfish area. Greater than 1.
        muscle_strain (double):          Strain of jellyfish in ? (dimensionless quantity).
        contraction_rate (double):       Number of jellyfish contractions per minute.
        max_dR (double):                 Maximum change in radius during contraction.
        dR_rate (double):                Increase in radius change with distance from anchored end (Figure S6c in paper).
        folder_save (string):            Name of folder in which simulation images is to be saved.
        datapath (string):               Directory in which folder specificed in folder_save can be found.

OUTPUT:
        None. 
%}

function visco_offset_SLM_newmus_dt_8min(elast0, elast1, vis, damping_coefficient, bulk_modulus, area0, muscle_strain, contraction_rate, max_dR, dR_rate, offset, folder_save, datapath)
    %% Set up parameters, everything in Pa(N/m^2) and s

    %Measured parameters
    contraction_duration = 0.8; %s
    contraction_strength = (elast0+elast1)*muscle_strain; %Pa
    relax_duration = (60-contraction_rate*contraction_duration)/(contraction_rate + 1); %seconds
    
    if contraction_rate > 75
        contraction_duration = 60/contraction_rate;
    end
    if relax_duration < 0
        relax_duration = 0;
    end

    %time
    time_step = 8; %minutes
    time_end = 2000; %hours
    time_steps = time_end*60/time_step;
    a_r = [];
    vel = [];

    % %% Calculate the Maxwell relaxation constants. Used to update relaxed length of spring 1 in series with dashpot.
    % relax_param = (1-exp(-1*elast1/vis * (time_step * 60)));

    %% Writing images
    path0 = datapath;

    if ~isempty(path0)
        Dr = dir([path0 '/' folder_save '/graft_reorganization']);
        if isempty(Dr)
            S = mkdir(path0,[folder_save, '/graft_reorganization']);
            if ~S, disp('Fail to make folder!'); return;  
            end
        else
        end                        
        path1 =[path0, '/', folder_save];         % the place to store images
    else
        path0 = pwd; 
    end


    %% Set up jellyfish geometry as an array
    %This creates the initial offset array, assigning coordinates for each
    %node and connecting nodes. It also designates muscle and non-muscle
    %nodes.

    %% Start with uncut jellyfish, find the balance of elastic and pressure forces
    %%This function assumes that maxwell completely relaxes, i.e. the Kelvin-Voigt model
    [jelly_eq, ~, d_uncut] = equilibrium_initial_KV(elast0, vis, damping_coefficient, bulk_modulus, 100, 1200, area0, contraction_strength);

    %% Make a matrix with the right shape for offset graft
    [jelly_initial, row_start, row_end] = offset_mesh(offset); 
    [muscle_outer, muscle_inner] = offset_muscle(offset);
    a = size(jelly_initial);

    %% Then cut the jellyfish (update the x_ and y_coordinates) 
    jelly_initial(6:a(1), 1: a(1), :) = jelly_eq(6:a(1), 1: a(1), :);
    jelly_initial(1:6, 1+offset: a(1)+offset, 1) = jelly_eq(1:6, 1: a(1), 1) + offset*d_uncut;
    jelly_initial(1:6, 1+offset: a(1)+offset, 2) = jelly_eq(1:6, 1: a(1), 2);
    jelly_initial(6, 1+offset:6+offset, 1) = jelly_eq(6, 1:6, 1) + offset*d_uncut;
    jelly_initial(6, 7:a(1), 1) = jelly_eq(6, 7:a(1), 1);

    %% Convert array to graph
    jelly = convert_jelly_graph(jelly_initial, row_start, row_end);

    %% Initialise velocity, net force, elastic stress, pressure, muscle stress, outer muscles, and inner muscles in jelly graph. 
    jelly.Nodes.velocity = zeros(numnodes(jelly), 2);
    jelly.Nodes.F_net = zeros(numnodes(jelly), 2);
    jelly.Nodes.stress_elastic = jelly.Nodes.F_net;
    jelly.Nodes.pressure = jelly.Nodes.F_net;
    jelly.Nodes.stress_muscle = jelly.Nodes.F_net;
    jelly.Nodes.outmus = zeros(numnodes(jelly),1);
    jelly.Nodes.inmus = zeros(numnodes(jelly),1);

    %% Define nodes that make up the edges (in order)
    edges = find_edges(muscle_outer);
    j_edges = zeros(numnodes(jelly),1);
    c = 1;
    for i = 1:length(edges)-1
        for j = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(j,1) == edges(1, i) && jelly.Nodes.Node_Names(j,2) == edges(2, i)
                j_edges(j) = c;
                c = c+1;
            end
        end
    end

    jelly.Nodes.edges = j_edges;
    m1 = zeros(length(muscle_outer),1);
    m2 = zeros(length(muscle_outer),1);
    m3 = zeros(length(muscle_inner),1);
    m4 = zeros(length(muscle_inner),1);

    %%More parameters
    [jelly_i, row_start, row_end] = offset_mesh(offset); 
    jelly_off_i = convert_jelly_graph(jelly_i, row_start, row_end);
    jelly_off_i.Nodes.edges = j_edges;
    [j_area, ~] = area(jelly_off_i);
    area_relax = area0*j_area;

    %% Initialise relaxed lengths of springs and calculate initial strains
    jelly.Edges.d_rel0 = jelly_off_i.Edges.d_current; %This is a weird one. The relaxed length is the length
    % %before equilibrium is found. So I'm just initializing another offset graft

    % % IF IGNORING KV INITIALISATION, test initializing the simulation with zero stresses 
    % jelly.Edges.d_current = jelly_off_i.Edges.d_current; 
    % jelly.Edges.d_rel0 = jelly.Edges.d_current; % No stress across spring 0.

    jelly.Edges.d_rel1 = jelly.Edges.d_current; %This is assumed to be fully relaxed
    jelly.Edges.strain0 = (jelly.Edges.d_current - jelly.Edges.d_rel0)./jelly.Edges.d_rel0;
    jelly.Edges.strain1 = (jelly.Edges.d_current - jelly.Edges.d_rel1)./jelly.Edges.d_rel1;
    jelly.Edges.strainviscous = jelly.Edges.strain0 - jelly.Edges.strain1; % If spring 1 is assumed to be completely relaxed, then strain1 is zero, and all of the strain in the maxwell arm is across the viscous dashpot.
    
    outcount = 1;
    incount = 1;
    %% Define edges that make up the muscles
    for j = 1:length(muscle_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_outer(1,j,1) && jelly.Nodes.Node_Names(i,2) == muscle_outer(2,j,1)
                m1(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end
    for j = 1:length(muscle_outer)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_outer(1,j,2) && jelly.Nodes.Node_Names(i,2) == muscle_outer(2,j,2)
                m2(j) = i;
                jelly.Nodes.outmus(i) = outcount;
                outcount = outcount+1;
            end
        end
    end
    for j = 1:length(muscle_outer)-1
        a = findedge(jelly, m1(j), m1(j+1));
        jelly.Edges.muscle(a) = 1;
        b = findedge(jelly, m2(j), m2(j+1));
        jelly.Edges.muscle(b) = 1;
    end
    for j = 1:length(muscle_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_inner(1,j,1) && jelly.Nodes.Node_Names(i,2) == muscle_inner(2,j,1)
                m3(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    for j = 1:length(muscle_inner)
        for i = 1:numnodes(jelly)
            if jelly.Nodes.Node_Names(i,1) == muscle_inner(1,j,2) && jelly.Nodes.Node_Names(i,2) == muscle_inner(2,j,2)
                m4(j) = i;
                jelly.Nodes.inmus(i) = incount;
                incount = incount+1;
            end
        end
    end
    for j = 1:length(muscle_inner)-1
        a = findedge(jelly, m3(j), m3(j+1));
        jelly.Edges.muscle(a) = 0.5;
        b = findedge(jelly, m4(j), m4(j+1));
        jelly.Edges.muscle(b) = 0.5;
    end

    muscle_length = sum(jelly.Edges.d_current(jelly.Edges.muscle == 1));
    [jelly, ~] = remesh_SLM_newmus(jelly, muscle_length);
            
    if max(jelly.Edges.d_current) > 10
        display('Edge length greater than 10mm, Line 194');
        return
    elseif any(isfinite(jelly.Nodes.x_coord)-1) == 1 || any(isfinite(jelly.Nodes.y_coord)-1) == 1
        display('Infinity in coordinates, Line 197');
        return
    end

    %% Image initial jellyfish graft
    figure_graft = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
    hold off
    xlim([0, 11 + offset]);
    ylim([0, 11 + offset]);
    axis square; 

    path_graft = [path1, '/graft_reorganization'];
    cd(path_graft);                                                            % write the image data
    saveas(figure_graft, '0.jpg')
    cd(path0);     
    pause(0.001)


  %========================================================================================================================
    %% Start the simulation
    for time = 1:time_steps

       hours = (time/(60/time_step));

        %% Update the current length of edges
        for i = 1:numedges(jelly)
            [node1, node2] = findedge(jelly, i);
            dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
            dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
            dist_current = (dx^2 + dy^2)^(1/2);
            jelly.Edges.d_current(i) = dist_current;
        end
        
        if max(jelly.Edges.d_current) > 10 || max(jelly.Edges.strain0) > 1 || min(jelly.Edges.strain0) < -0.5
            cd(path1)
            writematrix(a_r, 'a_r.xlsx');
            writematrix(vel, 'velocity.xlsx');
            cd(path0)
         
            return
        end
            
        %% Calculate on each node stress from muscle contraction
        % Call contraction5offset.m calculates coordinate and direction of contraction
        [jelly, done] = contraction5offset(jelly, contraction_strength, muscle_strain, max_dR, dR_rate);
        if done == 0 % If an unstability or error arose in contraction5offset, exit simulation.
            cd(path1)
            a_r = cat(1, a_r, [1, hours]);
            writematrix(a_r, 'a_r.xlsx');
            writematrix(vel, 'velocity.xlsx');
            cd(path0)
            return
        end
        
        %% Estimate elastic stress arising from springs for each node 
	% and direction of elastic force 
	% Elastic force is function of the dashpot viscosity  
 	% jelly = SLM_elastic(jelly, elast0, elast1);
        jelly = SLM_viscoelastic(jelly, elast0, elast1, vis, time_step); 
        
	%% Estimate stress from pressure for each node 
 	% and direction of pressure force   
        jelly.Nodes.pressure = find_f_pressure(jelly, area_relax, bulk_modulus);  

	%% Estimate total stress for each node
	% Vectorial addition using the direction information   
 	% Muscle stress is present only during contraction. Implemented using contraction duration, which is the fraction of time that jellyfish spends on average contracted. 
 	stress_contract = jelly.Nodes.stress_elastic + jelly.Nodes.pressure + jelly.Nodes.stress_muscle;
	stress_relax = jelly.Nodes.stress_elastic + jelly.Nodes.pressure;
        stress_net = (stress_contract*contraction_rate*contraction_duration + stress_relax*relax_duration*(contraction_rate+1))/60;

	%% Estimate force from stress  
	edge_crossectional_area = 1e-3*1e-3;                    % Characteristic cross-ectional area, use initial mesh size. Units in meters squared.
        jelly.Nodes.F_net = stress_net*edge_crossectional_area; % Force in Newtons.
	

        %% Force balance equation 
	% Compute displacement each node. Vis Pa*s = Ns/m2, damping Ns/m
        jelly.Nodes.velocity = 1e3*jelly.Nodes.F_net./damping_coefficient; % Converts meters per second to millimeters per second.
        contraction_displacement = jelly.Nodes.velocity*time_step*60;
        jelly.Nodes.x_coord = jelly.Nodes.x_coord + contraction_displacement(:,1);
        jelly.Nodes.y_coord = jelly.Nodes.y_coord + contraction_displacement(:,2);

        % Update relaxed length of spring 1 due to the presence of dashpot
        % in series. This is done for the remeshing code.
        % jelly.Edges.d_rel1 = -1*(jelly.Edges.d_current.*jelly.Edges.d_rel1)./((jelly.Edges.d_rel1 - jelly.Edges.d_current).*relax_param - jelly.Edges.d_rel1);
       
        
        %% Remesh every 5 hours OR if constraints are met
        %% TODO: remove requirements of every 5 hours?
        [trash, idx] = sortrows(jelly.Nodes, {'edges'});
        edge_idx = idx(trash.edges~=0);
        new_con = [];
        for i = 1:length(edge_idx)
            node1 = edge_idx(i);
            node2 = edge_idx(mod(i+1, length(edge_idx))+1);

            dx = jelly.Nodes.x_coord(node1) - jelly.Nodes.x_coord(node2);
            dy = jelly.Nodes.y_coord(node1) - jelly.Nodes.y_coord(node2);
            d = (dx^2 + dy^2)^(1/2);

            new_con = cat(1, new_con, d);
        end
        
        if mod(hours, 5) == 0 || min(jelly.Edges.d_current) < 0.35 || min(new_con) < 1
            [jelly, lim_reached] = remesh_SLM_newmus(jelly, muscle_length);
            if lim_reached == 1
                display('Remeshing Limit Reached, Line 292'); 
                cd(path1)
                a_r = cat(1, a_r, [2, hours]);
                writematrix(a_r, 'a_r.xlsx');
                writematrix(vel, 'velocity.xlsx');
                cd(path0)
                return
            end
        end
        
   

        %% Find the current aspect ratio
        if any(isfinite(jelly.Nodes.x_coord)-1) == 1 || any(isfinite(jelly.Nodes.y_coord)-1) == 1
            cd(path1)
            a_r = cat(1, a_r, [3, hours]);
            writematrix(a_r, 'a_r.xlsx');
            writematrix(vel, 'velocity.xlsx');
            cd(path0)
            display('Infinity in coordinates, Line 311'); 
            return
        end

        if mod(hours, 5) == 0
            %% Generate and save image of new relaxed jelly every 5 hours
            % Graft reorganization
            figure(1)
            figure_graft = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
            hold off
            xlim([0, 11 + offset]);
            ylim([0, 11 + offset]);
	    axis square; 

            % Need to explicitly create subfolders.
            path_graft = [path1, '/graft_reorganization'];
            cd(path_graft);
            i = hours;
            saveas(figure_graft, [num2str(i) '.jpg'])   
            pause(0.001)
            
            % % Contraction Stress
            % figure(2)
            % figure_contraction = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
            % hold on
            % quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.stress_muscle(:, 1), jelly.Nodes.stress_muscle(:, 2), 'AutoScale', 'on', 'LineWidth', 1.5);
            % 
            % hold off
            % xlim([0, 11 + offset]);
            % ylim([0, 11 + offset]);
            % axis square 

            % path_contraction = [path1, '/stress_contraction'];
            % cd(path_contraction);                                                            % write the image data
            % i = hours;
            % saveas(figure_contraction, [num2str(i) '.jpg'])    
            % pause(0.001)
            % 
            % % Pressure
            % figure(3)
            % figure_pressure = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
            % hold on
            % quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.pressure(:, 1), jelly.Nodes.pressure(:, 2), 'AutoScale', 'on', 'LineWidth', 1.5);
            % 
            % hold off
            % xlim([0, 11 + offset]);
            % ylim([0, 11 + offset]);
            % axis square; 
            % path_pressure = [path1, '/pressure'];
            % cd(path_pressure);                                                            % write the image data
            % i = hours;
            % saveas(figure_pressure, [num2str(i) '.jpg'])
            % pause(0.001)
            % 
            % % Elastic Stress
            % figure(4)
            % figure_elastic = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
            % hold on
            % quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.stress_elastic(:, 1), jelly.Nodes.stress_elastic(:, 2), 'AutoScale', 'on', 'LineWidth', 1.5);
            % 
            % hold off
            % xlim([0, 11 + offset]);
            % ylim([0, 11 + offset]);
            % axis square
            % path_elastic = [path1, '/stress_elastic'];
            % cd(path_elastic);                                                            % write the image data
            % i = hours;
            % saveas(figure_elastic, [num2str(i) '.jpg'])    
            % pause(0.001)
            % 
            % % Net Force
            % figure(5)
            % figure_net_force = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0, 'LineWidth', 1, 'NodeLabel', {});
            % hold on
            % quiver(jelly.Nodes.x_coord, jelly.Nodes.y_coord, jelly.Nodes.F_net(:, 1), jelly.Nodes.F_net(:, 2), 'AutoScale', 'on', 'LineWidth', 1.5);
            % 
            % hold off
            % xlim([0, 11 + offset]);
            % ylim([0, 11 + offset]);
            % axis square
            % path_stress_net = [path1, '/F_net'];
            % cd(path_stress_net);                                                            % write the image data
            % i = hours;
            % saveas(figure_net_force, [num2str(i) '.jpg'])    
            % pause(0.001)

            cd(path0); 
            
         
        end
        if mod(hours, 5) == 0
            aspect = aspect_ratio(jelly);
            a_r = cat(1, a_r, [aspect, hours]);
            [~, edge_idx] = area(jelly);
            vel_max = max((jelly.Nodes.velocity(edge_idx,1).^2) + jelly.Nodes.velocity(edge_idx,2).^2).^(1/2);
            vel = cat(1, vel, [vel_max, hours]);
        end  

    end

    %% final image
    figure_graft = plot(jelly, 'XData', jelly.Nodes.x_coord, 'YData', jelly.Nodes.y_coord, 'EdgeCData', jelly.Edges.strain0 + jelly.Edges.strain1, 'LineWidth', 1);
    hold off
    xlim([0, 11 + offset]);
    ylim([0, 11 + offset]);
    axis square 
    %% Plot aspect ratio
    cd(path1);                                                            % write the image data
    writematrix(a_r, 'a_r.xlsx');
    writematrix(vel, 'velocity.xlsx');
    writematrix([jelly.Nodes.x_coord, jelly.Nodes.y_coord], 'final_position.xlsx');
    cd(path0);  
