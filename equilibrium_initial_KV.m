%{
======================================================================
    Function that simulates jellyfish for a specified set of parameters. 
======================================================================
    INPUT:
        elast0 (scalar):                Elasticity of spring
        viscosity (scalar):             Viscosity of dashpot in Maxwell Model
        bulk_modulus (scalar):          Bulk modulus of jellyfish 
        dt (scalar):                    Time step 
        N_time (int):                   Number of time steps
        area_0 (scalar):                ?
        contraction_strength (scalar):  ?
%}

function [jelly_eq, jelly_area, d_uncut] = equilibrium_initial_KV(elast0, viscosity, bulk_modulus, dt, N_time , area_0, contraction_strength)

    %% Initialize an uncut jelly matrix
    [jelly, row_start, row_end] = offset_mesh(0);       % Initialise jellyfish with no offset
    jelly_converted = convert_jelly_graph(jelly, row_start, row_end);

    [muscle_outer, muscle_inner] = whole_muscle();      % Initialise muscle layers
    muscle_num = 2;                                       % Number of muscle layers
    
    f1 = figure; 
    plot(jelly_converted, 'XData', jelly_converted.Nodes.x_coord, 'YData', jelly_converted.Nodes.y_coord);
    title(sprintf("Time: %g hours" ,0))
    hold on
    % xlim([0, 11 + offset]);
    % ylim([-1, 12]);
    
    edges = find_edges(muscle_outer);
    dist_rel = mesh_dist(jelly, row_start, row_end);
    j_area = area_initial(jelly, row_start, row_end);
    jelly_area = sum(j_area, 'all');
    area_relax = area_0*jelly_area;

    %% Arrays
    mesh_size = size(jelly);
    dist_current = zeros(mesh_size(1), mesh_size(2), 3); %mm
    strain = zeros(mesh_size(1), mesh_size(2), 3);
    F_elastic = zeros(mesh_size(1), mesh_size(2), 3, 2);
    
    for time = 1:N_time
        %% Find the muscle force at each node
        % define the current curve and length of the muscle
        F_muscle = zeros(mesh_size(1), mesh_size(2), 2);
        for mus = 1:muscle_num
            outer_muscle_nodes_index = muscle_outer(:,:,mus);
            inner_muscle_nodes_index = muscle_inner(:,:,mus);

            %% For outer muscles
            for m = 1:length(muscle_outer)
                if m < length(muscle_outer)
                    node0 = outer_muscle_nodes_index(:, m);
                    node1 = outer_muscle_nodes_index(:, m+1);
                    plot(jelly(node1(1), node1(2), 1), jelly(node1(1), node1(2), 2), 'or'); 
                    plot(jelly(node0(1), node0(2), 1), jelly(node0(1), node0(2), 2), 'or');
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));
                    
                    % Calculate muscle forces.
                    % TO DO: why does the increment in muscle forces here
                    % not scale with dt? As written now, there doesn't seem
                    % to be a reason for looping through time steps.
                    F_muscle(node0(1), node0(2), 1) = F_muscle(node0(1), node0(2), 1) + dx*contraction_strength;
                    F_muscle(node0(1), node0(2), 2) = F_muscle(node0(1), node0(2), 2) + dy*contraction_strength;
                end

                if m > 1
                    node0 = outer_muscle_nodes_index(:, m);
                    node1 = outer_muscle_nodes_index(:, m-1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));
                    
                    % Calculate muscle forces.
                    F_muscle(node0(1), node0(2), 1) = F_muscle(node0(1), node0(2), 1) + dx*contraction_strength;
                    F_muscle(node0(1), node0(2), 2) = F_muscle(node0(1), node0(2), 2) + dy*contraction_strength;
                end
            end
            
            %% For inner muscles
            for m = 1:length(muscle_inner)
                if m < length(muscle_inner)
                    node0 = inner_muscle_nodes_index(:, m);
                    node1 = inner_muscle_nodes_index(:, m+1);
                    plot(jelly(node1(1), node1(2), 1), jelly(node1(1), node1(2), 2), 'om'); 
                    plot(jelly(node0(1), node0(2), 1), jelly(node0(1), node0(2), 2), 'om');
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));
                    
                    % Calculate muscle forces.
                    F_muscle(node0(1), node0(2), 1) = F_muscle(node0(1), node0(2), 1) + dx*contraction_strength/2;
                    F_muscle(node0(1), node0(2), 2) = F_muscle(node0(1), node0(2), 2) + dy*contraction_strength/2;
                end

                if m > 1
                    node0 = inner_muscle_nodes_index(:, m);
                    node1 = inner_muscle_nodes_index(:, m-1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));
                    
                    % Calculate muscle forces.
                    F_muscle(node0(1), node0(2), 1) = F_muscle(node0(1), node0(2), 1) + dx*contraction_strength/2;
                    F_muscle(node0(1), node0(2), 2) = F_muscle(node0(1), node0(2), 2) + dy*contraction_strength/2;
                end
            end
        end
    
        %% Find the Elastic force at each node
        for i = 1:mesh_size(1)
            for j = row_start(i):row_end(i)
                if j < row_end(i)
                    x1 = jelly(i, j+1, 1) - jelly(i, j, 1);
                    y1 = jelly(i, j+1, 2) - jelly(i, j, 2);
                    dist_current(i, j, 1) = (x1^2 + y1^2)^(1/2); %horiz node to right
                    strain(i,j,1) = (dist_current(i,j,1) - dist_rel(i,j,1))/dist_rel(i,j,1);
                    F_elastic(i,j,1,1) = strain(i,j,1)*elast0*x1/dist_current(i,j,1);
                    F_elastic(i,j,1,2) = strain(i,j,1)*elast0*y1/dist_current(i,j,1);
                end
                if i < mesh_size(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                    x2 = jelly(i+1, j+1, 1) - jelly(i, j, 1);
                    y2 = jelly(i+1, j+1, 2) - jelly(i, j, 2);
                    dist_current(i,j,2) = (x2^2 + y2^2)^(1/2); %node diag right
                    strain(i,j,2) = (dist_current(i,j,2) - dist_rel(i,j,2))/dist_rel(i,j,2);
                    F_elastic(i,j,2,1) = strain(i,j,2)*elast0*x2/dist_current(i,j,2);
                    F_elastic(i,j,2,2) = strain(i,j,2)*elast0*y2/dist_current(i,j,2);
                end
                if i < mesh_size(1) && j >= row_start(i+1) && j <= row_end(i+1)
                    x3 = jelly(i+1, j, 1) - jelly(i, j, 1);
                    y3 = jelly(i+1, j, 2) - jelly(i, j, 2);
                    dist_current(i,j,3) = (x3^2 + y3^2)^(1/2); %node diag left
                    strain(i,j,3) = (dist_current(i,j,3) - dist_rel(i,j,3))/dist_rel(i,j,3);
                    F_elastic(i,j,3,1) = strain(i,j,3)*elast0*x3/dist_current(i,j,3);
                    F_elastic(i,j,3,2) = strain(i,j,3)*elast0*y3/dist_current(i,j,3);
                end
            end
        end

        %Reset F_net
        F_e = zeros(mesh_size(1), mesh_size(2), 2);
        
        %Each node is connected to 6 possible nodes
        for i = 1:mesh_size(1)
            %upper
            for j = row_start(i):row_end(i)
                if j < row_end(i)
                    F_e(i,j,1) = F_e(i,j,1) + F_elastic(i,j,1,1); %horiz node to right
                    F_e(i,j,2) = F_e(i,j,2) + F_elastic(i,j,1,2);
                end
                if i < mesh_size(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                    F_e(i,j,1) = F_e(i,j,1) + F_elastic(i,j,2,1);
                    F_e(i,j,2) = F_e(i,j,2) + F_elastic(i,j,2,2); %node diag right
                end
                if i < mesh_size(1) && j >= row_start(i+1) && j <= row_end(i+1)
                    F_e(i,j,1) = F_e(i,j,1) + F_elastic(i,j,3,1);
                    F_e(i,j,2) = F_e(i,j,2) + F_elastic(i,j,3,2); %node diag left
                end
            end
            %lower, stress is "minus" here because pos strain is away from
            %node that strain is calculated on
            for j = row_start(i):row_end(i)
                if j > row_start(i)
                    F_e(i,j,1) = F_e(i,j,1) - F_elastic(i,j-1,1,1); %horiz node to left
                    F_e(i,j,2) = F_e(i,j,2) - F_elastic(i,j-1,1,2);
                end
                if i > 1 && j-1 >= row_start(i-1) && j-1 <= row_end(i-1)
                    F_e(i,j,1) = F_e(i,j,1) - F_elastic(i-1,j-1,2,1);
                    F_e(i,j,2) = F_e(i,j,2) - F_elastic(i-1,j-1,2,2); %lower diag left
                end
                if i > 1 && j >= row_start(i-1) && j <= row_end(i-1)
                    F_e(i,j,1) = F_e(i,j,1) - F_elastic(i-1,j,3,1);
                    F_e(i,j,2) = F_e(i,j,2) - F_elastic(i-1,j,3,2); % lower diag right
                end
            end       
        end
        
        %% Pressure force
        %current area
        j_area = sum(area_initial(jelly, row_start, row_end),'all');
        F_pressure = find_f_pressure_initial(jelly, area_relax, j_area, bulk_modulus, edges);


        F_net = F_e + F_pressure + F_muscle/3; % Calculate net force.
        v = F_net./viscosity;                  % Calculate velocity due to viscous response to net force.
        displacement = v*dt;                   % Calculate displacement from velocity.
        jelly = jelly+displacement;            % Add displacements to jellyfish.
        
         % 
         % if mod(time,10)==0
         %    %% image new relaxed jelly
         %    for n = 1:11       %beads
         %        figure1 = scatter(jelly(n,row_start(n):row_end(n), 1), jelly(n,row_start(n):row_end(n), 2), 20, 'filled');
         %        hold on
         %    end
         %    hold off
         %    xlim([-3, 13]);
         %    ylim([-3, 13]);
         %    pause(0.001)
         % end
    end
    d_uncut = dist_rel(6,6,1);
    jelly_eq = jelly;

    close(f1); % close figure