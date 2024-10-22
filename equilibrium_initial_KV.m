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

function [jelly_eq, jelly_area, d_uncut] = equilibrium_initial_KV(A0, viscosity, damping_coefficient, bulk_modulus, dt, t, area_0, contraction_strength)

    %% Initialize an uncut jelly matrix
    [jelly, row_start, row_end] = offset_mesh(0); 
    [muscle_outer, muscle_inner] = whole_muscle();
    muscle_num = 2;
    
    edges = find_edges(muscle_outer);
    dist_rel = mesh_dist(jelly, row_start, row_end);
    j_area = area_initial(jelly, row_start, row_end);
    jelly_area = sum(j_area, 'all');
    area_relax = area_0*jelly_area;

    %% Arrays
    mesh_size = size(jelly);
    dist_current = zeros(mesh_size(1), mesh_size(2), 3); %mm
    strain = zeros(mesh_size(1), mesh_size(2), 3);
    elastic_stress = zeros(mesh_size(1), mesh_size(2), 3, 2);
    
    
    for time = 1:t
        %% Find the muscle stress at each node
        %define the current curve and length of the muscle
        stress_muscle = zeros(mesh_size(1), mesh_size(2), 2);
        for mus = 1:muscle_num
            c = muscle_outer(:,:,mus);
            d = muscle_inner(:,:,mus);
            for m = 1:length(muscle_outer)
                if m < length(muscle_outer)
                    node0 = c(:, m);
                    node1 = c(:, m+1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));

                    stress_muscle(node0(1), node0(2), 1) = stress_muscle(node0(1), node0(2), 1) + dx*contraction_strength;
                    stress_muscle(node0(1), node0(2), 2) = stress_muscle(node0(1), node0(2), 2) + dy*contraction_strength;
                end

                if m > 1
                    node0 = c(:, m);
                    node1 = c(:, m-1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));

                    stress_muscle(node0(1), node0(2), 1) = stress_muscle(node0(1), node0(2), 1) + dx*contraction_strength;
                    stress_muscle(node0(1), node0(2), 2) = stress_muscle(node0(1), node0(2), 2) + dy*contraction_strength;
                end
            end
            
            for m = 1:length(muscle_inner)
                if m < length(muscle_inner)
                    node0 = d(:, m);
                    node1 = d(:, m+1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));

                    stress_muscle(node0(1), node0(2), 1) = stress_muscle(node0(1), node0(2), 1) + dx*contraction_strength/2;
                    stress_muscle(node0(1), node0(2), 2) = stress_muscle(node0(1), node0(2), 2) + dy*contraction_strength/2;
                end

                if m > 1
                    node0 = d(:, m);
                    node1 = d(:, m-1);
                    x = jelly(node1(1), node1(2), 1) - jelly(node0(1), node0(2), 1);
                    y = jelly(node1(1), node1(2), 2) - jelly(node0(1), node0(2), 2);
                    dx = x/((x^2 + y^2)^(1/2));
                    dy = y/((x^2 + y^2)^(1/2));

                    stress_muscle(node0(1), node0(2), 1) = stress_muscle(node0(1), node0(2), 1) + dx*contraction_strength/2;
                    stress_muscle(node0(1), node0(2), 2) = stress_muscle(node0(1), node0(2), 2) + dy*contraction_strength/2;
                end
            end
        end
    
        %% Find the Elastic stress at each node
        for i = 1:mesh_size(1)
            for j = row_start(i):row_end(i)
                if j < row_end(i)
                    x1 = jelly(i, j+1, 1) - jelly(i, j, 1);
                    y1 = jelly(i, j+1, 2) - jelly(i, j, 2);
                    dist_current(i, j, 1) = (x1^2 + y1^2)^(1/2); %horiz node to right
                    strain(i,j,1) = (dist_current(i,j,1) - dist_rel(i,j,1))/dist_rel(i,j,1);
                    elastic_stress(i,j,1,1) = strain(i,j,1)*A0*x1/dist_current(i,j,1);
                    elastic_stress(i,j,1,2) = strain(i,j,1)*A0*y1/dist_current(i,j,1);
                end
                if i < mesh_size(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                    x2 = jelly(i+1, j+1, 1) - jelly(i, j, 1);
                    y2 = jelly(i+1, j+1, 2) - jelly(i, j, 2);
                    dist_current(i,j,2) = (x2^2 + y2^2)^(1/2); %node diag right
                    strain(i,j,2) = (dist_current(i,j,2) - dist_rel(i,j,2))/dist_rel(i,j,2);
                    elastic_stress(i,j,2,1) = strain(i,j,2)*A0*x2/dist_current(i,j,2);
                    elastic_stress(i,j,2,2) = strain(i,j,2)*A0*y2/dist_current(i,j,2);
                end
                if i < mesh_size(1) && j >= row_start(i+1) && j <= row_end(i+1)
                    x3 = jelly(i+1, j, 1) - jelly(i, j, 1);
                    y3 = jelly(i+1, j, 2) - jelly(i, j, 2);
                    dist_current(i,j,3) = (x3^2 + y3^2)^(1/2); %node diag left
                    strain(i,j,3) = (dist_current(i,j,3) - dist_rel(i,j,3))/dist_rel(i,j,3);
                    elastic_stress(i,j,3,1) = strain(i,j,3)*A0*x3/dist_current(i,j,3);
                    elastic_stress(i,j,3,2) = strain(i,j,3)*A0*y3/dist_current(i,j,3);
                end
            end
        end

        %Reset F_net
        stress_e = zeros(mesh_size(1), mesh_size(2), 2);
        
        %Each node is connected to 6 possible nodes
        for i = 1:mesh_size(1)
            %upper
            for j = row_start(i):row_end(i)
                if j < row_end(i)
                    stress_e(i,j,1) = stress_e(i,j,1) + elastic_stress(i,j,1,1); %horiz node to right
                    stress_e(i,j,2) = stress_e(i,j,2) + elastic_stress(i,j,1,2);
                end
                if i < mesh_size(1) && j+1 >= row_start(i+1) && j+1 <= row_end(i+1)
                    stress_e(i,j,1) = stress_e(i,j,1) + elastic_stress(i,j,2,1);
                    stress_e(i,j,2) = stress_e(i,j,2) + elastic_stress(i,j,2,2); %node diag right
                end
                if i < mesh_size(1) && j >= row_start(i+1) && j <= row_end(i+1)
                    stress_e(i,j,1) = stress_e(i,j,1) + elastic_stress(i,j,3,1);
                    stress_e(i,j,2) = stress_e(i,j,2) + elastic_stress(i,j,3,2); %node diag left
                end
            end
            %lower, stress is "minus" here because pos strain is away from
            %node that strain is calculated on
            for j = row_start(i):row_end(i)
                if j > row_start(i)
                    stress_e(i,j,1) = stress_e(i,j,1) - elastic_stress(i,j-1,1,1); %horiz node to left
                    stress_e(i,j,2) = stress_e(i,j,2) - elastic_stress(i,j-1,1,2);
                end
                if i > 1 && j-1 >= row_start(i-1) && j-1 <= row_end(i-1)
                    stress_e(i,j,1) = stress_e(i,j,1) - elastic_stress(i-1,j-1,2,1);
                    stress_e(i,j,2) = stress_e(i,j,2) - elastic_stress(i-1,j-1,2,2); %lower diag left
                end
                if i > 1 && j >= row_start(i-1) && j <= row_end(i-1)
                    stress_e(i,j,1) = stress_e(i,j,1) - elastic_stress(i-1,j,3,1);
                    stress_e(i,j,2) = stress_e(i,j,2) - elastic_stress(i-1,j,3,2); % lower diag right
                end
            end       
        end
        
        %% Calculate Pressure
        %current area
        j_area = sum(area_initial(jelly, row_start, row_end),'all');
        pressure = find_f_pressure_initial(jelly, area_relax, j_area, bulk_modulus, edges);


        stress_net = stress_e + pressure + stress_muscle/3; % Net stress in Pascals
        edge_area = 1e-3*1e-3;                              % In meters squared. 
        F_net = stress_net * edge_area;                     % Net force in Newtons.
        
        v = 1e3*F_net./damping_coefficient;                           % Factor of 1e3 converts meters per second to millimeters per second.
        displacement = v*dt;                                % In mm.
        jelly = jelly+displacement;                         % Update positions in jelly.
        

%          if mod(time,10)==0
%             %% image new relaxed jelly
%             for n = 1:11       %beads
%                 figure1 = scatter(jelly(n,row_start(n):row_end(n), 1), jelly(n,row_start(n):row_end(n), 2), 20, 'filled');
%                 hold on
%             end
%             hold off
%             xlim([-3, 13]);
%             ylim([-3, 13]);
%             pause(0.001)
%          end
    end
    d_uncut = dist_rel(6,6,1);
    jelly_eq = jelly;