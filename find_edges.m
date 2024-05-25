function edges = find_edges(muscle_outer)

	% Edges must go clockwise around the jelly
    edges = cat(2, muscle_outer(:,:,1), muscle_outer(:,:,2), muscle_outer(:,1,1));
    insert = zeros(2,1);
    i = 1;
    a = edges(2,i);
    while i <= length(edges)-1
        if edges(1,i+1) == edges(1,i) && edges(2,i+1) == edges(2,i)
            edges(:,i) = [];
        elseif abs(edges(2, i+1) - a) <= 1
            i = i+1;
            a = edges(2,i);
        elseif edges(2, i+1) - a > 1
            a = a+1;
            insert(1,1) = edges(1,i);
            insert(2,1) = a;
            edges = cat(2, edges(:,1:i), insert, edges(:,i+1:length(edges)));
        elseif edges(2, i+1) - a < 1
            a = a-1;
            insert(1,1) = edges(1,i);
            insert(2,1) = a;
            edges = cat(2, edges(:,1:i), insert, edges(:,i+1:length(edges)));
        end
    end