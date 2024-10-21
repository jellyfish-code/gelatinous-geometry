%{
======================================================================
    Function that (?)
======================================================================

INPUT:
        jelly (?):
        row_start (?): 
        row_end (?): 

OUTPUT:
        jelly_area (?):
%}

function jelly_area = area_initial(jelly, row_start, row_end)
    a = size(jelly);
    jelly_area = zeros(a(1), a(2), 2);
    for i = 1:length(row_start)
        for j = row_start(i):row_end(i)
            %a1, upper right diagonal finite element
            if i<length(row_start) && j<row_end(i) && j+1>=row_start(i+1) && j+1<=row_end(i+1)
                %then a1 exists, find area with shoelace method
                %coordinates are (i,j), (i+1, j+1), (i,j+1)
                a1 = abs(0.5*(jelly(i,j,1)*jelly(i+1,j+1,2) + jelly(i+1,j+1,1)*jelly(i,j+1,2) + jelly(i,j+1,1)*jelly(i,j,2) ...
                    - jelly(i,j,2)*jelly(i+1,j+1,1) - jelly(i+1,j+1,2)*jelly(i,j+1,1) - jelly(i,j+1,2)*jelly(i,j,1)));
                jelly_area(i,j,1) = a1;
            end
            %a2, lower right diagonal finite element
            if i>1 && j<row_end(i) && j<=row_end(i-1) && j>=row_start(i-1)
                %then a2 exists
                %coordinates are (i,j), (i-1,j), (i,j+1)
                a2 = abs(0.5*(jelly(i,j,1)*jelly(i-1,j,2) + jelly(i-1,j,1)*jelly(i,j+1,2) + jelly(i,j+1,1)*jelly(i,j,2) ...
                    - jelly(i,j,2)*jelly(i-1,j,1) - jelly(i-1,j,2)*jelly(i,j+1,1) - jelly(i,j+1,2)*jelly(i,j,1)));
                jelly_area(i,j,2) = a2;
            end
        end
    end
    
            
        