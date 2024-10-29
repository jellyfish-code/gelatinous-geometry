%{
======================================================================
    Function that initialises a set of nodes for the muscle layers for
    butterfly graft.
======================================================================

INPUT:
        None.

OUTPUT:
        muscle_top_outer (array of doubles):      Matrix of indices identifying outer edge of top muscle layer.
        muscle_top_inner (array of doubles):      Matrix of indices identifying inner edge of top muscle layer.
        muscle_lowleft_outer (array of doubles):  Matrix of indices identifying outer edge of left side of bottom muscle layer. 
        muscle_lowleft_inner (array of doubles):  Matrix of indices identifying inner edge of left side of bottom muscle layer.
        muscle_lowright_outer (array of doubles): Matrix of indices identifying outer edge of right side of bottom muscle layer.
        muscle_lowright_inner (array of doubles): Matrix of indices identifying inner edge of right side of bottom muscle layer.
%}

function [muscle_top_outer, muscle_top_inner, muscle_lowleft_outer,  muscle_lowleft_inner, muscle_lowright_outer, muscle_lowright_inner] = butterfly_muscle()
    muscle_toprow_outer1 = [1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 7];
    muscle_toprow_outer2 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
    muscle_toprow_inner1 = [2, 3, 4, 5, 6, 7, 7, 7, 8, 9, 8, 8, 8, 8, 8];
    muscle_toprow_inner2 = [1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 9, 10, 11, 12, 13];
    
    muscle_botleft_outer1 = [11, 11, 11, 11, 11, 11];
    muscle_botleft_outer2 = [1, 2, 3, 4, 5, 6];
    muscle_botleft_inner1 = [10, 10, 10, 10, 10, 9];
    muscle_botleft_inner2 = [1, 2, 3, 4, 5, 5]; 
    
    muscle_botright_outer1 = [12, 13, 14, 15, 16, 17];
    muscle_botright_outer2 = [8, 9 , 10, 11, 12, 13];
    muscle_botright_inner1 = [11, 12, 13, 14, 15, 16];
    muscle_botright_inner2 = [9, 9, 10, 11, 12, 13];
    
    muscle_top_outer = cat(1, muscle_toprow_outer1, muscle_toprow_outer2);
    muscle_top_inner = cat(1, muscle_toprow_inner1, muscle_toprow_inner2);
    muscle_lowleft_outer = cat(1, flip(muscle_botleft_outer1), flip(muscle_botleft_outer2));
    muscle_lowleft_inner = cat(1, flip(muscle_botleft_inner1), flip(muscle_botleft_inner2));
    muscle_lowright_outer = cat(1, flip(muscle_botright_outer1), flip(muscle_botright_outer2));
    muscle_lowright_inner = cat(1, flip(muscle_botright_inner1), flip(muscle_botright_inner2));