function result_list = L_system_rule_sequence_variable_length3( outfile, length_array, rule_array )


% modified from L-system_rule_sequence to take a vector of branch lengths
% and store them in the btsl file

% genereates a branching pattern using a set of rules for individual
% generations, given as an array of 1 and 0 where 1 means parallel to
% grandparent branch and 0 means orthogonal to grandparent branch.

% note that the branches are color coded so that the terminal branches are
% red and all the rest are blue.   this is very important and should not be
% changed, since this is how the analysis program will be able to recognize
% the terminal branches.....

% note also that length array has EIGHT elements, while rule array has SIX.
%  this is because we need to specify the length of B1 and B2 even though B1 does
%  not have a rule to define its orientation since its orientation is part
%  of the definition of the coordinate space, and B2 we are keeping fixed as orthogonal to the B0.

% rule array input only has SIX elements, because we start specifying
% orientations with B3.   B0 and B1 are used to determine the coordinate
% frame, and for now we're going ot always assume that B2 is orthogonal to
% B0, because otherwise dealing with the orthogonal vs paralle branches off
% of the central B2 trifurcation becomes extremely difficult

% example call:
%L_system_rule_sequence_variable_length('/Users/wallacemarshall/papers/Wei_branching_simulation/test_file_wholetree_010101.txt', [1 1 1 0.5 0.5 0.5 0.2], [0 1 0 1 0 1]);






%generates a btsl representation using a self-similar relation as an L
%system rule

% length1, length 2 are the lengths of the trunk branch and the two
% divergence main branch lengths that are 90 degrees to the trunk coming
% right off of the end.

% num_levels is how many levels of self-similar branching to do
% length_rescale is the ratio of the new branch to the parent branch length

%  at each level, generate two new branches oriented 180 degrees apart from each other and 90 degrees relative to the parent branch
% so the only issue is the dihedral angle and relative position along the
% parent branch
% dihedral is fixed to be either 0 or 90.   0 means take SAME ORIENTATION
% AS GRANDPARENT   90 means take an orientation mutually perpendicular to
% parent and grandparent.   however to make the simulation more
% interesting, we will make the choice of dihedral angles stochastic with a
% probability of having an orthogonal angle give by dihedral_prob




 branch_levels = length(rule_array);





fid1 = fopen(outfile,'w');

% write out the trunk information (i.e. B1)
fprintf(fid1, '%s %i %s \n', '{ L1 -L ', length_array(1), '  -D 10 -O [0, 0, 1] -C b -S 1.0 ');



% write out left main branch (the left, right and upper main branches are
% B2 in Wei's numbering scheme)

fprintf(fid1, '%s %i %s \n', '  {  L1.1  -L ', length_array(2), ' -D 9 -O [1, 0, 0] -C b -S 1.0 ');
generate_branch(fid1, branch_levels, length_array, 0, 0, 1, 1, 0, 0, rule_array, 'L1.1', branch_levels);





% finish left main branch
fprintf(fid1, '%s \n', '}');




% write out right main branch

fprintf(fid1, '%s %i %s \n', '  {  L1.2  -L ', length_array(2), ' -D 8 -O [-1, 0, 0] -C b -S 1.0 ');
generate_branch(fid1, branch_levels, length_array, 0, 0, 1, -1, 0, 0, rule_array, 'L1.2', branch_levels);





% finish right main branch
fprintf(fid1, '%s \n', '}');




% write out upper main branch

fprintf(fid1, '%s %i %s \n', '  {  L1.3  -L ', length_array(2), ' -D 8 -O [0, 0, 1] -C b -S 1.0 ');


generate_branch(fid1, branch_levels, length_array, 0, 0, 1, 0, 0, 1, rule_array, 'L1.3', branch_levels);





% finish upper main branch
fprintf(fid1, '%s \n', '}');




% finish the trunk
fprintf(fid1, '%s \n', '}');





end





function t = generate_branch(fid1, remaining_branches, length_array, grandparent_a, grandparent_b, grandparent_c, parent_a, parent_b, parent_c, list_of_rules, name_root, branch_levels)
% remaining_branches is how many more levels to call in the recursive call
% name_root is the current branch name
% parent a b c are the three direction cosines in the orientation vector of
% the parent branch
% grandparent a b c are the direction cosines of the grandparent branch
% which is needed to get the dihedral angle

calling_generate_branch_on = name_root;

branch_depth = branch_levels - remaining_branches; % so this should be zero the first time it is called

%new_length = length_array(branch_depth + 3);
new_remaining_branches = remaining_branches - 1;


if remaining_branches > 0  % put this before calculate new length to avoid accessing out of bounds on final recursive call

new_length = length_array(branch_depth + 3);



% test for the central branch of the trifurcation which has to be handled
% specially
if (grandparent_a == parent_a) && (grandparent_b == parent_b) && (grandparent_c == parent_c)
    special_case = 1;
else
    special_case = 0;
end


%if remaining_branches > 0
    
    dihedral = list_of_rules(1);
    rule_size = length(list_of_rules);
    new_rule_list = list_of_rules(2:rule_size);

    
    
    
        
        
    % call on the left sub-branch
    new_branchname = strcat(name_root, '.1');

    if dihedral==1
       new_a = grandparent_a;
       new_b = grandparent_b;
       new_c = grandparent_c;
        
    else % compute mutually orthogonal orientation vector
       new_a = not(logical(abs(parent_a)) | logical(abs(grandparent_a)));
       new_b = not(logical(abs(parent_b)) | logical(abs(grandparent_b)));
       new_c = not(logical(abs(parent_c)) | logical(abs(grandparent_c)));
        
        
    end
    
    
    % that procedure works except for the central branch of the
    % trifurcation, for which the parent and grandparent are parallel
    % so for that case we need to handle it specially.   in that case we
    % want the two branches to be parallel to B0 which is oriented 0 1 0
    if (special_case == 1) 
        if dihedral==0
        % if orthogonal we want it to have the same orientation as the
        % daughter branches from the two other B2 branches
            new_a = 0;
            new_b = 1;
            new_c = 0;
        else 
            % if parallel we want it to have the orientations of the B2
            % branches themselves
            new_a = 1;
            new_b = 0;
            new_c = 0;
        end
    end
    
    
    if remaining_branches == 1
        fprintf(fid1, '%s %s %s %i %s %i %s %i %s %i %s \n', '  {  ', new_branchname,  ' -L ', new_length, ' -D 8 -O [ ', new_a, ', ', new_b, ', ', new_c  ,'] -C r -S 1.0 ');
    else
        fprintf(fid1, '%s %s %s %i %s %i %s %i %s %i %s \n', '  {  ', new_branchname,  ' -L ', new_length, ' -D 8 -O [ ', new_a, ', ', new_b, ', ', new_c  ,'] -C b -S 1.0 ');
    end

    generate_branch(fid1, new_remaining_branches, length_array, parent_a, parent_b, parent_c, new_a, new_b, new_c, new_rule_list, new_branchname, branch_levels);

    fprintf(fid1, '%s \n', '}');




    % call on the right sub_branch
    new_branchname = strcat(name_root, '.2');

    new_a = -1*new_a; % the other branch is a reflection of the first one
    new_b = -1*new_b;
    new_c = -1*new_c;
    
    if remaining_branches == 1
        fprintf(fid1, '%s %s %s %i %s %i %s %i %s %i %s \n', '  {  ', new_branchname,  ' -L ', new_length, ' -D 8 -O [ ', new_a, ', ', new_b, ', ', new_c  ,'] -C r -S 1.0 ');
    else
        fprintf(fid1, '%s %s %s %i %s %i %s %i %s %i %s \n', '  {  ', new_branchname,  ' -L ', new_length, ' -D 8 -O [ ', new_a, ', ', new_b, ', ', new_c  ,'] -C b -S 1.0 ');
    end
        
        %fprintf(fid1, '%s %s %s %i %s \n', '  {  ', new_branchname,  ' -L ', new_length, ' -D 8 -O [-1, 0, 0] -C r -S 1.0 ');

    generate_branch(fid1, new_remaining_branches, length_array, parent_a, parent_b, parent_c, new_a, new_b, new_c, new_rule_list, new_branchname, branch_levels);
    %generate_branch(fid1, new_remaining_branches, new_length, length_rescale, dihedral_prob, new_branchname);
   
     fprintf(fid1, '%s \n', '}');
     
     
     
   
    
    
    
    
    
end  % exit and do nothing if down to zero remaining branches to generate

end

