**L_system_rule_sequence_variable_length3.m**

modified from L-system_rule_sequence to take a vector of branch lengths
and store them in the btsl file

genereates a branching pattern using a set of rules for individual
generations, given as an array of 1 and 0 where 1 means parallel to
grandparent branch and 0 means orthogonal to grandparent branch.

note that the branches are color coded so that the terminal branches are
red and all the rest are blue.   this is very important and should not be
changed, since this is how the analysis program will be able to recognize
the terminal branches.....

note also that length array has EIGHT elements, while rule array has SIX.
 this is because we need to specify the length of B1 and B2 even though B1 does
 not have a rule to define its orientation since its orientation is part
 of the definition of the coordinate space, and B2 we are keeping fixed as orthogonal to the B0.

rule array input only has SIX elements, because we start specifying
orientations with B3.   B0 and B1 are used to determine the coordinate
frame, and for now we're going ot always assume that B2 is orthogonal to
B0, because otherwise dealing with the orthogonal vs paralle branches off
of the central B2 trifurcation becomes extremely difficult

example call:
%L_system_rule_sequence_variable_length('/Users/wallacemarshall/papers/Wei_branching_simulation/test_file_wholetree_010101.txt', [1 1 1 0.5 0.5 0.5 0.2], [0 1 0 1 0 1]);






generates a btsl representation using a self-similar relation as an L
system rule

length1, length 2 are the lengths of the trunk branch and the two
divergence main branch lengths that are 90 degrees to the trunk coming
right off of the end.

num_levels is how many levels of self-similar branching to do
length_rescale is the ratio of the new branch to the parent branch length

 at each level, generate two new branches oriented 180 degrees apart from each other and 90 degrees relative to the parent branch
so the only issue is the dihedral angle and relative position along the
parent branch
dihedral is fixed to be either 0 or 90.   0 means take SAME ORIENTATION
AS GRANDPARENT   90 means take an orientation mutually perpendicular to
parent and grandparent.   however to make the simulation more
interesting, we will make the choice of dihedral angles stochastic with a
probability of having an orthogonal angle give by dihedral_prob