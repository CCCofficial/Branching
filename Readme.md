
# Author: Wallace Marshall, UCSF #

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



**btsl_draw_tree_different_angles_13.m**


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


modified from version 12 to clean up unused bits of code that had been
commented out as well as to modify the steric clash detection to take
into account Wei's measurements of branch diameter

modified from 11c_fixRotation to add back the other half of the kidney

modified from 11c to NOT accumulate rotations from parent branches, since
Wei's angle changes are based on a global orientation inherited directly
from an earlier branch and not propagated throuh intervening ancestors


testing-only version of 11c which only draws half the tree

 modified from version 11b to check if an angle change direction is
 impossible (meaning, the angles are being asked to collapse in a
 direction that is in-plane with the two branches but perpendicular to
 the parent.   in such a case, do not collapse the angles at all in that
 direction.   it is unclear how to handle cases like this so this
 solution is chosen because it is simplest.

modified from verison 11 to generate the trunk branch and add it to the
overall branch structure.

NOTE:   the collapse orientation vector entries tell how many branching
levels before a given branch the orientation is to be found, that is then
applied to the DAUGHTERS of the branch in question!!!   so ifyou want to
change the angles for a pair of branches relative to their grandparent,
that would be a value of 1 since 0 is the parent!!

modifiecation from verison 10 - now the collapse_angle_directions array
can contain non-negative integers, each of which specified how many
generations back from the parent the reference branch is that determines
the direction of angle collapse.  so, 0 means the direction of the
parent, 1 means direction of the grandparent, 2 means great-grandparent,
and so on.   this will let us explore a wider range of possibilities,
plus in the case that the actual kidney always refers to teh original
trunk branks for all angle changes, this will be able to handle that case
as well as visually indicate how many generationns of orientation data
the system must be keeping track of, if it is only using local
information, which givesn a measure of how much less parsimonious a local
only scheme would become relative to a global positional information
scheme.

note about branch levels.   level 1, called the main branch here, is
actually what Wei called the "secondary" branch, ie. the first set of
branches off of the main trunk.   level 2 is the tertiary branch, which
is the one that is a trifurcation.   so, the whole btsl input file to
this program assumes the main branch is ALREADY PRESENT but implicit.
so there is no main branch in the representation here.   the convention
we used in making the btsl file is that the secondary branch is parallel
to the Z axis (which is why we reflect all coordinates across the XY
plane in this program when we want to make the second half of the
kidney), and then the trifurcations happen parallel to X and -X, plus Z.
 so, the main trunk would be parallel to the Y axis!!!   so, in the case
 that a rotation needs to refer back to the implicit trunk, we need to
either handle that as an exceptional case or else add one more branch
which is the invisible trunk.  it can have as its base 0 -1 0 and as its
end 0 0 0 that would make everythign consistent!

modification from version 9 *** MAJOR MODIFICAITON ****
modified to allow the angle change to take place in either of three
possible direction, which we define by the numbers 1 or 2.
 0 means the angle change takes place so as to make the daughter branches
   more paralell with the parent branch, which is what was done in version
   9 and all previous versions with angle changes.
 1 means the angle change takes place so as to make the daughter branches
   more paralell to the grandparent branch

so, in addition to the angle array, the program will take another input
called collapse_direction_directions which contains values 1 or 2

note:   the first entry has to be 1 because there is no grandparent to
branch parallel to!!!!

note also that certain combinations of branching patterns and angle
collapse orientations are not defined, i..e whenever a branch is paralell to its
grandparent, it is not possible to collapse angles in that direction.
in those cases, the method will spit out zero for all results




modified from verison 8 to also output the skewness of the radial
distance distribution as well as the entropy, since we want to get an
idea of how much all the distances are bunched up at a particular
distance from the centroid.

modified from version 7 to calculate steric clash between all branches

**** NOT WORKING YET ******

usage version 7
%btsl_draw_tree_different_angles_7('/Users/wallacemarshall/papers/Wei_branching_simulation/test_file_wholetree_L_000000.txt', [180 180 90 60 60 75 90 120], 1)


modified from version 6 to plot and analyze the distribution of tips in
angle space in spherical coordinates

modified from btsl_draw_tree_different_angles_5 to calculate distance to
nearest neighbor distribution for endpoints as a measure of clash and to
calculate the correlation dimension of the endpoint distribution as a
measure of the fractal dimension of the nephron distribution

usage:
btsl_draw_tree_different_angles_6('/Users/wallacemarshall/papers/Wei_branching_simulation/test_file_wholetree_L_000000.txt', [180 180 90 60 60 75 90 120], 1)

draw_plots is a flag 1 means draw all the figures, 0 means suppress the
figures




modified from v4 to do more spatial statistics

modified from v3 to calculate spatial statistics etc

modified from v 2 to apply rotation to the first branch

modified to read the lengths from the input file

modified from btsl_draw_tree_change_angles5 to allow the angles to all be
different, so instead of taking a single parameter "new angle" it takes a
vector of angles one for each level of branching.   these angles are the
angles BETWEEN the two daughter branches at that level.   also we need to
specify the lengths of each branch.   *** also note that from now on we
should calculate and plot the WHOLE tree, not just starting with B3.   so
our angle and length vecotrs need to have eight elements not six.
 L_system_rule_sequence starts with B1 as the parent branch, in order to
 render and analyze the whole kidney we need to make a copy of the entire
 final structure and reflect it across the Z axis.   so we only need to
 do the calculate for one half of the kidney and then we just generate
 the other half by reflection.

****** from this version onwards, ALL ANGLES ARE IN DEGREES!!!!

 note about LENGTHS:   all lengths are specified in the L-system rule
 sequence program so we don't worry about them here.   they are
 hard-encoded in the btsl input file.



%modified from btsl_draw_tree_change_angles3 to put in the cmyk color
coding scheme from Wei

wei's color code
  branch     c    M    Y     K        converts to:    r      g       b
3            60   90   0     0                       102      25      255
  4            70    15    0    0                 77       217    255
  5         50  0   100    0                     128      255        0
    6         0   35    85     0                 255     166        38
converted using online calculator:
http://easycalculation.com/colorconverter/rgb-color-converter.php


modified from btsl_analyze_endpoints_convexhull_2 to allow the angles
between branches to change over time as Wei has seen.



only draw the hull if drawhull = 1

modified from btsl_analyze_endpoints_convexhull to only start plotting
with the third branch and then to color each branch according to the
color code to from wei   

wei's color code
  branch     c    M    Y     K
4            60   90   0     0
  5            70    15    0    0
  6         50  0   100    0
    7         0   35    85     0

so we want to keep track of colors

 since we are starting with the 3rd branch, we will only plot L1.1.1 and
 its children


modified from btsl_analyze_endpoints  to define distance to the surface
based on generating a convex hull containing all the endpoints 

modified from btsl_parse_render.   program to take a btsl specification
and analyze the spatial arrangement of the endpoints of all terminal
branches.   for the first siple analysis, take the centroid of all the
terminal points and then measure how far each one is away from the
centroid and then make a histogram of the distances.   as a final value, compute the average distance
of all points from the centroid.

 program to read a branching tissue specification language (btsl) file
 and parse it into a tree representation that can be used to generate
 simulated branch images and analyze the overall organ structure

language format:

 individual branch description:
 {Branch_ID  -L length -D diameter -O orientation (format [a,b,c]) -C
 color -S shift in position (as ratio of full length of parent branch)
    sub_branch_1
    sub_branch_2 ....
  }

example

 {L1 -L   5 -D 10 -O [0, 0, 1] -C r -S 0.5 {  L1.1  -L 4 -D 9 -O [0, 1, 0] -C r -S 0.5}
 {  L1.2  -L 3 -D 8 -O [0, -1, 0] -C r -S 0}
 {  L1.3  -L 2 -D 7 -O [0, 0, 1] -C r -S 0
      {  L1.3.1  -L 2 -D 7 -O [0, 1, 0] -C b -S 0 }
%
 }
 }



 tissue:
 Branch_description_x {Branch_description_x_1  }  {Branch_description_x_2
 }   {Branch_description_x3   } ......

  these can be nested infinitely.   a branch can have any number of child
  branches each of which is described by an entry in curly braces


**enumerate_branching_patterns_v10.m**
modified to call btsl_draw_tree_different_angles_13 which calculates
steric clash using measured value for branch diameter instead of
estimating it as a ratio of branch length as in previous versions

%enumerate_branching_patterns_v10('/Users/wallacemarshall/papers/Wei_branching_simulation/btsl_temp_file.txt', 8, [105 50 95 157 178 173 96 61], [176 180 81 59 65 83 103 125], [0 0 0 1 2 3 4 5])



Conventions for lenth and angle array:

length array the first entry is NOT the trunk branch (a.k.a. primary branch which is yellow in Wei's drawings), rahter it is the
first pairof branches that come out of the main trunk, also called
"secondary branch".  it is OK that we don't give a length for the primary
branch because it has no influence on any of the spatial statistics for
endpointdistribution.

the SECOND entry is the trifurcating branch, also called "tertiary
branch" which is colored light blue in the Figure.

conventions for rule sequences:
 each entry in the rule vector specifies the branching orientation of the
 CHILDREN of the corresponding branch.  so the first entry concerns the
 tertiary branch and tells which way its children (i.e. the fourth branch) will be oriented.
 the final entry describes the second-to-last branching level and
 specified which way the terminal branches will be oriented.
 the orientation of the tertiary branch itselfis always fixed to
 being perpendicular to the primary branch.  this is because the
 initial simulations were focused on half kidneys so that the dihedral
 angle was not defined up through the tertiary branching level.  this is
 a constraint that coudl be relaxed in the future.

as an example, if the 4th and 5th branches are orthogonal to their
grandparents, but then all the ones after that are parallel, as is the
case in the kidney, then the correct rule string should be 001111






modified from version 8 to use btsl_draw_tree_different_angles_12 in
which we have corrected the fact that earlier versions (11c) were
accumulating rotations generated by angle collapse of earlier branches.
this needed to be fixed since Wei's new scheme for angle collapse is
strictly in reference to a single branch, and is therefore really a
global cue not a local cue inherited from an immediate ancester.  this
new version (12) has been verified by examining the tree diagrams to
generate correct direction of angle collapse.   hence any results
generated using enumerate branching pattersn v8 are invalid.

modified from version 7 to use btsl_draw_tree_different_angles_11c which
in the case of an invalid direction change (i.e. a direction change
vector that is coplanar with the two branches but perpendicular to the
parent) simply doesnt' change the angles for that pair of branches.
this is necessary beause for certain collapse arrange and branching
patterns there will be invalid direciton changes during enumeration

example call
enumerate_branching_patterns_v7('/Users/wallacemarshall/papers/Wei_branching_simulation/btsl_temp_file.txt', 8, [105 50 95 157 178 173 96 61], [176 180 81 59 65 83 103 125], [0 0 0 1 2 3 4 5])


example call to the underlying program just for reference
%btsl_draw_tree_different_angles_11c('/Users/wallacemarshall/papers/Wei_branching_simulation/test_file_wholetree_001111b.txt', [180 180 90 60 60 75 90 120], [0 0 0 1 2 3 4 5], 1)


modified from version 6 o use btsl_draw_tree_different_angles_11b which
includes the trunk branch in the tree structure

modified from verison 5b to use btsl_draw_tree_different_angles_11,
whichi allows the collapse direction array to refer to any earlier
generation of branching, i.e. 0 is parent, 1 is grandparent, 2 is
great-grand-parent, and so on, as a way to represent the actual situation
in the kidney where it may be the case that angle changes are carried out
in directions parallel to the initial trunk, such that each new
generation of branches has to refer to an earlier and earlier generation
branch.    so now collapse array can have any value up to 8.


modified from version 5 in that it does not enumerate over all collapse
angle combinations, but just takes one collapse angle array as an input.
in this array, 0 means collapse parallel to the parent, and 1 is collapse
parallel to the grandparent branch.



 MAJOR MODIFICATION - calles btsl_draw_tree_different_angles_v10 which
 allows angles to change in either of two directions.   this means that
 now in addition to enumerating over branching patterns, for each
 branching pattern we also enumerate over all collapse directions.  for
 this enumeration, we define 1 means the angles change to make the
 daughters aligned to the parent, and 2 means the angles change to make
 the daughters aligned to the grandparent.   

modified rom verison 3 to use btsl_draw)tree_different_angles_v98 which
also calculates skewness and entropy of radial distance distribution

modified from version 2 to use btsl_draw_tree_different_angles_v8 which
calculates and reports steric clash.

modified from L-system_enumerate to use variable lengths and variable
angles.   

example  enumerate_branching_patterns_v2('/Users/wallacemarshall/papers/Wei_branching_simulation/btsl_temp_file.txt', 8, [100 50 200 200 200 200 100 60], [180 180 90 60 60 75 90 120])


length array based on wei and smythe data
[100 50 200 200 200 200 100 60]
note that this has eight values whereas the rule arrays that are
enumerated only have six because the first two branches never change and
so they are not enumerated over.

angle array based on Wei and Smythe data
%[180 180 90 60 60 75 90 120]
note this also has eight entries for the eight different angles
the first angle is the first branch off the main trunk, and the second is
the trifurcation.


modified from L-system_enumerate to use a convex hull to measure distance
of each end-point to the surface, where in this case the surface will be
defined as the convex hull containing all the endpoints

enumerate all outcomes for all possible rule lists of a given length
have to give the lengths of the first two branches which don't count in
the rule list, and then the length rescale factor, and the number of
levels (total - the first two count in this total)
