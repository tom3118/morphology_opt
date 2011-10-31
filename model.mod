option randseed 0;
option solver snopt;
param N =4; # number of joints
param M = 1; # number of task points
param DIM = 3; # dimensionality of task space;
param LINKMIN = 0.1; # the length of the shortest link, should be > 0;
set tasks = {1..M}; # task labels, arbitrary
set activetasks; #subset of tasks used in objective
set joints = {1..N}; # joint indices, must be integer
set axes = {1..3}; #
set onetwo = {1..2};

# boxes defining each "task set"
param taskbounds{axes,tasks,1..2};

var t{axes,tasks}; # task points

param taskweights{tasks}; # determines how important that task is to the total
     #objective, e.g. 0 for tasks put in to represent dynamic feasibility 

# The main search variable "x" is the (x,y,z) joint locations for each task
var x{joints,axes,tasks};# := Uniform(0,1); 

# The variables of interest are the link lengths and the joint orientations
# these remain constant across tasks
var lengths{joints} >= LINKMIN; 

# the relative orientation of joints that share a link
# the cos and sin (resp) of the angle between adjacent joints
var rel{joints,1..2};# := Uniform(0,1);

# the orientation of the origin is still underspecified
# We define this by the unit vector orthogonal to the axis of the
# first link: rel[1,1],rel[1,2],relz0
var relz0;# := Uniform(0,1);

# Helper variables
var links{joints,axes,tasks};#:= Uniform(0,1); # link vectors

#joint unit vectors
#COMPATIBILITY in other code jointaxis[i,:,k,1] is z-hat_i
# and [...2] is x-hat
var jointaxis{joints,axes,tasks,onetwo} :=Uniform(0,1); 

# cross product of joint unit vectors, (helper var) gives y-hat
var orientations{joints,axes,tasks} := Uniform(0,1);

#scalars in said axis going from the origin to the task point
#var jointdisplacements{joints,tasks,onetwo} := Uniform(0,1); 

# the jacobian (actually this is it's transpose!)
var J{joints,axes,tasks} := Uniform(0,1); 
var JJT{axes,axes,tasks} default 0;
# important to initialize to something positiv definite
let{j in axes, k in tasks} JJT[j,j,k] := 1;  



#(max(JJT[1,1,k]*JJT[2,2,k]-JJT[1,2,k]*JJT[2,1,k],1e-6)))^(1/DIM)/(sum{i in joints} lengths[i])^2;

# JJT = JJ'
#subject to jjt {j1 in axes, j2 in axes, k in tasks}:
#     JJT[j1,j2,k] = sum{i in joints} J[i,j1,k]*J[i,j2,k];

# start at the origin
subject to origin {j in axes, k in tasks}:
     x[1,j,k] = 0;

# the link displacements are consistent with joint locations
subject to member {i in joints, k in tasks, j in axes: i > 1}:
     links[i-1,j,k] = x[i,j,k] - x[i-1,j,k];
     # and the end effector must reach the task point
subject to member_last {k in tasks, j in axes}:
     links[N,j,k] = t[j,k] - x[N,j,k];

# the task point must be in the task set
subject to task_set{k in tasks, j in axes}:
     taskbounds[j,k,1] <= t[j,k] <= taskbounds[j,k,2];


# the length of each link is constant
subject to length_c {i in joints, k in tasks}:
     lengths[i]^2 = sum{ j in axes} links[i,j,k]^2;

# ARTIFICIAL constraint to avoid possible degenerate cases
# in which link lengths are LINKMIN
# with more task points than DOF, this should perhaps be removed.
#subject to shortening {i in joints: i > 1}:
#     lengths[i] <= lengths[i-1];

# the joint axes are unit and orthogonal
subject to unit_vec {i in joints, k in tasks,l in onetwo}:
     sum{j in axes} jointaxis[i,j,k,l]^2 = 1;
subject to orthogonal {i in joints, k in tasks}:
     sum{j in axes} jointaxis[i,j,k,1]*jointaxis[i,j,k,2] = 0;

# the first joint axis is parallel to the incoming link
# the first axis associated with the origin is left free
# this, plus the previous constraint defines the second axis
subject to axis1 {i in joints, k in tasks, j in axes}:
     jointaxis[i,j,k,1]*lengths[i]  = links[i,j,k];

# the the outgoing link is coplanar with the axis
# together with "unit_vec", "orthogonal", and "axis1", this defines the
# axis up to a sign change.
# (c \dot \hat{a})\hat{a} + (c \dot \hat{b})\hat{b} = c
subject to coplanar {i in joints, j in axes, k in tasks:i>1}:
     (sum{j1 in axes} links[i-1,j1,k]*jointaxis[i,j1,k,1])*jointaxis[i,j,k,1] +
     (sum{j2 in axes} links[i-1,j2,k]*jointaxis[i,j2,k,2])*jointaxis[i,j,k,2]
     = links[i-1,j,k];

# define the joint orientation vector
subject to orientation1 {i in joints, k in tasks}:
     orientations[i,1,k] = jointaxis[i,2,k,1]*jointaxis[i,3,k,2] - 
     jointaxis[i,3,k,1]*jointaxis[i,2,k,2];
subject to orientation2 {i in joints, k in tasks}:
     orientations[i,2,k] = jointaxis[i,3,k,1]*jointaxis[i,1,k,2] -
     jointaxis[i,1,k,1]*jointaxis[i,3,k,2];
subject to orientation3 {i in joints, k in tasks}:
     orientations[i,3,k] = jointaxis[i,1,k,1]*jointaxis[i,2,k,2] -
     jointaxis[i,2,k,1]*jointaxis[i,1,k,2];

# the joint axes on opposite ends of a link have fixed orientation relative
# to eachother 
# (computed by (e_{i,1} \cross e_{i,2}) \cdot (e_{i-1,1} \cross e_{i-1,2})
# that is, 
# incidentally, this is just an angle, so we don't need a three vector to hold
# it.  But since we're keeping things quadratic, this seems okay

#subject to fixed_rel1 {i in joints, k in tasks: i > 1}:
#     rel[i,1] = orientations[i,2,k]*orientations[i,3,k] - 
#     orientations[i-1,3,k]*orientations[i-1,2,k];
#subject to fixed_rel2 {i in joints, k in tasks: i > 1}:
#     rel[i,2] = orientations[i,3,k]*orientations[i,1,k] - 
#     orientations[i-1,1,k]*orientations[i-1,3,k];
#subject to fixed_rel3 {i in joints, k in tasks: i > 1}:
#     rel[i,3] = orientations[i,1,k]*orientations[i,2,k] - 
#     orientations[i-1,2,k]*orientations[i-1,1,k];
#subject to fixed_origin {k in tasks,j in axes}:
#     rel[1,j] = orientations[1,j,k];

# TT: june 21, I was having trouble figuring that last one out so I'm rewriting
# it
# What we want is a constant angle, simply insisting on constant dot-product
# opens the door to sign flip so we use
# cos and sin, by dotting the orientation with the previous two (non-z) axes
subject to fixed_rel_cos {i in joints, k in tasks: i > 1}:
     rel[i,1] = sum{j in axes} orientations[i,j,k]*orientations[i-1,j,k];
subject to fixed_rel_sin {i in joints, k in tasks: i > 1}:
     rel[i,2] = sum{j in axes} jointaxis[i-1,j,k,1]*orientations[i,j,k];

# require the origin axis to be stationary
subject to fixed_origin {k in tasks, j in onetwo}:
     rel[1,j] = orientations[1,j,k];
subject to fixed_originz {k in tasks}:
     relz0 = orientations[1,3,k];

# Determine the projection of the task point in the joint basis
#subject to jointdisplacement {i in joints, k in tasks,l in onetwo}:
#     jointdisplacements[i,k,l] = sum{j in axes}
#              (t[j,k] - x[i,j,k])*jointaxis[i,j,k,l];

# The jacobian is simple in the joint basis
#subject to jacobian {i in joints, k in tasks, j in axes}:
#     J[i,j,k] = jointdisplacements[i,k,1]*jointaxis[i,j,k,2] - 
#     jointdisplacements[i,k,2]*jointaxis[i,j,k,1];

# The following will replace the preceding once we have tested a bit more
# we could also use a cross-product for this, 
# avoiding having the jointdisplacements variable
#     J[i,:,k] = cross(yhat_i,x[-1,:,0]-x[i,:,0])
subject to jacobian1 {i in joints, k in activetasks}:
     J[i,1,k] = orientations[i,2,k]*(t[3,k]-x[i,3,k])-
     orientations[i,3,k]*(t[2,k]-x[i,2,k]);
subject to jacobian2 {i in joints, k in activetasks}:
     J[i,2,k] = orientations[i,3,k]*(t[1,k]-x[i,1,k])-
     orientations[i,1,k]*(t[3,k]-x[i,3,k]);
subject to jacobian3 {i in joints, k in activetasks}:
     J[i,3,k] = orientations[i,1,k]*(t[2,k]-x[i,2,k])-
     orientations[i,2,k]*(t[1,k]-x[i,1,k]);


#data;

