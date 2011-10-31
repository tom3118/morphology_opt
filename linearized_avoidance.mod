###
# For each obstacle we choose a single constraint to maintain
###
param Nobs;
set obstacles = {1..Nobs};
# triangle mesh describing obstacles
param obsts{axes,1..3,obstacles};

param lintol = 1e-2;

# linearized constraints A, b such that A [d,u,v].T < b
param A{joints,1..3,obstacles,tasks};
param b{joints,obstacles,tasks};

# the u,v,d intersection point of each triangle and each link
var int_v{obstacles,joints,tasks};
var int_u{obstacles,joints,tasks};
var int_d{obstacles,joints,tasks};

# compute where these are
subject to intersection_pts{o in obstacles,a in axes,j in joints, k in tasks:
     j > 1}:
     -lintol <= x[j-1,a,k]+int_d[o,j,k]*(x[j,a,k]-x[j-1,a,k]) - obsts[a,1,o] +
     int_u[o,j,k]*(obsts[a,2,o] - obsts[a,1,o]) +
     int_v[o,j,k]*(obsts[a,3,o] - obsts[a,1,o]) <= lintol;
subject to intersection_pts1{o in obstacles,a in axes,k in tasks}:
     -lintol <= int_d[o,j,k]*(x[1,a,k]) - obsts[a,1,o] +
     int_u[o,j,k]*(obsts[a,2,o] - obsts[a,1,o]) +
     int_v[o,j,k]*(obsts[a,3,o] - obsts[a,1,o]) <= lintol;


# maintain the non-intersection constraint
subject to non_intersection{o in obstacles, j in joints, k in tasks}:
     A[j,1,o,k]*int_d[o,j,k] +      
     A[j,2,o,k]*int_u[o,j,k] + 
     A[j,3,o,k]*int_v[o,j,k] <= b[j,o,k];
