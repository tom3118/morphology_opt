"""
This file provides functionality to create an ampl morphology optimization 
problem from a collection of feasible executions (for instance found by the
RRT* module).  This approach is preferable to the more "brute force" approach
in which we hook directly to SNOPT because AMPL provides extensive aids such as
symbolic differentiation.

Example useage: 
1) Given an initial morphology and a set of tasks, run RRT* to generate
trajectories, T, that accomplish them.
b) instantiate that morphology as a Morphology Object, M
2) construct an equivalent task spec, S, and obstacle environment, O, here 
(soon to be automated)
3) M.amplPrint(T,S,O) will create an ampl readible optimization problem such
that the return morphology
i) maximizes the weighted geometric manipulatability (
ii) can approach the task points with trajectories that stay on the "same side"
of obstacles as T.
4) Use RRT* to verify feasibility and compute new trajectories
5) iterate until M stops changing.

Note that in a strict sense, 3ii) does not gaurantee feasibility, but as a
practical matter it is extremely likely to.  There is a trade-off between the
number of constraints and the flexibility of the optimization problem: 
with many, the result is increasingly likely to be feasible, but less likely to
be far from (i.e., much better than) the initial morphology.  In practice it
seems best to keep these constraints as few as possible and re-check for 
feasibility using RRT*.  One can then iterate on the new trajectories which
will likely lead to different constraints.

In practice this is both effective and efficient for systems with about 6 DOF
 even with tight feasibility constraints (specifically, reaching into a small 
hole).  With 8-9 DOF it remains relatively efficient, but
the optimization is much more conservative, requiring many more iterations, and
depending more on initial conditions.
"""

from numpy import *
import pdb
import pyampl 
import models

from constraints import Trot, findLoosest, findSeparator

# set FIX to True to use AMPL to validate the initial conditions, 
# but not optimize
FIX=False
# Should the generated equality constraints (e.g. barycentric coords) be
# equality or within a tolerance?
EQUALITY=False
# Use the (newer) pure linearization model (nonintersection_linear.mod)
# or (older) nonintersection.mod
# The former is incompatible with findLoosest
# The latter is better suited to the findLoosest constraint generation
PURELIN=True


class Morphology(object):
  """
  This class encapsulates a manipulator morphology
  If we restrict our attention to rotary joints
  this consists of an initial joint axis, lengths and azimuths for each link

  The main use of this class is to create AMPL models, data, etc
  
  """
  def __init__(self, T0=eye(3), Nlinks=0, azimuths=[], lengths=[]):
    self.T0 = T0
    if azimuths <> []:
      self.azimuths=azimuths
      self.Nlinks=len(azimuths)
      if len(lengths)==len(azimuths):
        self.lengths = lengths
      else:
        self.lengths = ones(len(azimuths))
    elif lengths <> []:
      self.lengths = lengths
      self.Nlinks = len(lengths)
      self.azimuths = pi/2*ones(len(lengths))
    elif Nlinks:
      self.Nlinks = Nlinks
      self.lengths = ones(Nlinks)
      self.azimuths = pi/2*ones(Nlinks)
    else:
      print "expected at least one of {Nlinks,azimuths,lengths}"
      pdb.set_trace()
    self.ntaskpoints = 0
    self.azimuths[0] = 0
    self.activetasks = []

  def amplPrint(self,trajlist,taskboundarylist=None, obstacles=None):
    """
    @return an ampl string containing the complete model and data
    @param trajlist, a list of trajectories length L
    @param taskboundarylist, a list of task boundaries of length L
    @param obstacles ObstacleSet object
    """
    prefix = "param N = %s;\n"%self.Nlinks
    trajs = ''.join([self.toAmplTrajectory(trajlist[i],obstacles,
                                           taskboundarylist[i])
                     for i in range(len(trajlist))])
    prefix += "param M = %s;\n"%self.ntaskpoints
    prefix += pyampl.AmplDriver.writeSet('activetasks',self.activetasks)
    obs = obstacles and str(obstacles) or ''
    fix = FIX and models.fixall or ''
    cons = EQUALITY and models.constraints2 or models.constraints
    intersection_model = (PURELIN and models.nonintersection_model2 or 
                          models.nonintersection_model)
    return (prefix + models.model + cons + intersection_model + obs + 
            self.toAmplData() + models.obj_model + trajs + fix)

  def toAmplData(self):
    """
    @return a string representing this object in a data format consistent with
    a specific AMPL model
    """
    return self._toAmplStartPoint(pi/2*ones(self.Nlinks),onlyself=True)

  def toAmplState(self,state,obstacles=None,taskboundaries=None):
    """
     This function takes a set of joint angles and the values of all of the
    Ampl model variables, returning them as a string

    @param state: a list of joint angles
    @param obstacles: an ObstacleSet object
    @param taskboundaries: if present, this should be in [(min,max)] format,
    if omitted, we assume this task is artificial 
    (defaults to a wide boundaries)
    """
    return self._toAmplStartPoint(state,taskboundaries,obstacles)

  def toAmplTrajectory(self,trajectory,obstacles=None,taskboundaries=None):
    """
    Calls toAmplState on each state in the trajectory, putting a boundary
    on the LAST ONLY.

    @param trajectory: array N by self.Nlinks of joint states
    @param obstacles: ObstacleSet object
    @param taskboundaries: [(min,max)]*3 the shape of the goal set
    if omitted, it puts a unit box around the last state
    """
    if taskboundaries is None:
      taskboundaries = [(-.5,.5)]*3
    return ''.join([self._toAmplStartPoint(
          trajectory[i], (i == len(trajectory)-1 and taskboundaries or None),
          obstacles)
                    for i in range(len(trajectory))])


  def _toAmplStartPoint(self,state,taskboundaries=None,
                       obstacles=None,onlyself=False):
    """
    This function takes a set of joint angles and the values of all of the
    Ampl model variables, returning them as a string

    @param state: a list of joint angles
    @param taskboundaries: if present, this should be in [(min,max)] format,
    if omitted, we assume this task is artificial
    @param obstacles: an ObstacleSet object
    @param onlyself: if True, then only print the length, rel and relz0 vars
    """
    if not onlyself:
      self.ntaskpoints += 1
    elif obstacles:
      print ("This useage is not supported.  If onlyself is True, "
             "the obstacle avoidance constraints are undefined.")
      pdb.set_trace() 

    tb = zeros((3,1,2))
    if taskboundaries is not None:
      taskweights = 1     
      tb[:,0,:] = array(taskboundaries)
    else:
      taskweights = 0
      UB = sum(self.lengths)
      LB = -UB
      tb[:,0,:] = array([(LB,UB)]*3)


    zhat = array([0,0,1])
    yhat = array([0,1,0])
    xhat = array([1,0,0])

    #allocate space for modeling variables (and helpers)
    x = zeros((self.Nlinks+1,3,1))
    J = zeros((self.Nlinks,3,1))
    JJT = zeros((3,3,1))
    jointaxis = zeros((self.Nlinks,3,1,2))
    azimuths = zeros((self.Nlinks,2))
    azimuths[0,0] = 1.0 #arbitrary, non-zero for safety
    origin = zeros((3,2))
    origin[:,0] = dot(self.T0,xhat)
    origin[:,1] = dot(self.T0,yhat)
    # adjacent joint transforms
    Ts = [self.T0] + [Trot(az=self.azimuths[i],
                              incl=state[i])
                          for i in range(self.Nlinks)]
    #origin to joint transforms (i.e. forward kinematics)
    Tsout = [reduce(dot,Ts[:i]) for i in range(1,self.Nlinks+2)]
    for i in range(1,self.Nlinks+1):
      x[i,:,0] = x[i-1,:,0] + dot(Tsout[i],zhat)*self.lengths[i-1]
      jointaxis[i-1,:,0,0] = dot(Tsout[i-1],xhat)
      jointaxis[i-1,:,0,1] = dot(Tsout[i-1],yhat)
    for i in range(1,self.Nlinks):
      azimuths[i,0] = dot(jointaxis[i-1,:,0,0],jointaxis[i,:,0,0])
      azimuths[i,1] = dot(jointaxis[i-1,:,0,1],jointaxis[i,:,0,0])

    assert not any(abs(sum(azimuths**2,1)-1) > .001)
     
    # compute Jacobian and JJT elipsoid
    for i in range(self.Nlinks):
      J[i,:,0] = cross(jointaxis[i,:,0,1],x[-1,:,0]-x[i,:,0])
    JJT[:,:,0] = dot(J[:,:,0].T,J[:,:,0])

    #set up the argument lists
    joints = range(1,self.Nlinks+1)
    tasks = [self.ntaskpoints]
    axes = range(1,4)
    onetwo = range(1,3)
    if onlyself:
      amplvars = [('lengths',self.lengths,[joints]),
                  ('azimuths',azimuths,[joints,onetwo]),
                  ('origin',origin,[axes,onetwo])]
    else:
      amplvars = [('x',x[1:],[joints,axes,tasks]),
                  ('jointaxis',jointaxis,[joints,axes,tasks,onetwo]),
                  ]
    if taskboundaries is not None:
      amplvars += [('J',J,[joints,axes,tasks]),
                  ('JJT',JJT,[axes,axes,tasks]),
                  ('taskbounds', tb, [axes,tasks,onetwo])]
      self.activetasks += tasks
    if obstacles is not None:
      constraints = self.printConstraints(state, obstacles)
    else:
      constraints = ''
    return constraints + ''.join([pyampl.AmplDriver.writeValue(*av)
                                  for av in amplvars])
  

  def stateToPose(self,state):
    """
    Forward kinematics
    @param state, a list of joint angles
    @return positions of the joints in three-space (of the same length)
    """
    zhat = array([0,0,1])
    x = zeros((self.Nlinks+1,3))
    Ts = [self.T0] + [Trot(az=self.azimuths[i],
                              incl=state[i]) 
                          for i in range(self.Nlinks)]
    #origin to joint transforms
    Tsout = [reduce(dot,Ts[:i]) for i in range(1,self.Nlinks+2)]
    for i in range(1,self.Nlinks+1):
      x[i,:] = x[i-1,:] + dot(Tsout[i],zhat)*self.lengths[i-1]
    return x[1:,:]



  def printConstraints(self, state, obstacles):
    """
    @param x
    the particular structure of which is defined elsewhere
    see linearized_nonintersection.mod, smp::rotary.h
    """
    x = zeros((self.Nlinks+1,3))
    x[1:,:] = self.stateToPose(state)
    Nlinks = self.Nlinks
    A = zeros((Nlinks,3,obstacles.Nobs,1))
    b = zeros((Nlinks,obstacles.Nobs,1))

    for i in range(Nlinks):
      A[i,:,:,0],b[i,:,0] = obstacles.findSeparator(x[i,:],x[i+1,:])
    obs = range(1,obstacles.Nobs+1)
    joints = range(1,self.Nlinks+1)
    tasks = [self.ntaskpoints]
    axes = range(1,4)
    amplvars = [('A',A,[joints,axes,obs,tasks]),
                ('b',b,[joints,obs,tasks])]
    return ''.join([pyampl.AmplDriver.writeValue(*av) for av in amplvars])


class ObstacleSet(object):
  """
  THis class is responsible for holding, then printing
  a set of obstacles

  note that this printing is only usefull if we are using the
  integer (or integer-like) formulation involving the obstacles

  If we are using the linearized obstacle constraints, we should call
  findLoosest to figure out what these constraints should be

  param Nobs;
  set obstacles = {1..Nobs};
  param obsts{axes,1..3,obstacles};
  """

  
  def __init__(self,triangles):
    """
    @param mesh, array((3,3,Nobs))
    """
    self.mesh = triangles
    self.Nobs = len(self.mesh[0][0])

  def __str__(self):
    """    
    note that these are params not vars: writeValue might need optional
    args for that
    """
    return (pyampl.AmplDriver.writeValue('Nobs',self.Nobs) + 
            pyampl.AmplDriver.writeValue('obsts',self.mesh,[
          range(1,4),range(1,4),range(1,self.Nobs+1)]))

  def findLoosest(self,x0,x1,retduv=False):
    """
    @param x: array((N,3)) joint locations
    @return (A,b) A: array((3,nobs)),b: array((nobs))
    such that A[:,j][d,u,v] < b[j]
    for x0+d(x1-x0) = a0 + u(a1-a0) + v(a2-a0)
    """
    naxis,ndim,nobs = shape(self.mesh)
    retA,retb = zeros((3,self.Nobs)),zeros(self.Nobs)
    retd = zeros(self.Nobs)
    retu = zeros(self.Nobs)
    retv = zeros(self.Nobs)
    for i in range(self.Nobs):
      a0,a1,a2 = self.mesh[:,:,i]
      (retA[:,i],retb[i], (retd[i],retu[i],retv[i])) = \
          findLoosest(x0,x1,a0,a1,a2,True)
    if retduv:
      return retA,retb,(retd,retu,retv)
    return retA,retb    
    
  def findSeparator(self,x0,x1):
    naxis,ndim,nobs = shape(self.mesh)
    retA,retb = zeros((3,self.Nobs)),zeros(self.Nobs)
    for i in range(self.Nobs):
      a0,a1,a2 = self.mesh[:,:,i]
      retA[:,i],retb[i] = findSeparator(x0,x1,a0,a1,a2)
    return retA,retb


# test code follows using a default manipulator and pre-recorded
# AMPL output
if __name__ == "__main__":
  
  # the default arm has unit lengths and right angle twists
  m = Morphology(Nlinks=8)

  b,w,h = 10,1.5,2 #traj80
  obstacles = [[array([b,w,h]),array([-b,w,h]),array([0,w,b])],
               [array([b,w,h]),array([-b,w,h]),array([0,b,h])],
               [array([b,-w,h]),array([-b,-w,h]),array([0,-w,b])],
               [array([b,-w,h]),array([-b,-w,h]),array([0,-b,h])],
               [array([w,b,h]),array([w,-b,h]),array([w,0,b])],
               [array([w,b,h]),array([w,-b,h]),array([b,0,h])],
               [array([-w,b,h]),array([-w,-b,h]),array([-w,0,b])],
               [array([-w,b,h]),array([-w,-b,h]),array([-b,0,h])]]
  o = ObstacleSet(array(obstacles).transpose((1,2,0)))

  taskbounds = zeros((3,1,2))
  for i in range(3):
    taskbounds[i,0,0] = -1
    taskbounds[i,0,1] = 1
  taskbounds[2,0,:] = [3,5]
  tb = [(-1,1),(-1,1),(1.5,5.5)]

  # This was scraped from a run of standalone_rrtstar_rotary_arm by
  # pyampl.parsing.loadAmplOutputFile
  traj80 = array([[ -0.0669514, -0.526082, -2.52114, -2.10788, -0.986146, 1.19291, -0.573421, 1.54715, ],
  [ 0.680584, 0.0504292, -1.67355, -2.20698, -0.552153, 2.08762, -1.10055, 0.694865, ],
  [ 0.0680809, 0.881161, -1.07144, -1.46794, -0.302133, 1.22488, -0.501548, 0.888032, ],
  [ -0.194962, 1.19568, -1.01061, -1.38181, -0.730083, 1.13772, -0.764716, -0.259215, ],
  [ -0.214959, 0.927878, -0.823715, -1.33969, 0.959288, 1.2091, -0.64841, -0.998987, ],
  [ -0.179889, 0.909993, -0.924341, -1.35635, 1.02285, 1.25653, -0.628009, -1.00324, ],
  [ -0.00782369, 1.06282, -1.0429, -0.0612777, 0.305502, 0.975966, 1.22341, -1.16585, ],
  [ 0.0899392, 1.06559, -1.05859, 0.0437039, 0.346789, 0.96165, 1.32156, -1.26077, ],
  [ 0.808275, 0.871808, -1.54887, 0.246969, 0.742867, 0.765691, 1.32755, -1.26793, ],
  [ 1.22121, -0.0735615, -1.03729, 0.344724, 0.0808869, 0.373615, 0.302209, -0.836744, ],
  [ 1.12387, -0.425992, 0.216282, -0.512312, 0.875609, 0.359851, 0.259246, -1.27194, ],
  [ 1.43182, -0.35081, 0.48709, -0.464127, 0.657647, 0.179177, 0.768875, -1.15306, ],
  [ 1.53685, -0.32459, 0.585822, -0.442846, 0.64436, 0.0901271, 0.719349, -1.19644, ],
  [ 3.14159, 0, 0, 0, 0, 0, 0, 0, ],
  ])

  S = m.amplPrint([traj80[::-1]],[tb],o)
  with open('test.mod','w') as f:
    f.write(S)

