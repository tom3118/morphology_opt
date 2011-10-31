"""
This file contains helper methods for computing linearized constraints in
morphology.py
"""

def Trot(az=0,incl=0):
  """
  compute rotation matrix:
  """
  return  array([[cos(az), sin(az), 0],
                [-cos(incl)*sin(az),cos(incl)*cos(az),sin(incl)],
                [sin(az)*sin(incl), -sin(incl)*cos(az), cos(incl)]]).T

def findLoosest(x0,x1,a0,a1,a2,retduv=False):
  """  compute the barycentric coordinates
  x0+d(x1-x0) = a0 + u(a1-a0) + v(a2-a0)
  the non-intersection contraints are (d,u,v) < 0, (d,u,v) < 1,
  u+v < 0, and u+v > 1

  the "loosest" is defined as the one which, when scaled by the corresponding
  edge length (e.g. ||x1-x0|| for d) is furthest from the constraint.

  @param retduv: if true, returns (A,b,(d,u,v))

  @return A, b such that A [d,u,v].T < b


  NOTE: this method should be avoided in favor of findSeparator
  """
  dmax = -1 
  try:
    d,u,v = linalg.solve(array([x0-x1,a1-a0,a2-a0]).T,x0-a0)
    if max(map(abs,[d,u,v])) > 1e10:
      raise linalg.linalg.LinAlgError(
        "Poor conditioning: Constraint too sensitive, discarding")
    if d < 0:
      dmax = -d*linalg.norm(x0-x1)
      A,b = array([1,0,0]), 0
    else:
      dmax = (d-1)*linalg.norm(x0-x1)
      A,b = array([-1,0,0]), -1
    if u < 0:
      dm = -u*linalg.norm(a0-a1)
      if dm > dmax:
        dmax = dm
        A,b = array([0,1,0]), 0
    elif u > 1:
      dm = (u-1)*linalg.norm(a0-a1)
      if dm > dmax:
        dmax = dm
        A,b = array([0,-1,0]), -1    
    if v < 0:
      dm = -v*linalg.norm(a0-a2)
      if dm > dmax:
        dmax = dm
        A,b = array([0,0,1]), 0
    elif v > 1:
      dm = (v-1)*linalg.norm(a0-a2)
      if dm > dmax:
        dmax = dm
        A,b = array([0,0,-1]), -1
    if u+v < 0:
      dm = -(u+v)*linalg.norm(a2-a1)
      if dm > dmax:
        dmax = dm
        A,b = array([1,1,0]), 0
    elif u+v > 1:
      dm = (u+v-1)*linalg.norm(a2-a1)
      if dm > dmax:
        dmax = dm
        A,b = array([-1,-1,0]), -1
    if dmax < 0:
      #constraint violated
      A,b = array([0,0,0]),1
  except linalg.linalg.LinAlgError, e: # linalg singular
    # the link is parallel to the plane the return shouldn't matter
    A,b = array([0,0,0]), 1
    d,u,v = 0,0,0
  if retduv:
    return A,b,(d,u,v)
  return A,b

def findSeparator(x0,x1,a0,a1,a2):
  """
  find a plane separating the segment (x0,x1) and the triangle (a0,a1,a2)
  This finds an optimal separator and that touches the triangle.
  
  There are 5 cases:
  seg,  tri
  \hline
  pt,   pt
  pt,   edge
  pt,   plane
  edge, pt
  edge, edge

  @return A,b such that Ax = b is the separator, Ax_0 <= b and
  Aa_i >= b holding with equality at at least one i
  """
  # check pt,pt
  seps = reduce(list.__add__,[[(x,a) for x in x0,x1] for a in a0,a1,a2])
  dists = array([linalg.norm(x-a) for x,a in seps])
  mindist = min(dists)
  x,a = seps[argmin(dists)]
  A,b = a-x,dot(a-x,a) 

  # check pt, edge
  for t1,t2 in [(a0,a1),(a1,a2),(a2,a0)]:
    for x in x0,x1:
      if linalg.norm(cross(x-t1,t2-t1))/linalg.norm(t2-t1) < mindist:
        if 0 < dot(t2-t1,x-t1)/(linalg.norm(t2-t1)**2) < 1:
          a = t1+(t2-t1)*dot(t2-t1,x-t1)/(linalg.norm(t2-t1)**2)
          assert not abs(dot(a-x,t1-t2)) > .01
          temp = A,b,mindist
          A,b = a-x,dot(a-x,a)
          mindist = linalg.norm(cross(x-t1,t2-t1))/linalg.norm(t2-t1)

 # check pt, plane
  for x in x0,x1:
    if abs(dot(x-a0,cross(a1-a0,a2-a0))/
           linalg.norm(cross(a1-a0,a2-a0))) < mindist:
      n = cross(a1-a0,a2-a0)
      n /= linalg.norm(n)
      a = x + n*dot(a0-x,n)
      # check that it is in the triangle by computing barycentric coords
      # (the uu,vv mess is just a fast inverse)
      uu,uv,vv = (dot(a1-a0,a1-a0),dot(a2-a0,a1-a0),dot(a2-a0,a2-a0))
      denom = 1.0*(uu*vv - uv**2)
      u = (vv*dot(a-a0,a1-a0) - uv*dot(a-a0,a2-a0))/denom
      v = (uu*dot(a-a0,a2-a0) - uv*dot(a-a0,a1-a0))/denom
      assert not u+v == 0
      if 0 < v and 0 < u and u+v < 1:
        A,b = a-x,dot(a-x,a)
        mindist = linalg.norm(A)

  # check edge, pt (reverses pt, edge but A is negated and b differs)
  t1,t2 = x0,x1
  for x in a0,a1,a2:
    if linalg.norm(cross(x-t1,t2-t1))/linalg.norm(t2-t1) < mindist:
      #check that it is on the segment
      if 0 < dot(t2-t1,x-t1)/(linalg.norm(t2-t1)**2) < 1:
        a = t1 + (t2-t1)*dot(t2-t1,x-t1)/(linalg.norm(t2-t1)**2)
        A,b = x-a,dot(x-a,x)
        mindist = linalg.norm(cross(x-t1,t2-t1))/linalg.norm(t2-t1)

  # check edge, edge 
  # for explanation see
  # http://2000clicks.com/mathhelp/GeometryPointsAndLines3D.aspx
  for t0,t1 in [(a0,a1),(a1,a2),(a2,a0)]:
    if (abs(dot(x0-t0,cross(x1-x0,t1-t0)))
        /linalg.norm(cross(x1-x0,t1-t0))) < mindist:
      # check that the closest point is on both segments
      u = cross(x1-x0,t1-t0)/linalg.norm(cross(x1-x0,t1-t0))
      g = dot(x0-t0,u)
      ts =  linalg.lstsq(array([x0-x1,t1-t0]).T,g*u + x0-t0)[0]
      if 0 < ts[0] < 1 and 0 < ts[1] < 1:
        x = x0 + ts[0]*(x1-x0)
        a = t0 + ts[1]*(t1-t0)
        A,b = a-x,dot(a-x,a)
        mindist = abs(g)

  return A,b
