function v = logmap_sphere(xi,xj)
%LOGMAP_SPHERE calculates logmap(xi,xj).
%
%    v = LOGMAP_SPHERE(xi,xj)
%
%    xi, xj are unit vectors (points on the unit sphere).
%    v is a tangent vector.
%    
%   See also EXPMAP_SPHERE, LOGMAP_VECS_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:32:13 $

v = (xj-xi'*xj*xi)/zero2one(sqrt(1-(xi'*xj)^2))*zero2one(acos(xi'*xj));