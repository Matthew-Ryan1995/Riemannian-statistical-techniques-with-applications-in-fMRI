function exp_x_y = expmap_sphere(x,y)
%EXPMAP_SPHERE maps a tangent vector y onto the unit sphere.
%
%    exp_x_y = EXPMAP_SPHERE(x,y)
%
%    x, exp_x_y are unit vectors.
%    y is a tangent vector.
%    
%   See also LOGMAP_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:26:13 $

exp_x_y = cos(norm(y))*x+sin(norm(y))*y/zero2one(norm(y));