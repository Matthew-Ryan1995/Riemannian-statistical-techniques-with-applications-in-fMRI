function X = expmap_vecs_sphere(xi,V)
%EXPMAP_VECS_SPHERE calculats expmap(xi,V) for the unit sphere.
%
%   X = EXPMAP_VECS_SPHERE(xi, V)
%
%   X is dimX x N column vectors
%   xi is the base point (unit vector). 
%   V is dimY x N. N tangent vectors.
%
%   See also EXPMAP_VECS_SPHERE, EXPMAP_SPHERE, MGLM_SPHERE, MGLM_LOGEUC_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:37:37 $

X = zeros(size(V));
for j = 1:size(V,2)
    X(:,j) = expmap_sphere(xi,V(:,j));
end