function V = logmap_vecs_sphere(xi,X)
%LOGMAP_VECS_SPHERE calculates logmap(xi,X).
%
%    v = LOGMAP_SPHERE(xi,X)
%
%    xi is the base point.
%    X are a set of unit vectors (points on the unit sphere).
%    v is a tangent vector.
%    
%   See also LOGMAP_SPHERE, EXPMAP_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 18:16:13 $

if size(xi,2) == size(X,2)
    V = zeros(size(X));
    for j = 1:size(X,2)
        V(:,j) = logmap_sphere(xi(:,j),X(:,j));
    end
else
    V = zeros(size(X));
    for j = 1:size(X,2)
        xj = X(:,j);
        V(:,j) = logmap_sphere(xi,xj);
    end
end

