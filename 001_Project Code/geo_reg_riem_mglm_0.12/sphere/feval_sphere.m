function E = feval_sphere(p,V,X,Y)
%FEVAL_SPHERE evaluates the objective function value (the sum of squared geodesic errors) of MGLM on the unit sphere. 
%
%    E = FEVAL_SPHERE(p, V, X, Y)
%
%   !! make sure that X is centered if p, V are calculated by centered X !!
%
%   X is covariates, a set of column vectors. (dimX X N, where dimX = size(X,1))
%   Y is responses, a set of points on the unit sphere. (dim_Y x  N, where dim_Y = size(Y,1)).
%   p is the base point.
%   V is a set of tangent vectors (dimY x dimX).
%   E is the sum of squared geodesic errors.
%
%   See also LOGMAP_SPHERE, EXPMAP_SPHERE, GSQERR_SPHERE, PREDICTION_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:42:13 $

P_hat = prediction_sphere(p,V,X);
E = gsqerr_sphere(P_hat, Y);