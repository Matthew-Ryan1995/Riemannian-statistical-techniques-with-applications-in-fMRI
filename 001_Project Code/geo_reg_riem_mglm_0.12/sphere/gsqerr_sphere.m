function gsr = gsqerr_sphere(X, X_hat)
%GSQERR_SPHERE returns the sum of geodesic squared error on the unit sphere.
%
%    gsr = GSQERR_SPHERE(X, X_hat)
%
%
%   X, X_hat are a set of points on the unit sphere (dim_X x  N, where dim_X = size(X,1)).
%
%   See also FEVAL_SPHERE, R2STAT_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

[ndim, ndata] = size(X);
tanvec = zeros(ndim, ndata);
for idata = 1:ndata
    tanvec(:,idata) = logmap_sphere(X(:,idata), X_hat(:,idata));
end
%tanvec
gsr = trace(tanvec'*tanvec);