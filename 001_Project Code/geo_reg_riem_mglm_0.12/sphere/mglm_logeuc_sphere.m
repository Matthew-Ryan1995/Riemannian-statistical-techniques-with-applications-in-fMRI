function [p, V, E, Y_hat, U] = mglm_logeuc_sphere(X,Y,varargin)
%MGLM_SPHERE performs MGLM on the unit sphere by by Log Euclidean
%framework.
%
%   [p, V, E, Y_hat, U] = MGLM_LOGEUC_SPHERE(X, Y)
%
%   [p, V, E, Y_hat, gnorm] = MGLM_SPHERE(X, Y, niter)
%   has optional parameter niter for Karcher mean calculation.  
%
%   The result is in p, V, E, Y_hat.
%
%   X is dimX x N column vectors
%   Y is dimY x N column vectors (points on the unit sphere in R^dimY).
%
%   p is a base point.
%   V is a set of tangent vectors (dimY x dimX).
%   E is the sum of squared geodesic error.
%   Y_hat is the prediction.
%   U is the orthogonal basis of TpM.
%
%   See also MGLM_SPHERE, EXPMAP_SPHERE, FEVAL_SPHERE, PREDICTION_SPHERE,
%   LOGMAP_VECS_SPHERE, PARALLELTRANSLATEATOB_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:50:37 $

if nargin >=3
    niter = varargin{1};
else
    niter = 500;
end
% Linear transform

[ ndim ndata] = size(X);
Xc = X - repmat(mean(X,2),1,ndata);
p = karcher_mean_sphere(Y, ones(ndata,1)/ndata, niter);
logY = logmap_vecs_sphere(p, Y);

% Get orthogonal bases of TpM
U = null(ones(size(p,1),1)*p');


Yu = U'*logY; %logY is represented by U
% Yu = L*X
L = Yu/Xc;
V = U*L;
logY_hat = V*Xc;
Y_hat = expmap_vecs_sphere(p,logY_hat);
E = gsqerr_sphere(Y, Y_hat); % Geodesic squared error