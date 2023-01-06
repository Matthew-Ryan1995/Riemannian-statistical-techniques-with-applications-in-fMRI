function Yhat = prediction_sphere(p,V,X)
%PREDICTION_SPHERE predicts Yhat based on estimate p, V and covariate X.
%
%   See also MGLM_SPHERE, FEVAL_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:42:13 $

ndata = length(X);
Yhat = zeros(size(p,1),ndata);

for i = 1:ndata
    Yhat(:,i) = expmap_sphere(p,V*X(:,i));
end
