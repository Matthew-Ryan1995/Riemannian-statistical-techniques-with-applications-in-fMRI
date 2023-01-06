function xbar = karcher_mean_spd(X, W, niter)
%KARCHER_MEAN_SPD calculates the intrinsic mean with weight W on SPD manifolds.
%
%   xbar = KARCHER_MEAN_SPD(X, [], niter)
%   xbar = KARCHER_MEAN_SPD(X, W, niter)
%
%   W is weights.
%   X is a set of points on SPD manifolds.
%   xbar is the Karcher mean of X.
%   niter is the maximum iterations.
%
%   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:20:53 $

xbar = X(:,:,1);

if isempty(W)
    for iter = 1:niter
        phi = mean(logmap_pt2array_spd(xbar,X),3);
        xbar = expmap_spd(xbar, phi);
        if norm(phi) < 1e-18
            break
        end
    end
else
    W = W/norm(W,1);
    for iter = 1:niter
        tmp = logmap_pt2array_spd(xbar,X);
        wtmp = zeros(size(tmp));
        for i = 1:size(tmp,3)
            wtmp(:,:,i) = W(i)*tmp(:,:,i);
        end
        phi = sum(wtmp,3);
        xbar = expmap_spd(xbar, phi);
        if norm(phi) < 1e-18
            break
        end
    end
end
