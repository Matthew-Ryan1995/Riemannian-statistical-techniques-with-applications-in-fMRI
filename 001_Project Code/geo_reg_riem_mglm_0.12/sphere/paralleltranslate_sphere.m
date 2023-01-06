function W_new = paralleltranslate_sphere(p,V,W)
%PARALLELTRANSLATE_SPHERE transports a tangent vector w from TpM along V.
%
%   w_new = PARALLELTRANSLATE_SPHERE(p, V, w)
%
%   a, b are points on the unit sphere. 
%   w is a set of tangent vectors.
%   w_new is a set of transported tangent vectors.
%
%   See also MGLM_SPHERE, PARALLELTRANSLATEATOB_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $

if size(V,2) < size(W,2)
    V = V*ones(1,size(W,2));
end

vnorm = sqrt(diag(V'*V));
scale_factor = diag(V'*W)./vnorm.^2;
scale_factor(isnan(scale_factor)) = 0;
W_orth = W - V*diag(scale_factor);
if size(p,2) == 1
    V_par = p*((-sin(vnorm)).*vnorm)' + V*diag(cos(vnorm));
else
    V_par = p*diag((-sin(vnorm)).*vnorm) + V*diag(cos(vnorm));
end
W_new = V_par*diag(scale_factor) + W_orth;
