function [p, V, E, Y_hat, gnorm] = mglm_sphere(X, Y, varargin)
%MGLM_SPHERE performs MGLM on the unit sphere by interative method.
%
%   [p, V, E, Y_hat, gnorm] = MGLM_SPHERE(X, Y)
%   [p, V, E, Y_hat, gnorm] = MGLM_SPHERE(X, Y, MAXITER)
%   has optional parameter MAXITER.  
%
%   The result is in p, V, E, Y_hat.
%
%   X is dimX x N column vectors
%   Y is dimY x N column vectors (points on the unit sphere in R^dimY).
%   p is a base point.
%   V is a set of tangent vectors (dimY x dimX).
%   E is the history of the sum of squared geodesic error.
%   Y_hat is the prediction.
%   gnorm is the history of norm of gradients.
%
%   See also MGLM_LOGEUC_SPHERE, EXPMAP_SPHERE, FEVAL_SPHERE, PREDICTION_SPHERE,
%   LOGMAP_VECS_SPHERE, PARALLELTRANSLATEATOB_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:27:37 $

    
    ndimX = size(X,1);
    ndimY = size(Y,1);
    
    p = Y(:,1);
    V = zeros(ndimY,ndimX); % Random initialization
    V = proj_TpM(p,V); % To get the valid tangent vectors.
    
    if nargin >=3
        maxiter = varargin{1};
    else
        maxiter = 5000;
    end
    
    % Gradient Descent algorithm
    % Step size
    c1 = 1;
    
    % Safeguard parameter
    c2 = 1;

    E = [];
    gnorm = [];
    E = [E; feval_sphere(p,V,X,Y)];
    step = c1;
    for niter=1:maxiter
        Y_hat = prediction_sphere(p,V,X);
        J = logmap_vecs_sphere(Y_hat,Y);
        err_TpM = paralleltranslateAtoB_sphere(Y_hat, p, J);
        gradp = -sum(err_TpM,2);
        
        % v projection on to tanget space
        gradV = zeros(size(V));
        for iV = 1:size(V,2)
            gradV(:,iV) = -err_TpM*X(iV,:)';
        end
        gnorm_new = norm([gradV gradp]);
        
        % safeguard
        [gradp gradV] = safeguard(gradp, gradV,c2);
        
        moved = 0;
        for i = 1:200
            step = step*0.5;
            % Safegaurd for gradv, gradp
            p_new = unitvec(expmap_sphere(p,-step*gradp));
            V_new = V -step*gradV;
            V_new = paralleltranslateAtoB_sphere(p,p_new,V_new);
                        
            E_new = feval_sphere(p_new, V_new, X, Y);
            if E(end) > E_new
                p = p_new;
                V = V_new;
                E = [E; E_new];
                gnorm = [gnorm; gnorm_new];
                moved = 1;
                step = min(step*2,1);
                break
            end
            if step < 1e-20
                break
            end
        end
        if moved ~= 1
            break
        end
    end

    E = [E; feval_sphere(p,V,X,Y)];
    Y_hat = prediction_sphere(p,V,X);
end

function V = proj_TpM(p,V)
    for i = 1:size(V,2)
        v = V(:,i);
        V(:,i) = v-p'*v*p;
    end
end 

function [gradp, gradV] = safeguard(gradp, gradV, c2)
    vecnorms = @(A) sqrt(sum(A.^2,1));
    norms = [ vecnorms(gradV) norm(gradp)];
    maxnorm = max(norms);
    if maxnorm > c2
        gradV = gradV*c2/maxnorm;
        gradp = gradp*c2/maxnorm;
    end
    
end
