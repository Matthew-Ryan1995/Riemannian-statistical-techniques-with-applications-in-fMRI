function w_new = paralleltranslateAtoB_spd(a, b, w)
%PARALLELTRANSLATEATOB_SPD transports a set of tangent vectors w from TaM to
%TbM.
%
%   w_new = PARALLELTRANSLATEATOB_SPD(a, b, w)
%
%   a, b are points on SPD matrices. 
%   w is a set of tangent vectors.
%   w_new is a set of transported tangent vectors.
%
%   See also MGLM_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $

    if size(a,3) < size(b,3)
        a = repmat(a, [1 1 size(b,3)]);
    elseif size(a,3) > size(b,3)
        b = repmat(b, [1 1 size(a,3)]);
    end

    if size(b,3) ~= size(w,3)
        % a, b are fixed
        % This changes only w.
        fixab = 1;
        P1 = a;
        P2 = b;
    else
        fixab = 0;
    end

    w_new = zeros(size(w));
    for i = 1:size(w,3)
        if fixab == 0
            P1 = a(:,:,i);
            P2 = b(:,:,i);
        end

        if norm(P1-P2) < 1e-20
            w_new(:,:,i) = w(:,:,i);
            continue
        end
        
%       invP1 = inv(P1);
%       P12 = sqrtm(invP1*P2*P2*invP1);
%       T12 = P12\invP1*P2;
%       B = P1\w(:,:,i);
%       w_new(:,:,i) = P2*T12'*B*T12;
        w_new(:,:,i) = parallel(P1,P2,w(:,:,i));
        
    end

    %% symmetrization.
    for i = 1:size(w,3)
        w_new(:,:,i) = (w_new(:,:,i)+w_new(:,:,i)')/2;
    end
end

function w_new = parallel(p,q,w)
    rtp = sqrtm(p);
    invrtp = inv(rtp);
    v = logmap_spd(p,q);
    r = expm(invrtp*v/2*invrtp);
    w_new = rtp*r*invrtp*w*invrtp*r*rtp;
end
