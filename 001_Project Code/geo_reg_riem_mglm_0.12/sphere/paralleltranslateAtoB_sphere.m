function w_new = paralleltranslateAtoB_sphere(a, b, w)
%PARALLELTRANSLATEATOB_SPHERE transports a set of tangent vectors w from TaM to
%TbM.
%
%   w_new = PARALLELTRANSLATEATOB_SPHERE(a, b, w)
%
%   a, b are points on the unit sphere. 
%   w is a set of tangent vectors.
%   w_new is a set of transported tangent vectors.
%
%   See also MGLM_SPHERE, PARALLELTRANSLATE_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $


if size(a,2) < size(b,2)
    a = a*ones(1,size(b,2));
elseif size(a,2) > size(b,2)
    b = b*ones(1,size(a,2));
end

v = zeros(size(a));
for i = 1:size(a,2)
    v(:,i) = logmap_sphere(a(:,i),b(:,i));
end

w_new = paralleltranslate_sphere(a, v, w);

