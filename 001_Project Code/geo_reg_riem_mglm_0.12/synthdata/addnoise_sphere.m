function Anew = addnoise_sphere(A, maxerr,varargin)
if nargin ==3
    c = varargin{1};
else 
    c =1;
end
V = randn(size(A))*c;
V = V-V'*A*A;
if norm(V) > maxerr
    V = V/norm(V)*maxerr;
end
Anew = expmap_sphere(A, V);
