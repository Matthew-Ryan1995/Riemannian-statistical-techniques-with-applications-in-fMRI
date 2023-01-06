function r2  = r2stat_sphere(Y_bar, Y, Y_hat)
%R2STAT calculates R2 statistic.
%
%    r2  = r2stat(Y_bar, Y, Y_hat)
%
%    Example:
%        
%        Y_bar = karcher_mean_sphere(Y,[],1000);
%        r2  = r2stat_sphere(Y_bar, Y, Y_hat, logmap);
%
%   See also GSQERR_SPHERE, FEVAL_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

gvar = gsqerr_sphere(repmat(Y_bar,1,size(Y,2)), Y);
uvar = gsqerr_sphere(Y, Y_hat);
r2 = 1-uvar/gvar;
