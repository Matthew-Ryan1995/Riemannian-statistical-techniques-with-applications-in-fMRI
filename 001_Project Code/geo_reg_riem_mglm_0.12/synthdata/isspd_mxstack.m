function t = isspd_mxstack(Y)
t = 0 ;
for i=1:size(Y,3)
    t = t + (~isspd(Y(:,:,i)));
end
t = ~t;