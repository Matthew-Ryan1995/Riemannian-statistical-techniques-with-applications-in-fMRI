% DEMO on the unit sphere.
clear;
addpath(genpath(pwd));
rng(222)
disp('Start.');
for i = 1:10
    seed = randi(1000);
    rng(seed);
    
    synth_sphere_data
    
    Ybar = karcher_mean_sphere(Y,[],500);
    [p, V, E, Yhat, gnorm] = mglm_sphere(X,Y);
    r2_iterative = r2stat_sphere(Ybar,Y,Yhat);

    [ple, Vle, Ele, Yhatle, U ] = mglm_logeuc_sphere(X,Y);
    r2_logeuc = r2stat_sphere(Ybar,Y,Yhatle);
    
    fprintf('seed %-4d: r2_iterative %f r2_logeuc %f\n', seed, r2_iterative, r2_logeuc);
end
