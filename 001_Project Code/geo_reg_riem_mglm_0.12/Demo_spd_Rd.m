% DEMO on SPD manifolds
clear;
addpath(genpath(pwd))
rng(222);
disp('Start.')
for i = 1:10
    seed = randi(1000);
    rng(seed);
    
    synth_dti_data

    Ybar = karcher_mean_spd(Y,[],500);
    [p, V, E, Yhat, gnorm] = mglm_spd(X,Y);
    r2_iterative  = r2stat_spd(Ybar, Y, Yhat);

    [ple, Vle, Ele, Yhatle, logY, Yv ] = mglm_logeuc_spd(X,Y);
    r2_logeuc  = r2stat_spd(Ybar, Y, Yhatle);
    
    [ple_Rd, Vle_Rd, Ele_Rd, Yhatle_Rd, logY_Rd, Yv_Rd ] = mglm_logeuc_spd_Rd(X,Y);
    r2_logeuc_Rd  = r2stat_spd(Ybar, Y, Yhatle_Rd);
    
    fprintf('seed %-4d: r2_iterative %f r2_logeuc %f r2_logeuc_Rd %f\n', ...
        seed, r2_iterative, r2_logeuc, r2_logeuc_Rd)
end
