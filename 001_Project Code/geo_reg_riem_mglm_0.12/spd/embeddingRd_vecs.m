function Vnew = embeddingRd_vecs(p,V)
%EMBEDDINGR6_VECS embeds a set of tangnet vectors V in TpM onto R6. This transports V from TpM to T_{I}M.
%
%   This will be replaced later with the group action for n x n spd matrix and a vectorization. 
%    
%   Vnew = EMBEDDINGRD_VECS(p,V)
%
%   See also EMBEDDINGR6_VECS, LOGMAP_SPD

%   $ Hyunwoo J. Kim $  $ 2016/04/20 14:59:24 (CDT) $
    VwrI = grpaction_p2i(p, V);
    Vnew = symmx2vec(VwrI);
end