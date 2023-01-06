function Vnew = invembeddingRd_vecs(p,Vecs)
%EMBEDDINGR6_VECS embeds a set of tangnet vectors V in TpM onto R6. This transports V from TpM to T_{I}M.
%
%   This will be replaced later with the group action for n x n spd matrix and a vectorization. 
%    
%   Vnew = EMBEDDINGR6_VECS(p,V)
%
%   See also EMBEDDINGR6, LOGMAP_SPD

%   $ Hyunwoo J. Kim $  $ 2016/04/20 15:04:11 (CDT) $
%   $ Revision: 0.12 $  

    VatI = vec2symmx(Vecs);
    Vnew = grpaction_i2p(p, VatI);
end