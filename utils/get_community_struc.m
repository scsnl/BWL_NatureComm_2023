function [Q, Ci]=get_community_struc(conn_mtx)

numRep=100;

[Ci, Q] = modularity_louvain_und_sign(conn_mtx);
for irep = 1:numRep
	[Ci_next, Q_next] = modularity_louvain_und_sign(conn_mtx);
    if Q_next > Q
        Q = Q_next; Ci = Ci_next;
    end
end






