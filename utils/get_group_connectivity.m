function [grp_avg_fc, Q, Ci]=get_group_connectivity(pearson_corr_subj_z, numROI, pval)

numRep=1000;

[Conn_avg, Conn_pval]=matrix_element_ttest(pearson_corr_subj_z);
mask_idx = find(Conn_pval<pval);
pval_fdr = FDR(Conn_pval(mask_idx), pval);
Conn_fdr_dir = (Conn_pval <= pval_fdr) .* Conn_avg;
conn_mtx = Conn_fdr_dir;
[Ci, Q] = community_louvain(conn_mtx,[],[],'modularity');
grp_clust_rep = zeros(numROI, numRep);
for irep = 1:numRep
	[Ci_next, Q_next] = community_louvain(conn_mtx,[],[],'modularity');
    grp_clust_rep(:,irep) = Ci_next(:);
    if Q_next > Q
        Q = Q_next; Ci = Ci_next;
    end
end
grp_avg_fc=conn_mtx;






