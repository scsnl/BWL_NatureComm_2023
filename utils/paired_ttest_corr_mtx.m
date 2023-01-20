function [sig_t_mtx, p_mtx] = paired_ttest_corr_mtx(A, B, pval)

%% A and B are a groupwise matrix
% dimension: num_roi x num_roi x num_subj

t_mtx = zeros(size(A,1), size(A,2));
p_mtx = ones(size(A,1), size(A,2));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if i ~= j
            x = squeeze(A(i, j, :)) - squeeze(B(i,j,:));
            [h,p,ci,stats] = ttest(x);
            t_mtx(i, j) = stats.tstat;
            p_mtx(i,j) = p;
            
        end
    end
end

% find FDR threshold
mask=double(p_mtx<=pval);
linkidx = find(mask==1);
pthFDR = FDR(p_mtx(linkidx),pval);
if isempty(pthFDR)
    sig_t_mtx = zeros(size(p_mtx));
else
    sig_t_mtx = (p_mtx <= pthFDR) .* t_mtx .*mask;
end


