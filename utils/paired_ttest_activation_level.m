function [sig_t_mtx, p_mtx] = paired_ttest_activation_level(A, B, pval)

%% A and B are a groupwise matrix
% dimension: num_subj x num_roi

t_mtx = zeros(1, size(A,2));
p_mtx = ones(1, size(A,2));
for i = 1:size(A,1)
	x = squeeze(A(:, i)) - squeeze(B(:, i));
%     [h,p,ci,stats] = ttest(x);
    [h,p,ci,stats] = ttest(squeeze(A(:, i)), squeeze(B(:, i)));
    t_mtx(1, i) = stats.tstat;
    p_mtx(1, i) = p;
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


