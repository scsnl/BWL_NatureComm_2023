function [AvsB_avg, AvsB_pval] = matrix_element_ttest_pair(A, B)

%% A and B are a groupwise matrix
% dimension: num_roi x num_roi x num_subj
%%

AvsB_avg = zeros(size(A,1), size(A,2));
AvsB_pval = ones(size(A,1), size(A,2));

for i = 1:size(A,1)
    for j = 1:size(A,2)
        if i ~= j
            x = squeeze(A(i, j, :)) - squeeze(B(i,j,:));
            AvsB_avg(i, j) = mean(x);
            [~, ij_pval] = ttest(x);
            AvsB_pval(i,j) = ij_pval;
        end
    end
end

end