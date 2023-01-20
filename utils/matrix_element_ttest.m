function [A_avg, A_pval] = matrix_element_ttest(A)

%% A is a groupwise matrix
% dimension: num_roi x num_roi x num_subj
%%

A_avg = zeros(size(A,1), size(A,2));
A_pval = ones(size(A,1), size(A,2));

for i = 1:size(A,1)
    for j = 1:size(A,2)
        if i ~= j
            x = squeeze(A(i, j, :));
            A_avg(i, j) = mean(x);
            [~, ij_pval] = ttest(x);
            A_pval(i,j) = ij_pval;
        end
    end
end

end