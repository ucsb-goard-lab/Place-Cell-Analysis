function min_corr = subroutine_test_r_HPC(test_vec,data,plot_flag)
%SUBROUTINE_TEST_R_HPC Tests different r_neuropil values to find the value that minimizes the correlation.
%
%   inputs:
%   test_vec: a vector containing the r_neuropil values to test.
%   data: a structure containing the raw_F and neuropil_F data.
%   plot_flag: a binary flag indicating whether to plot the mean correlation (1) or not (0).
%
%   This function tests different r_neuropil values to find the value that minimizes the
%   correlation between the corrected response and the neuropil. It initializes a matrix of
%   correlation coefficients for each combination of cell and r_neuropil. Then, for each r_neuropil
%   value, it calculates the correlation coefficient using the subroutine_find_corr_HPC function and
%   adds it to the matrix. The function finds the r_neuropil value that minimizes the mean correlation
%   and returns a vector of these minimum correlation values for each cell. If plot_flag is set to 1,
%   the function also plots the mean correlation for each r_neuropil value.
%
%   Example:
%       min_corr = subroutine_test_r_HPC([0.1, 0.2, 0.3], data, 1)
%
%   See also ZEROS, LENGTH, SUBROUTINE_FIND_CORR_HPC, ISNAN, MEAN, MIN, PLOT, DISP.
%
%   Written by James Roney for Goard Lab, updated Oct 2016

% Initialize matrix of correlation coefficients of each combination (cell x r_neuropil) 
corr_mat = zeros(size(data.raw_F, 1),length(test_vec));

for i = 1:length(test_vec)
    corr_mat(:,i) = subroutine_find_corr_HPC(test_vec(i),data,0); % save vector from find_corr to one col  
end
corr_mat(isnan(corr_mat)) = 0;
[~,idx] = min(mean(corr_mat,1));
m = test_vec(idx);

if plot_flag == 1
    plot(test_vec, mean(corr_mat,1))
end    

[~,min_idx] = min(corr_mat,[],2);
min_corr = test_vec(min_idx)';
disp(['mean r_neuropil = ' num2str(mean(min_corr))])

