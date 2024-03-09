function final_corr = subroutine_find_corr_HPC(r_neuropil,data,print_flag)
%SUBROUTINE_FIND_CORR_HPC Finds the correlation coefficient for a specific r_neuropil value.
%
%   finputs:
%   r_neuropil: a scalar specifying the r_neuropil value.
%   data: a structure containing the raw_F and neuropil_F data.
%   print_flag: a binary flag indicating whether to print the average correlation (1) or not (0).
%
%   This function calculates the correlation coefficient for a specific r_neuropil value.
%   It first initializes a vector of correlation coefficients. Then, for each cell, it calculates
%   the squared correlation coefficient and adds it to the vector. The function returns this vector
%   of squared correlation coefficients. If print_flag is set to 1, the function also prints the
%   average correlation.
%
%   Example:
%       final_corr = subroutine_find_corr_HPC(0.5, data, 1)
%
%   See also ZEROS, SIZE, CORRCOEF, MEAN, DISP.
%
%   Written by James Roney for Goard Lab, updated Oct 2016

%Initialize vector of correlation coefficients
test_F = data.raw_F - r_neuropil*data.neuropil_F;
corr = zeros(size(data.raw_F, 1),1);

%Add squared correlation coefficient or each cell
for i = 1:size(data.raw_F, 1)
    j = corrcoef(test_F(i,:), data.neuropil_F(i,:)); % calculate correlation coeff matrix
    corr(i) = j(2).^2; % use j(1,2) to get correlation b/w 2 diferent series
end
final_corr = corr;  

if print_flag ==1
    disp(['avg corr = ' num2str(mean(corr))])
end    

