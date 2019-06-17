function [x,x_coeff]= OFlowSparseGenerate02(nOFlow,tau_max_x,S)
% Generate O-flows that are sparse under DCT transform
%   x is of the size nOFlow-by-tau_max_x
%   The DCT coefficients of each row of x has S+1 nonzeros
%       where the extra 1 nonzero is from DC component (x is non-negative)
if S<=0
    error('S should be positive integer');
elseif S>floor(tau_max_x/4)
    error('S is too large: S should be <= tau_max_x/4');
end

epsilon = 8e-2;
x_coeff = zeros(nOFlow,tau_max_x);
for nRow = 1:nOFlow
    % Take random S positions. Exclude the first position from these S positions 
    coeff_t = randn(1,S);
    coeff_t = sign(coeff_t)*2 + coeff_t;
    x_coeff(nRow,randsample(floor(tau_max_x/2),S)+1) = coeff_t;
end

x = idct(x_coeff');
x = x';

% make x non-negative
x_min = min(x,[],2); 
x_min_index = x_min < epsilon;
x(x_min_index,:) = x(x_min_index,:) ...
    - repmat(x_min(x_min_index),1,tau_max_x) + epsilon;
x_coeff = dct(x');
x_coeff = x_coeff';

end
