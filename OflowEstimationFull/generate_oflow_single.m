function [x]= generate_oflow_single(n_t,n_o)
%%%x is of the size nOFlow-by-tau_max_x
%%%x is non-negative
x=rand(1,n_o*n_t);
%%%make non-negative
x=x-min(x);
x=reshape(x,n_o,n_t);
end