function [second_deriv] = calculate_second_deriv(f_t,diff_filter,n_order,dv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate second derivative numerically removing linear components (setting both integral constants to 0)
% 
% Args:
%     f_t (array):                          function value
%     diff_filter (array):                  differentiation filter
%     n_order (int):                        order of differentiation filter
%     dv (double):                          discretization size
%
% Returns:
%     second_deriv (array):                 approximated second derivative values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_t_zeroed      = f_t - linspace(f_t(1),f_t(2),numel(f_t));

temp = imfilter([f_t_zeroed(1)-(2*n_order:-1:1)*(f_t_zeroed(2)-f_t_zeroed(1)),f_t_zeroed,f_t_zeroed(end)+(1:2*n_order)*(f_t_zeroed(end)-f_t_zeroed(end-1))],diff_filter,'replicate','same');

second_deriv    = temp(2*n_order+1:end-2*n_order)/(dv^2);
second_deriv    = second_deriv - linspace(second_deriv(1),second_deriv(end),numel(second_deriv));
second_deriv(abs(second_deriv) < 1e-7) = 0;

end 