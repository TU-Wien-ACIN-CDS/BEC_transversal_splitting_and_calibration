function [k,momentum_carpet] = calc_moment(x,psi,N_window,dim_x,k_thresh)

[N_x,N_t]       = size(psi);

k_window        = 2*pi*(1/(x(2)-x(1))/N_window)*(-N_window/2:N_window/2-1)';
selector        = abs(k_window)<k_thresh;
N_selected      = sum(selector);
k               = k_window(selector);

momentum_carpet = zeros(N_selected,N_t);

for iter_t = 1:N_t
    Y                           = fftshift(fft(psi(:,iter_t),N_window,dim_x));
    momentum_carpet(:,iter_t)   = Y(selector) / N_window;
end

end