function [d,sigma,x0] = calculate_difference_from_fit(x,psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates interwell distance, sigma and center of mass from psi
% 
% Args:
%     x (array):                            spatial values
%     psi (array):                          N_x x N_t wavefunctions of psi
%
% Returns:
%     d (array):                            N_x x 1 interwell distance values
%     sigma (array):                        N_x x 1 sigma in situ values
%     x0 (array):                           N_x x 1 center of mass values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[N_x,N_t]   = size(psi);


x0                      = trapz(x,x.*abs(psi).^2,1);
x_grid                  = repmat(x,[1,N_t]);
psi_plus                = psi;
psi_plus(x_grid<repmat(x0,[N_x,1]))     = 0;
psi_minus               = psi;
psi_minus(x_grid>repmat(x0,[N_x,1]))    = 0;
mue_plus                = 2*trapz(x,x.*abs(psi_plus).^2,1);
mue_minus               = 2*trapz(x,x.*abs(psi_minus).^2,1);

d               = mue_plus - mue_minus;
sigma_1         = abs(trapz(x,psi))./(max(abs(psi),[],1))/(pi*sqrt(2));
sigma           = [sigma_1];

sse             = zeros(size(d));

delta_disc      = 5;
ln_phi          = log(abs(psi));
ln_phi_plus     = log(abs(psi_plus));
ln_phi_minus    = log(abs(psi_minus));

[~,maxindx_plus]    = max(ln_phi_plus,[],1);
[~,maxindx_minus]   = max(ln_phi_minus,[],1);

d_ln            = d;
sigma_ln        = sigma;

for indx = 1:N_t
    x_taken         = x(maxindx_plus(indx)-delta_disc:maxindx_plus(indx)+delta_disc);
    psi_taken       = ln_phi(maxindx_plus(indx)-delta_disc:maxindx_plus(indx)+delta_disc,indx);
    poly_par        = polyfit(x_taken,psi_taken,2);
    par_plus        = [-poly_par(2)/2/poly_par(1),1/sqrt(-2*poly_par(1))];

    x_taken         = x(maxindx_minus(indx)-delta_disc:maxindx_minus(indx)+delta_disc);
    psi_taken       = ln_phi(maxindx_minus(indx)-delta_disc:maxindx_minus(indx)+delta_disc,indx);
    poly_par        = polyfit(x_taken,psi_taken,2);
    par_minus        = [-poly_par(2)/2/poly_par(1),1/sqrt(-2*poly_par(1))];

    d_ln(indx)      = par_plus(1)-par_minus(1);
    sigma_ln(indx)  = 1/2*(par_plus(2)+par_minus(2));
end

d = d_ln;
sigma = sigma_ln;

end

