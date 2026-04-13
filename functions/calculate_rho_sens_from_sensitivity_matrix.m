function [I_CL,I_Sw] = calculate_rho_sens_from_sensitivity_matrix(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates collinearity index and factor of max to min eigenvalue from Sensitivity matrix
% 
% Args:
%     S (array):                            number of measurements x number of parameters sensitivity matrix
%
% Returns:
%     I_CL (double):                        collinearity index (1/minimum eigenvalue)
%     I_Sw (double):                        factor of max to min eigenvalue from Sensitivity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_w         = S./vecnorm(S,2,1);                % normalizing collumns
S_w_quad    = S_w'*S_w;                         % quadratic matrix

eigv_S_w    = abs(eig(S_w_quad));               % absolute eigenvalue
I_Sw        = max(eigv_S_w) / min(eigv_S_w);    % factor of max to min eigenvalue
I_CL        = 1 / sqrt(min(eigv_S_w));          % collinearity index

end