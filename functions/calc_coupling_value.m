function [J] = calc_coupling_value(psi_store, A_t, par, f_proj_fun, a2_fun, a4_fun, kappa_fun, r2_fun, r4_fun, nonlin_fun, kappa_fun_conj, r2_fun_conj, r4_fun_conj)

a2_vals = zeros(size(A_t));
a4_vals = zeros(size(A_t));

for i = 1:numel(a2_vals)
    a2_vals(i) = a2_fun(A_t(i));
    a4_vals(i) = a4_fun(A_t(i));
end

f_proj = f_proj_fun(psi_store,A_t);
kappa  = 1/2/par.m*kappa_fun(A_t);
r2a2   = r2_fun(A_t).*a2_vals;
r4a4   = r4_fun(A_t).*a4_vals;
f_proj_conj = f_proj_fun(psi_store,A_t);
kappa_conj  = 1/2/par.m*kappa_fun_conj(A_t);
r2a2_conj   = r2_fun_conj(A_t).*a2_vals;
r4a4_conj   = r4_fun_conj(A_t).*a4_vals;
nonlin_val = nonlin_fun(psi_store,A_t);


J = f_proj .* (kappa + r2a2 + r4a4) + f_proj_conj .* (kappa_conj + r2a2_conj + r4a4_conj) + nonlin_val;

if all(diff(A_t)>0)
    x_vals = A_t;
else
    x_vals = 1:numel(A_t);
end


figure
subplot(7,1,1)
plot(x_vals,real(f_proj))
grid on
ylabel("f^2")
xlabel("time index")
subplot(7,1,2)
plot(x_vals,real(kappa))
ylabel("kappa")
xlabel("time index")
grid on
subplot(7,1,3)
plot(x_vals,real(r2a2))
ylabel("r2a2")
xlabel("time index")
grid on
subplot(7,1,4)
plot(x_vals,real(r4a4))
ylabel("r4a4")
xlabel("time index")
grid on
subplot(7,1,5)
plot(x_vals,real(nonlin_val))
ylabel("nonlin")
xlabel("time index")
grid on
subplot(7,1,6)
plot(x_vals,real(J))
ylabel("J")
xlabel("time index")
grid on
subplot(7,1,7)
plot(A_t)
ylabel("A")
xlabel("time index")
grid on

end