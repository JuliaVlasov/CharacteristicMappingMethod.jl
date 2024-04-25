function [Mjet, M] = RK4_jetstages(params, sgsz, u1j_fun, u2j_fun, t, dt, mgsz)

% output maps corresponding to evaluation locations of all stages
% optionally output evaluation of the velocities at those maps.

Uj_fun = @(s) cat(3, u1j_fun(s), u2j_fun(s));

ngm = prod(mgsz);

M0 = zeros(ngm, 4, 2);
U0 = HCompose(params, sgsz, Uj_fun(t), mgsz, M0);

M1 = 0.5*dt*U0;
U1 = HCompose(params, sgsz, Uj_fun(t+0.5*dt), mgsz, M1);

M2 = 0.5*dt*U1;
U2 = HCompose(params, sgsz, Uj_fun(t+0.5*dt), mgsz, M2);

M3 = dt*U2;
U3 = HCompose(params, sgsz, Uj_fun(t+dt), mgsz, M3);


Mjet = dt*(U0 + 2*U1 + 2*U2 + U3)/6;

M =  cat(4, M0, M1, M2, M3);

return