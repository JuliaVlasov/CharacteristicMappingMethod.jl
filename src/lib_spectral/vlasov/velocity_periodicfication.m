function [v_periodic,sigma] = velocity_periodicfication(params, bump_transition_width)

if (nargin <2 || isempty(bump_transition_width))
    b = 0.2*params.L(2);
else 
    b = bump_transition_width;
end

%%% first calculate "vorticity" such that global vorticity omega = 0 
% (u_1,u_2) = (g(v), d/dx \phi)
% omega = d/dx u_2 - d/dv u_1 = d^2/(dx)2 \phi(x) - d/dv g(v) == 0 
% assume phi(x)=0 
% if g(v) = v 
% omega = d/dv g(v) = 1 != 0 ):
% so we do the following:
% look for a function h(v) = d/dv g(v) that keeps global omega=0, but is 1 close to the
% origin and negative on the boundaries (= heaviside function)
a = 0.4;
sigma = 0.5 - 0.5 * tanh(2*pi*(abs(params.v)-a*params.L(2))/(b));
h = sigma;%tanh((-abs(params.v)+a*params.L(2)/2)/(b))+0.5;
h = h-mean(h);
h = h./max(h, [], 'all');

kv = params.kv;
kv2 = kv.^2;
kv2(1)=1;
intu_hat =(-fft(h)./kv2);
intu = ifft(intu_hat,"symmetric");
v_periodic = real(ifft(intu_hat.*kv,"symmetric"));