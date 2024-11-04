function phi_fft = give_potential(params,f)


    % laplacian is division -|k|^2
    K = -params.kx'.^2; 
    K ( abs(K) < 1.0e-11 ) = 1;
    % to avoid a division by zero, we set the zeroth wavenumber to one.
    % this leaves it's respective Fourier coefficient unaltered, so the
    % zero mode of Sk is conserved. This way, Sk's zero mode implicitly
    % defined the zero mode of the result
    % Note that the zero mode is NOT uniquely defined: in a periodic
    % setting, the solution of Laplace's (or Poisson's) equation is only
    % defined up to a constant! You can freely overwrite the zero mode,
    % therefore.
    f_int = params.dv*sum(f,2); % trapez rule for periodic system is just the sum
    b = fft(params.l*f_int-1); % convert FFT( 1-l \int f dv)
    phi_fft = b ./ K; % solves second equation of vlassov poisson

end