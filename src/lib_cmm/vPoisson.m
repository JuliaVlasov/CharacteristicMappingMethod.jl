function [dphi_dx_h,kx,ky] = vPoisson(f, defKFOp, gsz, dom,dv)
l=1;
[kx, ky] = genWaveNumbers(dom, gsz);

% technically curl of negative laplacian, cancel the signs)


    % laplacian is division -|k|^2
    K2 = kx(1,:)'.^2; 
    K2 ( abs(K2) < 1.0e-11 ) = 1;
    % to avoid a division by zero, we set the zeroth wavenumber to one.
    % this leaves it's respective Fourier coefficient unaltered, so the
    % zero mode of Sk is conserved.dphi_dx_h = 1i*phi_fft.*kx(1,:); This way, Sk's zero mode implicitly
    % defined the zero mode of the result
    % Note that the zero mode is NOT uniquely defined: in a periodic
    % setting, the solution of Laplace's (or Poisson's) equation is only
    % defined up to a constant! You can freely overwrite the zero mode,
    % therefore.
    f_int = dv*sum(f,1); % trapez rule for periodic system is just the sum
    b = fft(l*f_int-1); % convert FFT( 1-l \int f dv)
    b(1) = 0;
    phi_fft = - b' ./ K2; % solves second equation of vlassov poisson
    dphi_dx_h = phi_fft.*kx(1,:)';
return