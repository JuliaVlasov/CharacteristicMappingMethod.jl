function defKFOp = setupKelvinFilteredOp_gauss_2p(R, p)
% the point is to generate a filter op from a single param

% use biLaplacian to smooth, conv with gaussian might be too smooth
% R should be the parameter for kernel = 1e-7 or smt of the type, some equivalent to truncation
% the power on the Laplacian must be high enough so that the transition from O(1) to O(1e-15) be
% fast enough

defKFOp = @(varargin) genKFOp(R, p, varargin{:});
% leave in general form, supply k2, but ppossibly kx, ky if kernel is not radial

end

function KFOp = genKFOp(R, p, varargin)
% genKFOp will be called once, this kernel will be defined once
% all subsequent calls of KFOp will use this kernel without redef

% chosen so that at k = R, filter = 2^-26
sig = R^(2*p)/(26*log(2));

k2 = varargin{1};

KFOp = @(wh) wh.*exp(-(abs(k2).^p)./sig);

end



% note: for given p, k=R is 2e-26 approx 1e-8
% k which gives 2e-3 approx 1e-1 is such that R^2p - k^2p = 23log(2)*sig = 3R^2p
% i.e. k = (R^2p - 23log(2)*sig)^(1/2p) = (R^2p (3/26) )^(1/2p) = R*((3/26)^(1/2p))

% rule of thumb, want k0 to be within 20% of R
% approx 3/26 = 1/9, k0 = log(1/3)/log(0.8)
% so about 5, i.e. we take exp(k^10)