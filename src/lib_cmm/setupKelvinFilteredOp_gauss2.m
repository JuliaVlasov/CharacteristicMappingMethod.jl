function defKFOp = setupKelvinFilteredOp_gauss2(R)
% the point is to generate a filter op from a single param

% use biLaplacian to smooth, conv with gaussian might be too smooth
% R should be the parameter for kernel = 1e-7 or smt of the type, some equivalent to truncation
% the power on the Laplacian must be high enough so that the transition from O(1) to O(1e-15) be
% fast enough

defKFOp = @(varargin) genKFOp(R, varargin{:});
% leave in general form, supply k2, but ppossibly kx, ky if kernel is not radial

end

function KFOp = genKFOp(R, varargin)
% genKFOp will be called once, this kernel will be defined once
% all subsequent calls of KFOp will use this kernel without redef

% chosen so that at k = R, filter = 2^-26
sig = R^4/(26*log(2));

k2 = varargin{1};

KFOp = @(wh) wh.*exp(-(abs(k2).^2)./sig);

end

