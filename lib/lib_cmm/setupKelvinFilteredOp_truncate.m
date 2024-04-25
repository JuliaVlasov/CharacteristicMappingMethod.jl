function defKFOp = setupKelvinFilteredOp_truncate(R)
% the point is to generate a filter op from a single param

defKFOp = @(varargin) genKFOp(R, varargin{:});
% leave in general form, supply k2, but ppossibly kx, ky if kernel is not radial

end

function KFOp = genKFOp(R, varargin)
% genKFOp will be called once, this kernel will be defined once
% all subsequent calls of KFOp will use this kernel without redef

% input as k2, kx, ky

k2 = varargin{1};
ifkeep = ( abs(k2) <= R^2 );

KFOp = @(wh) wh.*ifkeep;

end

