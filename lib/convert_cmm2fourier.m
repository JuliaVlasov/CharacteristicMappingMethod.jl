function [params] = convert_cmm2fourier(params)

params.nx = params.Nfine;
params.nv = params.Nfine;
params.x  = params.Lx*(0:params.nx-1)/(params.nx);
params.v  = 2*params.Lv*((0:params.nv-1)/(params.nv)-0.5); 
params.inicond = params.case;

params.grid.dx = params.x(2) - params.x(1);
params.grid.dv = params.v(2) - params.v(1);
[Xgrid, Vgrid] = meshgrid(params.x, params.v);
params.grid.Vgrid = Vgrid;
params.grid.size = size(Vgrid);
params.dealias      = zeros(params.nx,params.nv);