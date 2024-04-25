n=2^5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Nmap         =n; % grid of the submaps
params.Nfine        =2^9; % grid dimension for measuring and hist
params.Nsampling    =2^8; % grid for upsampling the velocity after solving laplace
params.Nplotting    =1024; % plotting grid
params.nv           = n;
params.Lv           = 4*pi;
params.detTol       = 1e-2; % incomp threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
params.l            = 1;
params.k            = 0.5;
params.eps          = 5e-2;
params.v0           = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain size
params.Lx           = 2*pi/params.k;
params.L = [params.Lx, params.Lv*2];                                        % domain size
params.dom = [0, 0, params.Lx, 2*params.Lv];                                % domain boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.CFL          = 0.8;
params.T_end        = 80;
params.dt           = 0.001/10;
params.iplot        = 100; % plot every iplot time steps
params.ihist        = 100; % only used if dt_hist not set
params.ilog         = 10;  %only used if dt_log not set
params.dt_hist      = 1;
params.dt_log       = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.bump_transition_width = 0.1*params.Lv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options:
params.case            = 'landau_damping';
params.filter          = 'gauss';