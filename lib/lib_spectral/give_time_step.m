function dt = give_time_step(params, Efield, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function selects dt using cfl criteria 
% and furthermore adapts dt to respect logging and
% backup time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first check cfl criteria
dt = params.dt;
if (params.CFL>0)
        vmax = params.Lv;
        dt = min(params.dx/vmax,params.dv/max(abs(Efield),[],'all'))^(4/3)*params.CFL;
end


% now decide if the time step should be adapted to match the logging and
% hist time step size
dt_log = params.dt_log;
dt_hist = params.dt_hist;

N_log = floor(t/dt_log);
time2_log = (N_log+1)*dt_log;
N_hist = floor(t/dt_hist);
time2_hist = (N_hist+1) * dt_hist;

% if the time step will exceed the next time to log then we make it smaller
if (t+dt > (time2_log) && abs(time2_log-t) > 0)
    dt = time2_log - t;
end

% if the time step will exceed the next time for hist, we make it smaller
if (t+dt > (time2_hist) && abs(time2_hist-t) > 0)
    dt = time2_hist - t;
end

assert(dt>1e-15,"time step to small")
