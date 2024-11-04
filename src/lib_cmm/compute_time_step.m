function dt = compute_time_step(params, f,u1j ,u2j , t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function selects dt using cfl criteria
% and furthermore adapts dt to respect logging and
% backup time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first check cfl criteria
maxu = max(abs(sqrt(u1j(:,1).^2+u2j(:,1).^2)), [], 'all');
maxw = max(abs(f), [], "all");
dt0 = params.dt0;
if (isfield(params, 'dt_constant')) && params.dt_constant
    dt = dt0;
else
    dt = dt0/max([maxu maxw 1]);
    % now decide if the time step should be adapted to match the logging and
    % hist time step size
    dt_log = params.dt_log;
    dt_hist = params.dt_hist;

    N_log = floor(t/dt_log);
    time2_log = (N_log+1)*dt_log;
    if abs(time2_log-t)<1e-14 % to catch roundoff errors in floor
        time2_log = (N_log+2)*dt_log;
    end

    N_hist = floor(t/dt_hist);
    time2_hist = (N_hist+1) * dt_hist;
    if abs(time2_hist-t)<1e-14
        time2_hist = (N_hist+2)*dt_hist;
    end

    % if the time step will exceed the next time to log then we make it smaller
    if (t+dt > (time2_log) && abs(time2_log-t) > 1e-12)
        dt = time2_log - t;
    end

    % if the time step will exceed the next time for hist, we make it smaller
    if (t+dt > (time2_hist) && abs(time2_hist-t) > 1e-12)
        dt = time2_hist - t;
    end
end


% if the time step will exceed the total time, we make it smaller
if (t+dt > params.T_end && abs(params.T_end-t) > 1e-12)
    dt = params.T_end - t;
end

assert(dt>1e-12,"time step to small")
