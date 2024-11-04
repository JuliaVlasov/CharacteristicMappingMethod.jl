function u1h_fun = defExt_L3(u1h, u1h_m, u1h_mm, tc, tm, tmm)

c = @(t) ((t-tm).*(t-tmm))./((tc-tm).*(tc-tmm));
m = @(t) ((t-tc).*(t-tmm))./((tm-tc).*(tm-tmm));
mm = @(t) ((t-tc).*(t-tm))./((tmm-tc).*(tmm-tm));

u1h_fun = @(t) c(t).*u1h + m(t).*u1h_m + mm(t).*u1h_mm;

end
