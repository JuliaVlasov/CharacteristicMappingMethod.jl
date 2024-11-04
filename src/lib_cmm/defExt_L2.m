function u1h_fun = defExt_L2(u1h, u1h_m, tc, tm)

c = @(t) (t-tm)./(tc-tm);
m = @(t) (t-tc)./(tm-tc);

u1h_fun = @(t) c(t).*u1h + m(t).*u1h_m;

end
