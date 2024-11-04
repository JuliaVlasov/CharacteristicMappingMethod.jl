function fjet = RK4_Lie1form(dom, sgsz, fj_fun, mgsz, M, t, dt)

f1 = HPullback1formjet(dom, sgsz, fj_fun(t), mgsz, M(:,:,:,1));
f2 = HPullback1formjet(dom, sgsz, fj_fun(t+0.5*dt), mgsz, M(:,:,:,2));
f3 = HPullback1formjet(dom, sgsz, fj_fun(t+0.5*dt), mgsz, M(:,:,:,3));
f4 = HPullback1formjet(dom, sgsz, fj_fun(t+dt), mgsz, M(:,:,:,4));

fjet = dt*(f1 + 2*f2 + 2*f3 + f4)/6;

return