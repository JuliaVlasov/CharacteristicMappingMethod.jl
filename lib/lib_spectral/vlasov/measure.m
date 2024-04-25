function [Mass,Momentum,Epot,Ekin,Etot,spectr_fft,rk] = measure(f,Efield,dx,dv,Vgrid)
       Mass = sum(f,"all")*dx*dv;
       Momentum = sum(f.*Vgrid,"all")*dx*dv;
       Epot = 0.5*sum(Efield.^2)*dx;
       Ekin = 0.5*sum(f.*(Vgrid.^2),"all")*dx*dv;
       Etot = Epot + Ekin;
       %Epot(it)
       %Ekin(it)
       [rk,spectr_fft] = spectrum(f);
end