function [Mass,Momentum,Epot,Ekin,Etot,L2norm,spectr_fft,rk,Emodes] = measure(f,Efield,grid)
       Sigma = 1;%grid.Sigma;
       Vgrid = grid.Vgrid;
       if (min(Vgrid,[],"all") > -1e-14)
            Vgrid = Vgrid - grid.L(2)/2;
       end
       Mass = sum(f.*Sigma,"all")*grid.dx*grid.dv;
       Momentum = sum(f.*Vgrid.*Sigma,"all")*grid.dx*grid.dv;
       Epot = 0.5*sum(Efield.^2)*grid.dx;
       Ekin = 0.5*sum(f.*(Vgrid.^2).*Sigma,"all")*grid.dx*grid.dv;
       Etot = Epot + Ekin;
       Emodes = fourier_modes(Efield, [1,2,3,4], grid);
       L2norm = sum(abs(f).^2.*Sigma,"all")*grid.dx*grid.dv;
       [rk,spectr_fft] = spectrum(f);
end


function [Emode_abs] = fourier_modes(Efield, k_list, grid)
% this function implements (63-65) of https://arxiv.org/pdf/1009.3046.pdf
% for odd fourier numbers (k=0.5 etc)
    ik = 1;
    Ek = fft(Efield);
    Emode_abs = abs(Ek(2:5))/grid.size(1);
%     for k = k_list
%         Esin_k = sum(Efield.*sin(2*pi*k*grid.x/grid.L(1)))*grid.dx;
%         Ecos_k = sum(Efield.*cos(2*pi*k*grid.x/grid.L(1)))*grid.dx;
%         Emode_abs(ik) = sqrt(abs(Esin_k)^2+abs(Ecos_k)^2)/grid.L(1);
%         ik = ik + 1;
%     end

end