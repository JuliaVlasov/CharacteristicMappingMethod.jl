
function [phi,f] = inicond_vlasov(params)

    
    switch params.inicond
        case 'landau_damping'
     
            f = 1/params.l*(1+params.eps*cos(params.k*params.X))/sqrt(2*pi).*exp(-params.V.^2/2);
            phi = cofitxy(give_potential(params,f));  % note that for vlassov we return the distribution function not the velocity vecto
        case 'two_stream'
            f = 1/params.l*(1+params.eps*cos(params.k*params.X))/(2*sqrt(2*pi)).*(exp(-(params.V-params.v0).^2/2)+exp(-(params.V+params.v0).^2/2));
            phi = cofitxy(give_potential(params,f));
        otherwise
            error('Initial condition unkown...')
    end
end