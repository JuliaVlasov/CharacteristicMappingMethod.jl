function [indorig, indtar] = FourierResampInd1D(norig, ntar)


ifftshiftmap = @(ind, size) mod( ind-1 - floor(size/2), size ) +1;

par = mod(norig, 2);


diff = ntar - norig;
head = abs(floor( (diff+par)/2 ));   % so that if diffy is even then always diffy/2, but if odd, head > tail
tail = abs(ceil( (diff-par)/2 ));    % and so if positive, the add more to head, if negative then tail is bigger in magnitude


switch sign(diff)
    case 1  % padding with zeros
        
        indtar = (head+1:head+norig);
        indorig = 1:norig;
        
        indtar = ifftshiftmap(indtar, ntar);
        indorig = ifftshiftmap(indorig, norig);
        
    case -1
        
        indtar = 1:ntar;
        indorig = head+1:norig-tail;
        
        indtar = ifftshiftmap(indtar, ntar);
        indorig = ifftshiftmap(indorig, norig);
        
    case 0
        indtar = 1:ntar;
        indorig = 1:norig;
end


return