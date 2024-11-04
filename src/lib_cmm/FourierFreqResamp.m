function fhat_tar = FourierFreqResamp(fhat, origsize, tarsize)


if any((origsize ~= tarsize))
    [indorig1, indtar1] = FourierResampInd1D(origsize(1), tarsize(1));
    [indorig2, indtar2] = FourierResampInd1D(origsize(2), tarsize(2));
    
    fhat_tar = zeros(tarsize);
    
    scale = prod(tarsize./origsize);
    
    fhat_tar(indtar1, indtar2) = fhat(indorig1, indorig2).*scale;
else
   fhat_tar = fhat; 
end

return