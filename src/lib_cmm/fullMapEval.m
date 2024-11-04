function [x, y, Jac] = fullMapEval(dom, gsz, BMList, x, y)

sz = size(x);
Jac = zeros([sz 4]);
Jac(:,:,[1 4]) = 1;
% J = Jac;

nmaps = numel(BMList);

for i =1:nmaps
    
    k = nmaps -i + 1;

    if(any(isnan(x.*y), "all"))
        disp("nan");
        keyboard;
        
        ifexit = false;
        if(ifexit)
            break
        end
    end

    [x, y, jac] = HMapInterp(dom, gsz, 'cubic', BMList{k}, x, y);
    Jac(:,:,[1 3]) = MatTimes(jac, Jac(:,:,[1 3]));
    Jac(:,:,[2 4]) = MatTimes(jac, Jac(:,:,[2 4]));

    
end


return