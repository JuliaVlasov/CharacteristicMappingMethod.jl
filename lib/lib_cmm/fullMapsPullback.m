function [xb, yb, JB,  F0, F1, F2, F3] = fullMapsPullback(dom, gsz, BMList, Form0List, Form1List, Form2List, Form3List, xb, yb)
% general backward map composition and pullabcks for all forms
nf0 = size(Form0List, 2);
nf1 = size(Form1List, 2);
nf2 = size(Form2List, 2);
nf3 = size(Form3List, 2); % volume forms

F0 = cell(1,nf0);
F1 = cell(1,nf1);
F2 = cell(1,nf2);
F3 = cell(1,nf3);

for i0 = 1:nf0
    F0{i0} = 0;
end
for i1 = 1:nf1
    F1{i1} = cat(3, 0, 0);
end
for i2 = 1:nf2
    F2{i2} = cat(3, 0, 0);
end
for i3 = 1:nf3
    F3{i3} = 0;
end


sz = size(xb);
JB = zeros([sz 4]);
JB(:,:,[1 4]) = 1;


nmaps = numel(BMList);

for i =1:nmaps
    
    k = nmaps -i + 1;

    if(any(isnan(xb), "all") || any(isnan(yb), "all"))
        disp("nan");
        keyboard;
        
        ifexit = false;
        if(ifexit)
            break
        end
    end


    for i0 = 1:nf0
        F0{i0} = F0{i0} + HFunInterp(dom, gsz, 'cubic', Form0List{k,i0}, [0 0], xb, yb);
    end
    for i1 = 1:nf1
        src = HFunInterp(dom, gsz, 'cubic', Form1List{k,i1}, [0 0], xb, yb);
        F1{i1} = F1{i1} + MatTimes(JB, src, true);
    end
    for i2 = 1:nf2
        src = HFunInterp(dom, gsz, 'cubic', Form2List{k,i2}, [0 0], xb, yb);
        F2{i2} = F2{i2} + AdjTimes(JB, src);
    end
    for i3 = 1:nf3
        src = HFunInterp(dom, gsz, 'cubic', Form3List{k,i3}, [0 0], xb, yb);
        F3{i3} = F3{i3} + src.*(JB(:,:,1).*JB(:,:,4) - JB(:,:,2).*JB(:,:,3));
    end


    [xb, yb, jac] = HMapInterp(dom, gsz, 'cubic', BMList{k}, xb, yb);
    JB(:,:,[1 3]) = MatTimes(jac, JB(:,:,[1 3]));
    JB(:,:,[2 4]) = MatTimes(jac, JB(:,:,[2 4]));

    
end



return