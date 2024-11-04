function [x, res] =  NewtonsMethod1D(func, dfunc, fval, x, tol, NIT)


res = (func(x)-fval);
for i = 1:NIT
   if max(abs(res),[], 'all') < tol
       break
   end
   
   dx = res./dfunc(x);
   
   x = x-dx;
   res = (func(x)-fval);
   
end

return