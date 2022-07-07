function [hermf,legf]= polyapprox(f,M,gM)
%f -- function that is to be approximated
%M -- order of polynomial that approximates f
%gM -- interval chopping for Legendre case

Max = @(v) pi^(-1/2).*exp(-(abs(v)).^2); % 1-dim maxwellian
Maxinv = @(v) sqrt(pi).*exp(abs(v).^2); %1/M

herma = zeros(1,M+1);
lega = zeros(1,M+1);

hermf = @(v) 0;
legf = @(v) 0;

for i=0:M 
   % -----------------------Calc coefficients Hermite polynomials-------------------------
   fun = @(v) ((2^i*factorial(i)))^(-1/2).*hermiteH(i,v).*f(v); 
   %f*phi*1/Max: maxwellian cancels with part of phi
   herma(i+1) = integral(fun,-inf,inf); %calc coefficients (scalar product)
   
   % ----------------------calc coefficients Legendre polynomials---------------------------
   fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*f(v).*sqrt(Maxinv(v)); %does not cancel
   lega(i+1) = integral(fun,-gM,gM); %calc coefficients (scalar product)
   
   % ---------------------calc projection of f onto approx space------------------
   hermf = @(v) hermf(v)+herma(i+1).*(((2^i)*factorial(i)))^(-1/2).*hermiteH(i,v).*Max(v);
   legf = @(v) legf(v) + lega(i+1).*sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt(Max(v));
   %one over two coefficients are zero, due to choice of f.
end

end