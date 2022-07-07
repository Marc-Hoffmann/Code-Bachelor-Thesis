clear
tic
M =2; %choose order of polynomials -> (M+1)x(M+1) matrices
Max = @(v) pi^(-1/2).*exp(-(abs(v)).^2); % 1-dim maxwellian
sig = @(v,w) 1; % choose sigma from application

A = zeros(M+1,M+1);
B = zeros(M+1,M+1);
C = zeros(M+1,M+1);

%Derivatives of first few Hermite polynomials calculated analytically
dHermiteP0 = @(x)0.*x;
dHermiteP1 = @(x)2.*x.^0;
dHermiteP2 = @(x)8.*x;
dHermiteP3 = @(x)24.*x.^2-12;
dHermiteP4 = @(x)64.*x.^3-96.*x;
dHermiteP5 = @(x)160.*x.^4-480.*x.^2+120;
dHermiteP6 = @(x)384.*x.^5-1920.*x.^3+1440.*x;
dHermiteP7 = @(x)896.*x.^6-6720.*x.^4+10080.*x.^2-1680;
dHermiteP8 = @(x)2048.*x.^7-21504.*x.^5+53760.*x.^3-26880.*x;
dHermiteP9 = @(x)4608.*x.^8-64512.*x.^6+241920.*x.^4-241920.*x.^2+30240;
dHermiteP10 = @(x)10240.*x.^9-184320.*x.^7+967680.*x.^5-1612800.*x.^3+604800.*x;
%max M=10 possible, otherwise we need to precompute more derivatives...



%------------------calc hermite polynomials---------------------
for i = 0:M
    for j = 0:M
        
        %integrant for A
       fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*v.*hermiteH(j,v).*hermiteH(i,v).*pi^(-1/2);
       %without .*Max(v) due to usage of gauss hermite quad
       
       A(i+1,j+1) = round(gausshermi(fun,20)*10^(10))/10^(10); % cut to 10 decimals
         switch j
           
             %integrant for B
           case  0
               fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP0(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case 1
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP1(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  2
               fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP2(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case 3
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP3(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  4
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP4(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  5
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP5(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
            case  6
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP6(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  7
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP7(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  8
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP8(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  9
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP9(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           case  10
                fun = @(v) (2^i*factorial(i))^(-1/2).*(2^j*factorial(j))^(-1/2).*1/2.*hermiteH(i,v).*(dHermiteP10(v)-2.*v.*hermiteH(j,v)).*pi^(-1/2); 
           otherwise 
               error('not enough derivatives available');
         end   
       
       B(i+1,j+1) = round(gausshermi(fun,20)*10^(10))/10^(10);

       

       %integrants for C
       int1 = @(v)integral(@(w)sig(v,w).*hermiteH(j,w).*Max(w).*(2^j*factorial(j))^(-1/2),-inf,inf); % int of sig*phi
       int2 = @(v)integral(@(w)Max(w).*sig(v,w),-inf,inf); % int of sig*Max
       fun = @(v) (2^i*factorial(i))^(-1/2).*hermiteH(i,v).*(int1(v)-(2^j*factorial(j))^(-1/2).*hermiteH(j,v).*int2(v)).*pi^(-1/2); 
 
       C(i+1,j+1) = round(gausshermi(fun,20)*10^(10))/10^(10);
    end
end
toc
