clear
tic
M = 2; %choose order of polynomials -> (M+1)x(M+1) matrices
gM = (2*M)^(2/7);
Max = @(v) pi^(-1/2).*exp(-(abs(v)).^2); % 1-dim maxwellian

sig = @(v,w) 1; %choose sigma from application

A = zeros(M+1,M+1);
B = zeros(M+1,M+1);
C = zeros(M+1,M+1);

%------------------calc legendre polynomials---------------------

%Derivatives of first few Legendre polynomials calculated analytically
dLegendreP0 = @(v)0.*v;
dLegendreP1 = @(v)v.^0;
dLegendreP2 = @(v)3.*(v);
dLegendreP3 = @(v)1/2.*(15.*(v).^2-3);
dLegendreP4 = @(v)1/8.*(140.*(v).^3-60.*(v));
dLegendreP5 = @(v)1/8.*(315.*(v).^4-210.*(v).^2+15);
dLegendreP6 = @(v)693/8.*v.^5 - 315/4.*v.^3+105/8.*v;
dLegendreP7 = @(v)3003/16.*v.^6-3465/16.*v.^4+945/16.*v.^2-35/16;
dLegendreP8 = @(v)6435/16.*v.^7-9009/16.*v.^5+3465/16.*v.^3-315/16.*v;
dLegendreP9 = @(v)109395/128.*v.^8-45045/32.*v.^6+45045/64.*v.^4-3465/32.*v.^2+315/128;
dLegendreP10 = @(v)230945/128.*v.^9-109395/32.*v.^7+135135/64.*v.^5-15015/32.*v.^3+3465/128.*v;
%max M=10 possible, otherwise we need to precompute more derivatives...

for i = 0:M
    for j = 0:M
        
        %integrant for A
       fun = @(v) v.*sqrt((2*j+1)/(2*(gM))).*legendreP(j,(v./(gM))).*sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM)));
       A(i+1,j+1) = round(integral(fun,-gM,gM)*10^(10))/10^(10);
       switch j
           %integrant for B
           case  0
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP0((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  1
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP1((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  2
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP2((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  3
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP3((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  4
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP4((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  5
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP5((v./(gM)))-legendreP(j,(v./(gM))).*v);
            case  6
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP6((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  7
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP7((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  8
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP8((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  9
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP9((v./(gM)))-legendreP(j,(v./(gM))).*v);
           case  10
               fun = @(v) sqrt((2*i+1)/(2*(gM))).*legendreP(i,(v./(gM))).*sqrt((2*j+1)/(2*(gM))).*(1/(gM).*dLegendreP10((v./(gM)))-legendreP(j,(v./(gM))).*v);
           otherwise 
               error('nicht genug Ableitungen von Legendre Polynom berechnet');
       end   
       
       B(i+1,j+1) = round(integral(fun,-gM,gM)*10^(10))/10^(10);
       
       %integrants for C
       int1 = @(v)integral(@(w)sig(v,w).*legendreP(j,w./(gM)).*sqrt(Max(w)).*sqrt((2*j+1)/(2*(gM))),-gM,gM);
       int2 = @(v)integral(@(w)Max(w).*sig(v,w),-inf,inf);
       
       fun = @(v) sqrt((2*i+1)/(2*(gM))).*(legendreP(i,v./(gM)).*(sqrt(Max(v)).*int1(v)-sqrt((2*j+1)/(2*(gM))).*legendreP(j,v./(gM)).*int2(v)));  
       C(i+1,j+1) = round(integral(fun,-gM,gM)*10^(10))/10^(10);
    end
end
toc
