% --------------------initilizing general stuff------------------------- 
tic
Max = @(v) pi^(-1/2).*exp(-(abs(v)).^2); % 1-dim maxwellian
Maxinv = @(v) sqrt(pi).*exp(abs(v).^2);

f = @(v) 1./(v.^2+1).*Max(v); %regularity k = 1, decaying n = 7 


maxM=10; %choose max order of polynomials

normH = zeros(1,maxM);
normLin = zeros(1,maxM);
normLout = zeros(1,maxM);
normL = zeros(1,maxM);

qH = zeros(1,maxM-2);
qL = zeros(1,maxM-2);

n= 7; % measure for decay behaviour
k= 1; % measure for regularity
for M=1:maxM
    
    gM = (2*M)^((2*k)/(n)); % determine interval length subject to decay behaviour and regularity of f
    
    [hermf,legf]=polyapprox(f,M,gM); % get approximation functions
    
   % ----------------------Error calculation----------------------------   
   
   errH = @(v) (f(v)-hermf(v)).^2.*Maxinv(v); % remember: Maxinv = 1/M
   errLin = @(v) (f(v)-legf(v)).^2.*Maxinv(v); 
   errLout = @(v) (f(v)).^2.*Maxinv(v);

   normH(M) = integral(errH,-20,20); % [-20,20] as otherwise errH explodes
   normLin(M) = integral(errLin,-gM,gM); % calc error inside I
   normLout(M) = integral(errLout,-20,-gM)+integral(errLout,gM,20); % and outside I
   normL(M) = normLin(M)+normLout(M); %total error
   
end


%---------------------------------plots-------------------------------------
xi =linspace(-gM,gM,100);

% figure
% hold on
% plot(xi,hermf(xi),'r'); %plotting approximations and f
% plot(xi,legf(xi),'b');
% plot(xi,f(xi),'k');

%plot convergence behaviour
yi= 1:maxM;
figure 
loglog(yi,normH(yi),'r',yi,normL(yi),'b',yi,((2.*yi).^-1)./exp(3/2),'k');
legend('Hermite','Legendre','(2M)^{-1}');

toc

        