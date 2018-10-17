clc;
close all;
clear all;
% This purpose of this code is:
% 1) to compare  ALO and out-of-sample error
% 2) to reproduce Figure 2 in the paper

% IMPORTANT:
% Change the path below, depending on the location of your glmnet_maptlab 
% addpath('/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/code/glmnet_matlab');
rng(0)
rho    =   0.2;
delta  =   0.2;
p      =   10000;% variable number
n      =   p*delta; % number of data
k      =   n*rho; % sparsity
%design_distribution='iid';
design_distribution='spiked';
switch design_distribution
    case 'spiked'
        a_  =  0.3;
        C   =  (1-a_)*spdiags(ones(p,1),0,p,p) + a_*ones(p,1)*ones(1,p);
        F   =  [sqrt(1-a_)*spdiags(ones(p,1),0,p,p)  sqrt(a_)*ones(p,1)];
        %C is equal to F*F'
        %check: sum(sum(C- F*F'))
    case 'iid'
        C   =  spdiags(ones(p,1),0,p,p);
        F   =  spdiags(ones(p,1),0,p,p);
end
X            =    F* randn(size(F,2),n);
X            =    X';
snr          =    1;
o            =    1;
beta         =    zeros(p,1); 
beta(1:k,1)  =    ones(k,1);
X            =    X/sqrt(beta'*C*beta/o^2);
F            =    F/sqrt(beta'*C*beta/o^2);
C            =    C/(beta'*C*beta/o^2);
Xb           =    X * beta;
y            =    Xb + randn(n,1)*o; 
disp('-----------------------------------------------------------------------------');
fprintf('snr=%0.2f| rho=%0.2f| delta=%0.2f \n',snr,rho,delta); 


%fit a preliminary cv to get lambda range
m               =    30;
options         =    glmnetSet;
options.nlambda =    m;
options.thresh  =    1E-9;
options.dfmax   =    ceil(0.65*n);
fit             =    glmnet(X, y,'gaussian',options);
lambda          =    fit.lambda;
lambda_min      =    min(lambda);
lambda_max      =    max(lambda);
lambda          =    logspace(log(lambda_min)/log(10),log(lambda_max)/log(10),m);
options.lambda  =    lambda;
options.dfmax   =    [];
fit             =    glmnet(X, y,'gaussian',options);
lambda          =    fit.lambda;
B               =    fit.beta;
a0              =    fit.a0;
ALO             =    zeros(m,1);
Err_extra       =    zeros(m,1);
df              =    fit.df;

fprintf('lambda min=%g | lambda max=%g \n',min(lambda),max(lambda)); 

for i=1:m
    [Bs, Is]      =   sort(abs(B(:,i)),'descend');
    S             =   Is(1:df(i));
    Xs            =   X(:,S);
    y_hat         =   Xs*B(S,i)+a0(i);% n x 1
    dh            =   diag(Xs*((Xs'*Xs)\Xs'));
    ALO(i)        =   norm((y - y_hat)./(1-dh))^2/n;
    Err_extra(i)  =   o^2 + (beta-B(:,i))'*C*(beta-B(:,i)); 
    r             =   y-y_hat;
    fprintf('iter=%d| df=%d| n=%d| lambda=%d \n',i,df(i),n,lambda(i));
end

disp('-----------------------------------------------------------------------------');

%close all;
%formatSpec = '/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/code/matfiles/extra_vs_in_sample_error_gaussian_%s_p%d';
%save(sprintf(formatSpec,design_distribution,p));

 
close all;
h=figure('units','normalized','outerposition',[0 0 0.7 0.7]);
semilogx(lambda,Err_extra,'b','linewidth',1);hold on;grid on;
semilogx(lambda,ALO,'-.r','linewidth',2);
xlabel('$\lambda$','interpreter','latex','fontsize',24,'fontweight','bold');
leg                 =    legend('$\mbox{Err}_{\mbox{extra},\lambda}$',...
    '$\mbox{ALO}_{\lambda}$', 'location','north');
leg.Interpreter     =    'latex';
leg.FontSize        =    22;
leg.FontName        =    'Courier';
leg.Location        =    'southeast';
title(sprintf(['  p =', num2str(p), ', n =', num2str(n)]),'interpreter','latex','fontsize',20,'fontname','courier')
ylabel('Prediction error','fontsize',28,'fontname','Courier');

axis tight;
yLimits             =     get(gca,'YLim');
set(gca,'fontsize',16);


 set(gcf,'PaperPositionMode','auto');
 formatSpec = '/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/extra_vs_alo_gaussian_lasso_%s_p%d';
 print(h,sprintf(formatSpec,design_distribution,p),'-depsc','-r300');
