clear;
clear;
close all;
clc;
% module 1
%cellName    = 'T4C3'; % good grid cell
%cellName    = 'T4C4'; % good grid cell
%cellName    = 'T6C3'; % relatively good
%cellName    = 'T6C6'; % relatively good

% module 2
%cellName   = 'T2C1'; % good grid cell
%cellName   = 'T2C3';  % good grid cell
%cellName   = 'T3C7'; % good grid cell
%cellName   = 'T4C7'; % relatively good

% module 3
%cellName   = 'T9C3';% good grid cell
%cellName   = 'T7C4'; % good grid cell
%cellName   = 'T12C2'; % good grid cell
%cellName   = 'T12C4';% good grid cell
%cellName   = 'T12C3';% good grid cell
%cellName   = 'T7C1';

%Module 1:
%6      8    9   10   11   13   14   15
%T3C5 T3C8 T4C3 T4C4 T4C5 T4C8 T6C3 T6C6

%Module 2: 
%1     2    3    4    5    7    12
%T2C1 T2C3 T2C5 T3C1 T3C3 T3C7 T4C7

%Module 3:
%16    17    18  19   20   21   22    23    24    25    26    27
%T7C1 T7C3 T7C4 T7C5 T7C6 T9C3 T12C2 T12C3 T12C4 T12C5 T12C6 T12C7




% IMPORTANT 
% IMPORTANT:
% Change the path below, depending on the location of your glmnet_maptlab 
% addpath('/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/code/glmnet_matlab');

cellName   = 'T9C3';
formatSpec = ['BEN_' cellName '.mat'];
load(sprintf(formatSpec));
load 'BEN_pos.mat';

% post, posx and posy are the recorded times and corresponding x and y locations
% cellTS contains the spike times
% note that the spike recording sampling is much higher than the xy
% location sampling rate
% in cellTS we observe spikes as close as 1ms
cellTS          =    cellTS(cellTS < max(post)); % we are not interested in spikes beyond the time we recorded the locations
dt              =    0.4; 
time            =    min(post):dt:max(post);
x_              =    interp1(post,posx,time)';
y_              =    interp1(post,posy,time)';
x               =    (x_ - min(x_))/(max(x_) - min(x_));  
y               =    (y_ - min(y_))/(max(y_) - min(y_));  
n               =    50;
%d               =    n*n;
x               =    ceil(x * n);
x               =    min(max(x,1),n);% x is in {1,...,n}
y               =    ceil(y * n);
y               =    min(max(y,0),n-1); % y is in {0,...,n-1}
T               =    length(time);
r               =    sparse(T,1);
spiked.bins     =    ceil((dt+cellTS - min(post))/dt); %times at which spikes were observed
spiked.bins     =    min(max(spiked.bins,1),T);
L               =    length(spiked.bins); %total number of observed spikes

for l=1:L
    r(spiked.bins(l)) = r(spiked.bins(l)) + 1;
end;
mean_log_fac =    mean(log(factorial(r)));
TN           =    zeros(T,1);
I            =    zeros(T,1);
R            =    zeros(n,n);

xx           =    repmat(1:n,n,1);
yy           =    xx';

o_list       =   [n/16 n/8 n/4 n/2];

basis_loc    =    1:n;


p            =    length(basis_loc) * length(basis_loc)*length(o_list);
Z            =    sparse(n*n,p);
k            =    0;
for o=o_list
    for i = basis_loc
        for j = basis_loc
            basis        =   exp(-0.5 *((xx-i).^2 + (yy-j).^2)/o^2);
            basis        =   basis .* (basis > 0.05); 
            k            =   k+1;
            Z(:,k)       =   basis(:);
        end;
    end;
end;


for t=1:T
    I(t)           =   (x(t) -1)*n + n-y(t);
    TN(t)          =   (I(t)-1)*T+t;
    R(n-y(t),x(t)) =   R(n-y(t),x(t))+r(t);
end;
XX       =    sparse(T,n*n);
XX(TN)   =    1; 
XXX      =    XX*Z;
Z        =    Z* spdiags( sqrt(1./diag(XXX'*XXX)),0,p,p);
X        =    XX*Z;
Ip       =    spdiags(ones(p,1),0,p,p);

z               =   ones(p,1);
m               =   10;
ALO             =   zeros(m,1);
LO              =   zeros(m,1);
TRAIN           =   zeros(m,1);
min_alo         =   inf; 
min_lo          =   inf; 
ALO_time        =   0;
LO_time         =   0;

tstart          =   tic;
options         =   glmnetSet;
options.thresh  =   1E-8;
options.lambda  =   logspace(-3,-1,m);
fit             =   glmnet(X, r,'poisson',options);
lambda          =   fit.lambda;
z               =   fit.beta;
a0              =   fit.a0;
FIT_time        =   toc(tstart);

% ALO
% This part is fast. 

for i=1:m
    tstart      = tic;
    S           = find(abs(z(:,i))> 1e-3 );
    Xs          = X(:,S);
    Xz          = Xs*z(S,i)+a0(i);% T x 1
    ez          = exp(Xz);% T x 1
    Hs          = Xs'*spdiags(ez,0,T,T)*Xs;
    invHs       = inv(Hs);
    h           = diag(Xs*(invHs*Xs')) .* ez; 
    ld          = ez-r;
    ldd         = ez;
    Xz_alo      = Xz +  ld./ldd .* h./(1-h);
    ez_alo      = exp(Xz_alo);
    ALO(i)      = -mean(r.*log(ez_alo)-ez_alo)-mean_log_fac;
    ALO_time    = ALO_time+toc(tstart);
    TRAIN(i)    = -mean(r.*log(ez)-ez)-mean_log_fac;
    if ALO(i) < min_alo
        z_alo       =  z(:,i);
        Zx_alo      =  Z*z_alo + a0(i);
        min_alo     =  ALO(i);
        lambda_alo  =  lambda(i);
    end;   
    if i==1
        z_lambda_max    =  z(:,i);
        Zx_lambda_max   =  Z*z_lambda_max + a0(i);
    end;
    if i==m
        z_lambda_min    =  z(:,i);
        Zx_lambda_min   =  Z*z_lambda_min + a0(i);
    end;
    fprintf(' ALO = %g |   TRAIN = %g | lambda=%g     \n',ALO(i),TRAIN(i),lambda(i));
end;




% LO
% This part takes a lot of time
tstart      =  tic;
for t       =  1:T
    r_      = r;
    r_(t)   = [];
    X_      = X;
    X_(t,:) = [];
    fit_    = glmnet(X_, r_,'poisson',options);
    z_      = fit_.beta;
    a0_     = fit_.a0';
    ez_     = exp(X(t,:)*z_+a0_)';
    LO      = LO  + (ez_ - r(t)*log(ez_))/T;
    if ~mod(t,100)
        fprintf('index=%g|  \n',t);
    end;
end;
LO          = LO-mean_log_fac;
LO_time     = toc(tstart);

[min_lo, min_lo_i]  =   min(LO);
z_lo                =   z(:,min_lo_i);
Zx_lo               =   Z*z_lo + a0(min_lo_i);
lambda_lo           =   lambda(min_lo_i);
fprintf('--------------------------------- \n');

for i=1:m
    fprintf('ALO = %g | LO = %g |  TRAIN = %g | lambda=%g     \n',ALO(i),LO(i),TRAIN(i),lambda(i));
end;

fprintf('ALO time(sec) = %g | LO time (sec) = %g |      \n',ALO_time,LO_time);


%formatSpec = ['./../matfiles/lasso_' cellName];
%save(sprintf(formatSpec));






%%
close all;
h=figure('units','normalized','outerposition',[0 0 1 1]);
fprintf('--------------------------------- \n');

subplot 231; 
plot(x_,y_,'linewidth',1);
hold on; 
plot(x_(spiked.bins),y_(spiked.bins),'.r','markersize',10);%grid on;
set(gca, 'XTick', []);
ylabel('y(cm)');
title(sprintf(cellName));
axis tight;


subplot 232; 
imagesc(reshape((exp(Zx_alo)/dt),n,n));  
title('ALO optimized firing rate (Hz)');
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %  set(gca, 'YTick', []);
set(gca, 'YTick', []);
clb=colorbar('manual','position',[0.63 0.584 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);

subplot 233; 
bar(r);
ylabel('action potetionals per bin');
xlabel(sprintf('bins (%g sec)',dt));
grid on;
axis tight;


subplot 234; 
im=imagesc(reshape((exp(Zx_lambda_min)/dt),n,n));  
title('undersmooth firing rate (Hz)');
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
yticks=ceil((linspace(-150,150,7)-min(y_))/(max(y_)-min(y_))*n);
xlabel('x(cm)');
ylabel('y(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %
set(gca, 'YTick', yticks, 'YTickLabel', linspace(150,-150,7)) % 
clb=colorbar('manual','position',[0.35 0.11 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);

subplot 235; 
semilogx(lambda,ALO,'.-','linewidth',2,'markersize',20);
hold on;
%semilogx(lambda,LO,'.-','linewidth',2,'markersize',20);
xlabel('$\lambda$','interpreter','latex','fontsize',15);
title('cross validation risk');
%leg=legend(sprintf('ALO, time = %g sec',ceil(ALO_time)),sprintf('LO, time = %g sec',ceil(LO_time)));%,'extra-sample error');
%leg.Interpreter='latex';
%leg.FontSize=13;
%leg.FontName='Courier';
%leg.Location='north';
%set(gca,'YAxisLocation','right');
grid on; 




subplot 236; 
imagesc(reshape((exp(Zx_lambda_max)/dt),n,n));  
title('oversmooth firing rate (Hz)');

xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %  
set(gca, 'YTick', []);
clb=colorbar('manual','position',[0.91 0.11 0.01 0.34]);




%%
close all;
h=figure('units','normalized','outerposition',[0 0 1 1]);
fprintf('--------------------------------- \n');

subplot 231; 
plot(x_,y_,'linewidth',1);
hold on; 
plot(x_(spiked.bins),y_(spiked.bins),'.r','markersize',10);%grid on;
set(gca, 'XTick', []);
ylabel('y(cm)');
title(sprintf(cellName));
axis tight;


subplot 232; 
imagesc(reshape((exp(Zx_alo)/dt),n,n));  
title('ALO optimized firing rate (Hz)');
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %  set(gca, 'YTick', []);
set(gca, 'YTick', []);
clb=colorbar('manual','position',[0.63 0.584 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);

subplot 233; 
bar(r);
ylabel('action potetionals per bin');
xlabel(sprintf('bins (%g sec)',dt));
grid on;
axis tight;

% imagesc(reshape(exp(z_lo)/dt,n,n));  
% axis normal; colormap jet;
% title('LO optimized firing rate (Hz)');
% xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
% set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %  
% set(gca, 'YTick', []);
% set(gca, 'XTick', []);
% clb=colorbar('manual','position',[0.91 0.584 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);


subplot 234; 
im=imagesc(reshape((exp(Zx_lambda_min)/dt),n,n));  
title('undersmooth firing rate (Hz)');
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
yticks=ceil((linspace(-150,150,7)-min(y_))/(max(y_)-min(y_))*n);
xlabel('x(cm)');
ylabel('y(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %
set(gca, 'YTick', yticks, 'YTickLabel', linspace(150,-150,7)) % 
clb=colorbar('manual','position',[0.35 0.11 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);

subplot 235; 
semilogx(lambda,ALO,'.-','linewidth',2,'markersize',20);
hold on;
semilogx(lambda,LO,'.-','linewidth',2,'markersize',20);
xlabel('$\lambda$','interpreter','latex','fontsize',15);
title('cross validation risk');
leg=legend(sprintf('ALO, time = %g sec',ceil(ALO_time)),sprintf('LO, time = %g sec',ceil(LO_time)));%,'extra-sample error');
leg.Interpreter='latex';
leg.FontSize=13;
leg.FontName='Courier';
leg.Location='north';
set(gca,'YAxisLocation','right');
grid on; 




subplot 236; 
imagesc(reshape((exp(Zx_lambda_max)/dt),n,n));  
title('oversmooth firing rate (Hz)');

xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %  
set(gca, 'YTick', []);
clb=colorbar('manual','position',[0.91 0.11 0.01 0.34]);
% ticklabels=exp(str2double(clb.TickLabels));
% ticklabels=ceil(ticklabels*100)/100;
% clb.TickLabels=num2str(ticklabels);

set(h,'PaperPositionMode','auto');
formatSpec = ['/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/lasso_' cellName];
print(gcf,sprintf(formatSpec),'-depsc','-r300');
%%

% the rat's trajectory superimposed with spikes 
% and the number spikes per bin


close all;
h=figure('units','normalized','outerposition',[0 0 1 0.7]);
fprintf('--------------------------------- \n');


subplot 121; 
plot(x_,y_,'linewidth',1);
hold on; 
plot(x_(spiked.bins),y_(spiked.bins),'.r','markersize',10);%grid on;
set(gca, 'XTick', []);
ylabel('y(cm)');
title(sprintf(cellName));
axis tight;


subplot 122; 
bar(r);
ylabel('action potetionals per bin');
xlabel(sprintf('bins (%g sec)',dt));
grid on;
axis tight;

        

set(h,'PaperPositionMode','auto');
formatSpec = ['/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/rawData' cellName];
print(gcf,sprintf(formatSpec),'-depsc','-r300');

%%
% The 4 truncated Gaussian bumps

k            =    0;
h            =    figure('units','normalized','outerposition',[0 0 1 0.4]);
for o=o_list
    k   = k+1;
    phi = exp(-0.5 *((xx-(n+1)/2).^2 + (yy-(n+1)/2).^2)/o^2);
    phi = phi .* (phi > 0.05); 
    subplot(1,4,k);
    imagesc(phi);
end;

subplot 141;
yticks=ceil((linspace(-150,150,7)-min(y_))/(max(y_)-min(y_))*n);
ylabel('y(cm)');
set(gca, 'YTick', yticks, 'YTickLabel', linspace(150,-150,7));
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7));
%clb=colorbar('manual','position',[0.47 0.584 0.01 0.34]);


subplot 142;
%set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7));
%clb=colorbar('manual','position',[0.91 0.11 0.01 0.34]);

subplot 143;
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'YTickLabel', []);
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)); 
%set(gca, 'YTick', yticks, 'YTickLabel', linspace(150,-150,7)) % 
%clb=colorbar('manual','position',[0.47 0.11 0.01 0.34]);

subplot 144;
xticks=ceil((linspace(-150,150,7)-min(x_))/(max(x_)-min(x_))*n);
xlabel('x(cm)');
set(gca, 'XTick', xticks, 'XTickLabel', linspace(-150,150,7)) %
set(gca, 'YTickLabel', []);
%clb=colorbar('manual','position',[0.91 0.11 0.01 0.34]);

set(h,'PaperPositionMode','auto');
formatSpec = ['/Users/krad/Dropbox/Fast LOOCV/ALO_JRSSB/figures/basis'];
print(gcf,sprintf(formatSpec),'-depsc','-r300');


