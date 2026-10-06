%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

%%Default graphing layout from Denis
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 25); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%% functions

%(Mexican hat) Wavelet Competition function
function alpha=alpha(x,y,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).* ...
(1.-((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2));
end

%dalpha(x-y)/dx, dalpha(y,x,...) = dalpha(y-x)/dy = -dalpha(y-x)/dx
function dalpha=dalpha(x,y,sigma,beta)
    dalpha = exp(sigma.^2.*beta.^2/2).* ...
(-3+((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2)).*(x-y+sigma.^2.*beta)./sigma.^2;
end
%d^2alpha(x-y)/dx^2, ddalpha(y,x,...) = d^2alpha(y-x)/dy^2
function ddalpha=ddalpha(x,y,sigma,beta)
A = (x-y+sigma^2*beta);
    ddalpha = exp(sigma.^2.*beta.^2/2-A.^2./(2*sigma^2)).* ...
((-3+(A./sigma).^2)+ ...
(A./sigma^2).^2.*(5-(A./sigma).^2));
end

%Gaussian centred on x0
function Gauss=G(sigma,x0,x)
    Gauss = exp(-(x-x0).^2./(2*sigma^2));
end
%derivative of gaussian
function dGauss=dG(sigma,x0,x)
    dGauss = -(x-x0)./(sigma^2).*exp(-(x-x0).^2./(2*sigma^2));
end
%second derivative of gaussian
function ddGauss=ddG(sigma,x0,x)
    ddGauss = 1/(sigma^2)*(-1+(x-x0).^2/(sigma^2)).*exp(-(x-x0).^2./(2*sigma^2));
end
%K is sum of 2 gaussians
function K=K(C1,C2,sigma1,sigma2,x0,x)
    K = C1*G(sigma1,x0,x)-C2*G(sigma2,x0,x);
end
%derivative of K (when alpha peaks at y=x, the solutions to this are x*)
function dK=dK(C1,C2,sigma1,sigma2,x0,x)
    dK = C1*dG(sigma1,x0,x)-C2*dG(sigma2,x0,x);
end

%second derivative of K
function ddK=ddK(C1,C2,sigma1,sigma2,x0,x)
    ddK = C1*ddG(sigma1,x0,x)-C2*ddG(sigma2,x0,x);
end

%% Parameters
C1 = 1;%2.5;%5*rand;
C2 = 0;%2;%C1*rand;
sigma = 0.5;%1;%needs to be less than 2.950390509 for branching
sigma1 = 0.5;%0.5*rand;
sigma2 = 0.35;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 1.5;%10;%1+5*rand;
beta=0;
r=rand+0.001;
m=0.01; %max size of mutation
eps = 1e-6; %fraction of population that mutates
xmin=0;
xmax=3;

%% Equilibria
N1 = @(x) K(C1,C2,sigma1,sigma2,x0,x)/alpha(0,0,sigma,beta); %1 resident equilibrium
%equilibrium for 2 residents 
%The equilibrium pop density of resident with trait x1, swap x1 and x2 to
%get pop dens of resident with trait x2
function N = N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)
    N = (K(C1,C2,sigma1,sigma2,x0,x1).*alpha(0,0,sigma,beta)- ...
        K(C1,C2,sigma1,sigma2,x0,x2).*alpha(x1,x2,sigma,beta))./( ...
        alpha(0,0,sigma,beta).^2-alpha(x1,x2,sigma,beta).*alpha(x2,x1,sigma,beta));
end
Nfun = @(x,xb) N(x,C1,C2,sigma1,sigma2,x0,sigma,beta,xb);
%Invasion fitness
f = @(x,y) r*(1-alpha(y,x,sigma,beta).*K(C1,C2,sigma1,sigma2,x0,x) ...
    ./(alpha(0,0,sigma,beta).*K(C1,C2,sigma1,sigma2,x0,y)));

%fitness for 2 residents
function fit12 =fit12(r,sigma,beta,C1,C2,sigma1,sigma2,x0,y,x1,x2)
    fit12 = r*(1-(alpha(y,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(y,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,y));
end
f2 = @(x1,x2,y) fit12(r,sigma,beta,C1,C2,sigma1,sigma2,x0,y,x1,x2);
%selection gradient for branch with resident x1
function secgrad1 = secgrad1(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad1 = r*(-(dalpha(0,0,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        dalpha(x1,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x1)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x1).*(alpha(0,0,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(x1,x2,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x1).^2);
end

%selection gradient for branch with resident x2
function secgrad2 = secgrad2(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2)
    secgrad2 = r*(-(dalpha(x2,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        dalpha(0,0,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x2)+ ...
        dK(C1,C2,sigma1,sigma2,x0,x2).*(alpha(x2,x1,sigma,beta).*N(x1,C1,C2,sigma1,sigma2,x0,sigma,beta,x2)+ ...
        alpha(0,0,sigma,beta).*N(x2,C1,C2,sigma1,sigma2,x0,sigma,beta,x1))./K(C1,C2,sigma1,sigma2,x0,x2).^2);
end

sx1 = @(x1,x2) secgrad1(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2);
sx2 = @(x1,x2) secgrad2(r,sigma,beta,C1,C2,sigma1,sigma2,x0,x1,x2);
%Evolutionary singular strategies
x1=x0;
%x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
%x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

%% Numerical Simulation

numberofrealisations = 1; % number of sample paths
T = 2000; % end time for simulation
xinitial = 1; %initial resident trait value
% matrix to hold trait value at each timestep, for each realisation
xplot = zeros(T+1,numberofrealisations);
Nplot = zeros(T+1,numberofrealisations);
Nbplot = zeros(T+1,numberofrealisations);
xbplot = zeros(T+1,numberofrealisations);
timeplot = zeros(T+1,numberofrealisations);
for i = 1:numberofrealisations %Run numberofrealisations simulations
    x=xinitial;
    xb=xinitial; %resident of the branch
    Nr=N1(x);
    Nb=Nr;        %initially there is no branch population
    time=0;
    j=1;
    xplot(j,i)=x;       %To keep track of trait values for simulation i at "time" j
    xbplot(j,i)=xb;     %To keep track of branched trait values for simulation i at "time" j
    Nplot(j,i)=Nr;       %To keep track of population density for simulation i at "time" j
    Nbplot(j,i)=Nb;
    timeplot(j,i)=0;    %To keep track of time that trait x is obtained for simulation i
    branch = false;     %means the trait hasn't branched yet

    while time<T        %goes until specified time T
        %only runs for the trait value sufficiently far from x*
        if branch==false
            y = x+2*m*rand-m; %mutant value is resident +-random number <= m
            if f(x,y) > 0
                if f(y,x)>0
                    xb=y; %update branch value if mutally invasive
                    Nr=Nfun(x,xb);  %branch point affects both equilibrium densitities
                    Nb=Nfun(xb,x);
                    branch=true;
                else
                    x=y; %update trait value if y has positive invasion fitness
                    xb=x;
                    Nr=N1(x); %update new population
                    Nb=Nr;
                end
            end
        %This runs once we have got sufficiently close to "branching" point
        else
            if rand<=Nr/(Nr+Nb) %pick which subpopulation is mutating weighted by its frequency in the population
                y = x+2*m*rand-m;
                if f2(x,xb,y)>0
                    x=y;
                    Nr=Nfun(x,xb);
                    Nb=Nfun(xb,x);
                    %need to check mutations don't push traits outside of
                    %mutualinvasibility zone
                    if Nr<=0 %branch with trait x dies, only xb remains
                        x=xb;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    elseif Nb<=0
                        xb=x;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    end
                end
            else
                yb = xb+2*m*rand-m;
                if f2(x,xb,yb)>0
                    xb=yb;
                    Nr=Nfun(x,xb);
                    Nb=Nfun(xb,x);
                    if Nr<=0 %branch with trait x dies, only xb remains
                        x=xb;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    elseif Nb<=0
                        xb=x;
                        Nr=N1(x);
                        Nb=Nr;
                        branch = false;
                    end
                end
            end
        end
        time=time+1;%update timestep
        %update matrices
        j=j+1;
        xplot(j,i)=x;
        xbplot(j,i)=xb;
        Nplot(j,i)=Nr;
        Nbplot(j,i)=Nb;
        timeplot(j,i)=time;
    end
end

%Plotting
figure;
C = orderedcolors("gem");
%{
%green colour theme for poster
C=[0.757 0.820 0.122
    0.596 0.788 0.169
    0.431 0.753 0.027
    0.318 0.694 0.020
    0.204 0.635 0.012
    0.224 0.580 0.027
    0.239 0.525 0.043
    0.122 0.455 0.051
    0.000 0.380 0.055];
%}
hold on;
%plot trait value by time for each realisation
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),xplot(:,i),'Color',C(mod(i,7)+1,:));
    k=stairs(timeplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
    set(k,'Linewidth',2);
end
mean(mean(abs(xplot(end-200:end,:)-xbplot(end-200:end,:)))) %finds mean separation between branches
hold on;
%plots singular strategies
hx1=yline(x1,'--',LineWidth=2);
%hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
%hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
legend(Location="southeast")
%legend([hx2 hx1 hx3], {sprintf('x_0+x_1 = %.2f', x2), sprintf('x_0 = %.2f', x1), sprintf('x_0-x_1 = %.2f', x3)});
%legend([hx1 hx3], {sprintf('x_0   = Repeller'), sprintf('x_0-x_1 = Branch point')});
legend(hx1,sprintf('x* = %.2f', x1))
xlabel('Time','FontSize',25);
ylabel('Trait value','FontSize',25);
axis([0 T 0 max(x0+1,xinitial)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%pop dens
figure;
Cb = orderedcolors("glow");
hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),Nplot(:,i),'Color',C(mod(i,7)+1,:));
    k=stairs(timeplot(:,i),Nbplot(:,i),'Color',Cb(mod(i,7)+1,:));
    set(h,'Linewidth',2);
    set(k,'Linewidth',2);
end
hold on;
xlabel('Time','FontSize',25);
ylabel('Population Density','FontSize',25);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%total pop dens
figure;
hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),Nplot(:,i)+Nbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
end
hold on;
xlabel('Time','FontSize',25);
ylabel('Population Density','FontSize',25);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

%trait evolution on a TEP

z=linspace(xmin,xmax,2000);
[X1,X2] = meshgrid(z, z);

F1 = f(X1,X2); %x1 resident, x2 mutant
F2 = f(X2,X1); %x2 resident, x1 mutant
S1 = sx1(X1,X2);
S2 = sx2(X1,X2);
%selection gradient
sx = @(x) r*(-dalpha(0,0,sigma,beta0)./alpha(0,0,sigma,beta0)+dK(C1,C2,sigma1,sigma2,x0,x) ...
    ./K(C1,C2,sigma1,sigma2,x0,x));
xstar=roots(chebfun(sx,[1e-6 5])); %singular strategies

F = zeros(size(X1));
F(F1>0)=1;   %fx(y)>0
F(F2>0)=2;   %fy(x)>0
F(F1>0&F2>0)=3;   %mutual invasibility

%don't care about direction outside of coexistence region
%S1(f(X1,X2)<0|f(X2,X1)<0)=0;
%S2(f(X1,X2)<0|f(X2,X1)<0)=0;

figure;
hold on
% TEP
contourf(X1,X2,F,'LineStyle','none');
contour(X1,X2,F,[0.5 1.5 2.5],'k','LineWidth',2); %black line for f=0
colormap([1 1 1; 0.8 0.9 0.9; 1 0.9 0.7; 0.6 0.9 0.3]); %blue for fx(y)>0, yellow for fy(x)>0, green for both
xlabel('$\textsf{Resident trait value}$ ($x_1$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Resident trait value}$ ($x_2$)','Interpreter','latex','FontSize',20);
axis equal
%selection gradient nullclines
contour(X1,X2,S1,[0 0],'r:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_1)=0');
contour(X1,X2,S2,[0 0],'b:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_2)=0');
%singular strategies
xline(xstar(1),'LineStyle','--','LineWidth',1.5,'Color',C(3,:),'DisplayName',['B1=',num2str(xstar(1),3)])
%xline(xstar(2),'LineStyle','--','LineWidth',1.5,'Color',C(2,:),'DisplayName',['D=',num2str(xstar(2),3)])
xline(xstar(end),'LineStyle','--','LineWidth',1.5,'Color',C(1,:),'DisplayName',['B2=',num2str(xstar(end),3)])

hold on;
%plots population density of each branch over time
for i=1:numberofrealisations
    h=stairs(xplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
end
hold on;

legend("off");
axis equal
box on

%% Branching video
%{
 vidfile = VideoWriter('testmovie.mp4(1)','MPEG-4');
 vidfile.FrameRate=40;
 open(vidfile);
 f=figure;
 for j = 1:T-3 %because I plot 3 steps at a time
    f.Name = ['Simulation time: t = ', num2str((j-1))];
    %Plotting
    hold on;
    for i=1:numberofrealisations
        h=stairs(timeplot([j:j+3],i),xplot([j:j+3],i),'Color',C(mod(i,7)+1,:));
        k=stairs(timeplot([j:j+3],i),xbplot([j:j+3],i),'Color',C(mod(i,7)+1,:));
        set(h,'Linewidth',2);
        set(k,'Linewidth',2);
    end
    hold on;
    hx1=plot([0,T],x1*ones(1,2),'k-',LineWidth=2);
    hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
    hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
    legend(Location="southeast")
    %legend([hx2 hx1 hx3], {sprintf('x_0+x_1 = %.2f', x2), sprintf('x_0 = %.2f', x1), sprintf('x_0-x_1 = %.2f', x3)});
    legend([hx1 hx3], {sprintf('x_0   = repeller'), sprintf('x_0-x_1 = branch point')});
    xlabel('time','FontSize',25);
    ylabel('Trait value','FontSize',25);
    axis([0 T 0 x0+1]);
    set(gca,'linewidth',1.5);
    set(gca,'FontSize',20);
    grid on;
    axis square;
    writeVideo(vidfile, getframe(gcf));%takes screenshot of the 3 new steps
    each time
 end
close(vidfile)
%}