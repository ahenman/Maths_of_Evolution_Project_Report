%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

%%Default graphing layout from Dennis
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 25); % 25 for 1X3, 20 for 1X2 figures
%set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','none');
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesXGrid','off');
set(0,'defaultAxesYGrid','off');
set(0,'defaultAxesTickDir','out');
set(0,'defaultAxesLineWidth',1.5);

%(Mexican hat) Wavelet Competition function
function alpha1=alpha1(x,y,sigma)
    alpha1 = (1.-((x-y)./sigma).^2).*exp(-(x-y).^2./(2*sigma^2));
end
%z=x-y

%Double Gaussian Carrying capacity K(x)
function Gauss=G(sigma,x0,x)
    Gauss = exp(-(x-x0).^2./(2*sigma^2));
end

function K=K(C1,C2,sigma1,sigma2,x0,x)
    K = C1*G(sigma1,x0,x)-C2*G(sigma2,x0,x);
end

%Parameters
C1 = 2.5;%5*rand;
C2 = 2;%C1*rand;
sigma = 1;%needs to be less than 2.950390509 for branching
sigma1 = 3;%0.5*rand;
sigma2 = 1.5;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 10;%1+5*rand;
r=rand+0.001;
m=0.1; %max size of mutation

%Equilibria
N1 = @(x) K(C1,C2,sigma1,sigma2,x0,x); %1 resident equilibrium
Nfun = @(x,xb) (K(C1,C2,sigma1,sigma2,x0,x)-alpha1(x,xb,sigma).*K(C1,C2,sigma1,sigma2,x0,xb))./(1-alpha1(x,xb,sigma).^2);
Nbfun = @(x,xb) (K(C1,C2,sigma1,sigma2,x0,xb)-alpha1(x,xb,sigma).*K(C1,C2,sigma1,sigma2,x0,x))./(1-alpha1(x,xb,sigma).^2);

%Invasion fitness
f = @(x,y) r*(1-alpha1(y,x,sigma)*K(C1,C2,sigma1,sigma2,x0,x) ...
    /K(C1,C2,sigma1,sigma2,x0,y));

f2 = @(x,xb,y) r*(1-(alpha1(y,x,sigma)*(K(C1,C2,sigma1,sigma2,x0,x)-alpha1(x,xb,sigma)*K(C1,C2,sigma1,sigma2,x0,xb)) ...
    +alpha1(y,xb,sigma)*(K(C1,C2,sigma1,sigma2,x0,xb)-alpha1(x,xb,sigma)*K(C1,C2,sigma1,sigma2,x0,x))) ...
    /(K(C1,C2,sigma1,sigma2,x0,y)*(1-alpha1(x,xb,sigma)^2)));

%Evolutionary singular strategies
x1=x0;
x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

%% Numerical Simulation

numberofrealisations = 7; % number of sample paths
T = 700; % end time for simulation
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
    N=N1(x);
    Nb=0;        %initially there is no branch population
    time=0;
    j=1;
    xplot(j,i)=x;       %To keep track of trait values for simulation i at "time" j
    xbplot(j,i)=xb;     %To keep track of branched trait values for simulation i at "time" j
    Nplot(j,i)=N;       %To keep track of population density for simulation i at "time" j
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
                    N=Nfun(x,xb);  %branch point affects both equilibrium densitities
                    Nb=Nbfun(x,xb);
                    branch=true;
                else
                    x=y; %update trait value if y has positive invasion fitness
                    xb=x;
                    N=N1(x); %update new population         
                end
            end
        %This runs once we have got sufficiently close to "branching" point
        else
            if rand<=N/(N+Nb) %pick which subpopulation is mutating weighted by its frequency in the population
                y = x+2*m*rand-m;
                if f2(x,xb,y)>0
                    x=y;
                    N=Nfun(x,xb);
                    Nb=Nbfun(x,xb);
                end
            else
                yb = xb+2*m*rand-m;
                if f2(x,xb,yb)>0
                    xb=yb;
                    N=Nfun(x,xb);
                    Nb=Nbfun(x,xb);
                end
            end
        end
        time=time+1;%update timestep
        %update matrices
        j=j+1;
        xplot(j,i)=x;
        xbplot(j,i)=xb;
        Nplot(j,i)=N;
        Nbplot(j,i)=Nb;
        timeplot(j,i)=time;
    end
end

%Plotting
figure;
C = orderedcolors("gem");
%{
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
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),xplot(:,i),'Color',C(mod(i,7)+1,:));
    k=stairs(timeplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:));
    set(h,'Linewidth',2);
    set(k,'Linewidth',2);
end
mean(mean(abs(xplot-xbplot)))
hold on;
hx1=plot([0,T],x1*ones(1,2),'k-',LineWidth=2);
hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
legend(Location="southeast")
legend([hx2 hx1 hx3], {sprintf('x_0+x_1 = %.2f', x2), sprintf('x_0 = %.2f', x1), sprintf('x_0-x_1 = %.2f', x3)});
%legend([hx1 hx3], {sprintf('x_0   = Repeller'), sprintf('x_0-x_1 = Branch point')});
xlabel('Time','FontSize',25);
ylabel('Trait value','FontSize',25);
axis([0 T 0 max(x2+1,xinitial)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
axis square;

figure;
Cb = orderedcolors("glow");
hold on;
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

%% Branching video
%{
 vidfile = VideoWriter('testmovie.mp4(1)','MPEG-4');
 vidfile.FrameRate=40;
 open(vidfile);
 f=figure;
 for j = 1:T-3
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
    writeVideo(vidfile, getframe(gcf));
 end
close(vidfile)
%}