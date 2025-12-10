%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

%%Default graphing layout from Dennis
set(0,'defaultTextFontName', 'Arial');
set(0,'defaultaxesfontsize', 20); % 25 for 1X3, 20 for 1X2 figures
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
sigma = 1;
sigma1 = 3;%0.5*rand;
sigma2 = 1.5;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 10;%1+5*rand;
r=rand+0.001;
m=0.1; %max size of mutation

%Invasion fitness
f = @(x,y) r*(1-alpha1(y,x,sigma)*K(C1,C2,sigma1,sigma2,x0,x) ...
    /K(C1,C2,sigma1,sigma2,x0,y));

%Evolutionary singular strategies
x1=x0;
x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

%Numerical Simulation

numberofrealisations = 5; % number of sample paths
T = 500; % end time for simulation
xinitial = 1; %initial resident trait value
% matrix to hold trait value at each timestep, for each realisation
xplot = zeros(T+1,numberofrealisations);
timeplot = zeros(T+1,numberofrealisations);
for i = 1:numberofrealisations %Run numberofrealisations simulations
    x=xinitial;
    time=0;
    j=1;
    xplot(j,i)=x;       %To keep track of trait values for simulation i at "time" j
    timeplot(j,i)=0;    %To keep track of time that trait x is obtained for simulation i
    branch = false;     %means the trait hasn't branched yet

    while time<T        %goes until specified time T
        %only runs for the trait value sufficiently far from x*
        if abs(x-x2)>m && abs(x-x3)>m && branch==false
            y = x+2*m*rand-m; %mutant value is resident +-random number <= m
            if f(x,y) > 0
                x=y; %update trait value if y has positive invasion fitness
            end
            time=time+1;%update timestep
            %update matrices
            j=j+1;
            xplot(j,i)=x;
            timeplot(j,i)=time;
        %This runs once we have got sufficiently close to "branching" point
        else
            branch=true;
        end
    end
end

%Plotting
figure;
hold on;
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),xplot(:,i));
    set(h,'Linewidth',2);
end
hold on;
hx1=plot([0,T],x1*ones(1,2),'k-',LineWidth=2);
hx2=plot([0,T],x2*ones(1,2),'r-.',LineWidth=2);
hx3=plot([0,T],x3*ones(1,2),'b--',LineWidth=2);
legend(Location="southeast")
legend([hx1 hx2 hx3], {sprintf('x^* = %.2f', x1), sprintf('x^* = %.2f', x2), sprintf('x^* = %.2f', x3)});
xlabel('time');
ylabel('Resident trait value');
axis([0 T 0 (x2+1)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;

%imagesc
%plot of sigma1 by 2



