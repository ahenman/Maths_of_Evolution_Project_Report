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
sigma = 0.5;
sigma1 = 0.3;%0.5*rand;
sigma2 = 0.15;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 1.5;%1+5*rand;
r=rand+0.001;

%Invasion fitness
f = @(x,y) r*(1-alpha1(y,x,sigma)*K(C1,C2,sigma1,sigma2,x0,x) ...
    /K(C1,C2,sigma1,sigma2,x0,y));

%Numerical Simulation

numberofrealisations = 500; % number of sample paths
T = 30; % end time for simulation
xinitial = 1; %initial resident trait value
% matrix to hold trait value at each timestep, for each realisation
xplot = zeros(T+1,numberofrealisations);
timeplot = zeros(T+1,numberofrealisations);
for i = 1:numberofrealisations
    x=xinitial;
    time=0;
    j=1;
    xplot(j,i)=x;
    timeplot(j,i)=0;

    while time<T
        y = x+0.4*rand-0.2; %mutant value is resident +-random number <= 0.2
        if f(x,y) > 0
            x=y; %update trait value
        end
        time=time+1;%update timestep
        %update matrices
        j=j+1;
        xplot(j,i)=x;
        timeplot(j,i)=time;
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
xlabel('time');
ylabel('Resident trait value');
axis([0 T 0 (xinitial+1)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;