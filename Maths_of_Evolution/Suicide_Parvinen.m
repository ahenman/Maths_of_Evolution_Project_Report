%%Alice Henman%%
%%Chapter 3 Evolutionary Suicide%%
%%Classifying singular strategies%%

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

%% Resource consumer model N by R
tspan=[0 500];
function dNdt = dNdt(N,R,d,f,s)
    dNdt = N*(f*s*R-d);
end

function dRdt = dRdt(N,R,a,b,c,s)
    dRdt = R*(a*R/(1+R)-b-c*R-s*N);
end

function ddt = odefcn1(n,a,b,c,d,f,s)
N = n(1);
R = n(2);
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = dNdt(N,R,d,f,s); % resident population
    ddt(2) = dRdt(N,R,a,b,c,s); % mutant population density
end
a=25;
b=4;
c=4;
d=7;
f=1;
s0=14/3;
R1=(a-b-c-sqrt((b+c-a)^2-4*b*c))/(2*c);
R2=(a-b-c+sqrt((b+c-a)^2-4*b*c))/(2*c);
Rbar=d/(f*s0);
Nbar=(a*Rbar-(b+c*Rbar)*(1+Rbar))/(s0*(1+Rbar));

[t,n1] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 1]);
[t,n2] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 2]);
[t,n3] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 3]);
[t,n4] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [1 5]);
[t,n5] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [2 5]);
[t,n6] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [3 5]);
[t,n7] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [0.1 R2]);

figure;
plot(abs(n1(:,1)),abs(n1(:,2)),LineWidth=2.5);
hold on;
plot(abs(n2(:,1)),abs(n2(:,2)),LineWidth=2.5);
plot(abs(n3(:,1)),abs(n3(:,2)),LineWidth=2.5);
plot(abs(n4(:,1)),abs(n4(:,2)),LineWidth=2.5);
plot(abs(n5(:,1)),abs(n5(:,2)),LineWidth=2.5);
plot(abs(n6(:,1)),abs(n6(:,2)),LineWidth=2.5);
plot(abs(n7(:,1)),abs(n7(:,2)),LineWidth=2.5);
plot(0,R1,'.',MarkerSize=20);
plot(0,R2,'.',MarkerSize=20);
plot(Nbar,Rbar,'.',MarkerSize=20);

text(0,R1,'\boldmath$(0,R_1)$',Interpreter='latex',FontSize=17);
text(0,R2,'\boldmath$\mathbf(0,R_2)$',Interpreter='latex',FontSize=17);
text(Nbar,Rbar,'\boldmath$\mathbf(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);

xlim([0 5])
xlabel('Consumer population density (N)');
ylabel('Resource population density (R)');
%{
legend('(N_0,R_0)=(5,1)','(N_0,R_0)=(5,2)','(N_0,R_0)=(5,3)','(N_0,R_0)=(1,5)' ...
    , '(N_0,R_0)=(2,5)','(N_0,R_0)=(3,5)','(N_0,R_0)=(0.1,R_2)');
legend('Location','westoutside')
%}
%legend(sprintf('s=%.2f',s0))

%% Bifurcation diagram

tspan=[0 5000];
sspan=linspace(2,6,50);
Nspan = sspan;
g=0;
R=d/(f*2);
N=(a*Rbar-(b+c*Rbar)*(1+Rbar))/(2*(1+Rbar));
for s = sspan
    g=g+1;
    [t,n] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s) , tspan, [N R]);
    Nspan(g)=n(end,1);
    N=n(end,1);
    R=n(end,2);
end
figure;
plot(sspan,Nspan,'LineWidth',1.2);
hold on;
xline(14/3,'LineWidth',1.2);
xlabel('s')
ylabel('N')