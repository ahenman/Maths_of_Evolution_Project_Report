%%Alice Henman%%
%%Chapter 3 Evolutionary Suicide%%
%%Classifying singular strategies%%

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
set(0,'defaultLineLineWidth',2)

%% Resource consumer model N by R
tspan=[0 500];%linspace(0,500,100000);%
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
sstar=d*(c+sqrt(a*c))/(f*(a-c));
s0=2;
R1=(a-b-c-sqrt((b+c-a)^2-4*b*c))/(2*c);
R2=(a-b-c+sqrt((b+c-a)^2-4*b*c))/(2*c);
Rbar=d/(f*s0);
Nbar=(a*Rbar-(b+c*Rbar)*(1+Rbar))/(s0*(1+Rbar));

[t,n11] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 1]);
[t,n12] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 2]);
[t,n13] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 3]);
[t,n14] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [3 5]);
[t,n15] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [2 5]);
[t,n16] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [1 5]);
%[t,n17] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [0.1 R2]);

s0=5;

[t,n21] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 1]);
[t,n22] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 2]);
[t,n23] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [5 3]);
[t,n24] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [3 5]);
[t,n25] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [2 5]);
[t,n26] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [1 5]);
%[t,n27] = ode45(@(t,n) odefcn1(n,a,b,c,d,f,s0) , tspan, [0.1 R2]);

figure;
G=[0.757 0.820 0.122
    0.596 0.788 0.169
    0.431 0.753 0.027
    0.318 0.694 0.020
    0.204 0.635 0.012
    0.224 0.580 0.027
    0.239 0.525 0.043
    0.122 0.455 0.051
    0.000 0.380 0.055];
%{
plot(abs(n1(:,1)),abs(n1(:,2)),LineWidth=2.5,Color=G(9,:));
hold on;
plot(abs(n2(:,1)),abs(n2(:,2)),LineWidth=2.5,Color=G(7,:));
plot(abs(n3(:,1)),abs(n3(:,2)),LineWidth=2.5,Color=G(5,:));
plot(abs(n4(:,1)),abs(n4(:,2)),LineWidth=2.5,Color=G(4,:));
plot(abs(n5(:,1)),abs(n5(:,2)),LineWidth=2.5,Color=G(3,:));
plot(abs(n6(:,1)),abs(n6(:,2)),LineWidth=2.5,Color=G(2,:));
%plot(abs(n7(:,1)),abs(n7(:,2)),LineWidth=2.5,Color=G(1,:));
%equilibria:
plot(0,0,'.',MarkerSize=20);
%plot(0,R1,'.',MarkerSize=20);
%plot(0,R2,'.',MarkerSize=20);
plot(Nbar,Rbar,'.',MarkerSize=20);

%text(0,R1,'\boldmath$(0,R_1)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');
%text(0,R2,'\boldmath$(0,R_2)$',Interpreter='latex',FontSize=17);
text(Nbar,Rbar,'\boldmath$(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);
text(0,0,'\boldmath$(0,0)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');

axis equal
xlim([0 5])
ylim([0 5])
xlabel('Consumer population density (N)','FontSize',22);
ylabel('Resource population density (R)','FontSize',22);
%{
legend('(N_0,R_0)=(5,1)','(N_0,R_0)=(5,2)','(N_0,R_0)=(5,3)','(N_0,R_0)=(1,5)' ...
    , '(N_0,R_0)=(2,5)','(N_0,R_0)=(3,5)','(N_0,R_0)=(0.1,R_2)');
legend('Location','westoutside')
%}
%legend(sprintf('s=%.2f',s0))
%}
%% Create plots side-by-side
tile = tiledlayout(1,2);

ax1 = nexttile;
plot(ax1,abs(n11(:,1)),abs(n11(:,2)),LineWidth=2.5,Color=G(9,:));
hold on;
plot(ax1,abs(n12(:,1)),abs(n12(:,2)),LineWidth=2.5,Color=G(7,:));
plot(ax1,abs(n13(:,1)),abs(n13(:,2)),LineWidth=2.5,Color=G(5,:));
plot(ax1,abs(n14(:,1)),abs(n14(:,2)),LineWidth=2.5,Color=G(4,:));
plot(ax1,abs(n15(:,1)),abs(n15(:,2)),LineWidth=2.5,Color=G(3,:));
plot(ax1,abs(n16(:,1)),abs(n16(:,2)),LineWidth=2.5,Color=G(2,:));
%plot(ax1,abs(n17(:,1)),abs(n17(:,2)),LineWidth=2.5,Color=G(1,:));
%equilibria:
plot(ax1,0,0,'.',MarkerSize=20);
%plot(ax1,0,R1,'.',MarkerSize=20);
%plot(ax1,0,R2,'.',MarkerSize=20);
plot(ax1,Nbar,Rbar,'.',MarkerSize=20);

%text(ax1,0,R1,'\boldmath$(0,R_1)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');
%text(ax1,0,R2,'\boldmath$(0,R_2)$',Interpreter='latex',FontSize=17);
text(ax1,Nbar,Rbar,'\boldmath$(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);
text(ax1,0,0,'\boldmath$(0,0)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');

axis equal
xlim(ax1,[0 5])
ylim(ax1,[0 5])
title(ax1,'Stable coexistence','FontSize',20);

ylabel(ax1,'Resource density (R)','FontSize',22);

ax2 = nexttile;
plot(ax2,abs(n21(:,1)),abs(n21(:,2)),LineWidth=2.5,Color=G(9,:));
hold on;
plot(ax2,abs(n22(:,1)),abs(n22(:,2)),LineWidth=2.5,Color=G(7,:));
plot(ax2,abs(n23(:,1)),abs(n23(:,2)),LineWidth=2.5,Color=G(5,:));
plot(ax2,abs(n24(:,1)),abs(n24(:,2)),LineWidth=2.5,Color=G(4,:));
plot(ax2,abs(n25(:,1)),abs(n25(:,2)),LineWidth=2.5,Color=G(3,:));
plot(ax2,abs(n26(:,1)),abs(n26(:,2)),LineWidth=2.5,Color=G(2,:));
%plot(ax2,abs(n27(:,1)),abs(n27(:,2)),LineWidth=2.5,Color=G(1,:));
%equilibria:
plot(ax2,0,0,'.',MarkerSize=20);
%plot(ax2,0,R1,'.',MarkerSize=20);
%plot(ax2,0,R2,'.',MarkerSize=20);
plot(ax2,Nbar,Rbar,'.',MarkerSize=20);

%text(ax2,0,R1,'\boldmath$(0,R_1)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');
%text(ax2,0,R2,'\boldmath$(0,R_2)$',Interpreter='latex',FontSize=17);
text(ax2,Nbar,Rbar,'\boldmath$(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);
text(ax2,0,0,'\boldmath$(0,0)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');

axis equal
xlim(ax2,[0 5])
ylim(ax2,[0 5])
title(ax2,'Mutual extinction','FontSize',20);

% Link the axes
linkaxes([ax1,ax2],'y');

% Add shared title and axis labels
xlabel(['Consumer population density (N)                                    ' ...
    '            '],'FontSize',22);

% Move plots closer together
yticklabels(ax2,{});
xticks(ax1,[0:5]);
xticks(ax2,[0:5]);
tile.TileSpacing = 'compact';

%% EvoSuicide video

%{
 vidfile = VideoWriter('Suicidespiral(1).mp4','MPEG-4');
 vidfile.FrameRate=10;
 open(vidfile);
 f=figure;
 for j = 1:100
    f.Name = ['Simulation time: t = ', num2str((j-1))];
    %Plotting
    tile = tiledlayout(1,2);

    ax1 = nexttile;
    plot(ax1,abs(n11([1:j+1],1)),abs(n11([1:j+1],2)),LineWidth=2.5,Color=G(9,:));
    hold on;
    plot(ax1,abs(n12([1:j+1],1)),abs(n12([1:j+1],2)),LineWidth=2.5,Color=G(7,:));
    plot(ax1,abs(n13([1:j+1],1)),abs(n13([1:j+1],2)),LineWidth=2.5,Color=G(5,:));
    plot(ax1,abs(n14([1:j+1],1)),abs(n14([1:j+1],2)),LineWidth=2.5,Color=G(4,:));
    plot(ax1,abs(n15([1:j+1],1)),abs(n15([1:j+1],2)),LineWidth=2.5,Color=G(3,:));
    plot(ax1,abs(n16([1:j+1],1)),abs(n16([1:j+1],2)),LineWidth=2.5,Color=G(2,:));
    %equilibria:
    plot(ax1,0,0,'.',MarkerSize=20);
    plot(ax1,Nbar,Rbar,'.',MarkerSize=20);

    text(ax1,Nbar,Rbar,'\boldmath$(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);
    text(ax1,0,0,'\boldmath$(0,0)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');
    
    axis equal
    xlim(ax1,[0 5])
    ylim(ax1,[0 5])
    
    ylabel(ax1,'Resource density (R)','FontSize',22);
    
    ax2 = nexttile;
    plot(ax2,abs(n21([1:j+1],1)),abs(n21([1:j+1],2)),LineWidth=2.5,Color=G(9,:));
    hold on;
    plot(ax2,abs(n22([1:j+1],1)),abs(n22([1:j+1],2)),LineWidth=2.5,Color=G(7,:));
    plot(ax2,abs(n23([1:j+1],1)),abs(n23([1:j+1],2)),LineWidth=2.5,Color=G(5,:));
    plot(ax2,abs(n24([1:j+1],1)),abs(n24([1:j+1],2)),LineWidth=2.5,Color=G(4,:));
    plot(ax2,abs(n25([1:j+1],1)),abs(n25([1:j+1],2)),LineWidth=2.5,Color=G(3,:));
    plot(ax2,abs(n26([1:j+1],1)),abs(n26([1:j+1],2)),LineWidth=2.5,Color=G(2,:));

    %equilibria:
    plot(ax2,0,0,'.',MarkerSize=20);
    plot(ax2,Nbar,Rbar,'.',MarkerSize=20);
    
    text(ax2,Nbar,Rbar,'\boldmath$(\bar{N},\bar{R})$',Interpreter='latex',FontSize=17);
    text(ax2,0,0,'\boldmath$(0,0)$',Interpreter='latex',FontSize=17,VerticalAlignment='bottom');
    
    axis equal
    xlim(ax2,[0 5])
    ylim(ax2,[0 5])
    
    % Link the axes
    linkaxes([ax1,ax2],'y');
    
    % Add shared title and axis labels
    xlabel(['Consumer population density (N)                                    ' ...
        '            '],'FontSize',22);
    
    % Move plots closer together
    yticklabels(ax2,{})
    xticks(ax1,[0:5]);
    xticks(ax2,[0:5]);
    tile.TileSpacing = 'compact';
    writeVideo(vidfile, getframe(gcf));
 end
close(vidfile)
%}
%% PIP

%Invasion fitness
fit = @(x,y) d*(y./x-1);

z = linspace(0,6,500);
[X,Y] = meshgrid(z, z);
F = fit(X,Y);

figure;
hold on
% PIP
contourf(X,Y,F>0,1,'LineStyle','none')
contour(X,Y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',25);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',25);
p1 = patch(NaN, NaN, [0.6 0.9 0.3]);  % green block for s>0
p2 = patch(NaN, NaN, [1 1 1]);  % white block for s<0
hZero = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % black line for s=0

xline(d/(R2*f));
xline(sstar);
z = ones(size(Y));
area([0 d/(R2*f)],[6;6],"FaceColor",G(9,:),"FaceAlpha",0.5);
area([sstar 6],[6;6],"FaceColor",G(9,:),"FaceAlpha",0.5);
text(d/(R2*f),6,'$\partial V$','Interpreter','latex','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',25);
text(sstar,6,'$\partial V$','Interpreter','latex','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',25);

axis equal
box on

%% Resource consumer model N by R
tspan=[0 500];%linspace(0,500,100000);%
x=2;
y=x+0.2*rand-0.1;

function dRdt2 = dRdt2(N,M,R,a,b,c,x,y)
    dRdt2 = R*(a*R/(1+R)-b-c*R-x*N-y*M);
end

function ddt2 = odefcn2(n,a,b,c,d,f,x,y)
N = n(1);
M = n(2);
R = n(3);
ddt2 = zeros(3,1); % Initialize the derivative vector
    ddt2(1) = dNdt(N,R,d,f,x); % resident population
    ddt2(2) = dNdt(M,R,d,f,y); % mutant population density
    ddt2(3) = dRdt2(N,M,R,a,b,c,x,y); %Resource set by resident
end
R0=d/(f*x);
N0=[(a*R0-(b+c*R0)*(1+R0))/(x*(1+R0)),0.1,R0];
[t,n] = ode45(@(t,n) odefcn2(n,a,b,c,d,f,x,y) , tspan, N0);
figure;
plot(t,n(:,1),'-o',t,n(:,2),'-o',t,n(:,3),'-o');
legend(['N,x=',num2str(x)],['M,y=',num2str(y)],'R');

%% simulation
%{
tspan=[0 50000];    %timespan for each simulation
T=300;              %Number of simulations
x=2;                %initial trait value
%initial densities
R0=d/(f*x);
N0=(a*R0-(b+c*R0)*(1+R0))/(x*(1+R0));
%Storage for data
Rplot=zeros(T+1,1);
Nplot=zeros(T+1,1);
splot=zeros(T+1,1);
Tplot=zeros(T+1,1);
%initial values
Rplot(1)=R0;
Nplot(1)=N0;
splot(1)=x;
Tplot(1)=0;
for time=2:T+1
    y=x+0.1*rand-0.05; %Generate mutant (+- 0.05 from resident)
    %run time evolution of mutant, resident and resource
    [t,n] = ode45(@(t,n) odefcn2(n,a,b,c,d,f,x,y) , tspan, [N0,0.1,R0]);
    %if unsuccesful population hasn't died
    if abs(min(n(end,1),n(end,2))) > 1e-13
        %extend simulation
        [t,n] = ode45(@(t,n) odefcn2(n,a,b,c,d,f,x,y) , tspan, ...
            [n(end,1),n(end,2),n(end,3)]);
    end
    %Make initial values for next simulation follow on from final values
    R0=n(end,3);
    N0=max(n(end,1),n(end,2));
    %add results of simulation to data
    Rplot(time)=R0;
    Nplot(time)=N0;
    Tplot(time)=time-1;
    %update trait value to be the fittest
    if fit(x,y)>0
        x=y;
    end
    splot(time)=x;
end
%to find the time that s is closest to s*
Tstar = find(abs(splot-sstar)==min(abs(splot-sstar)),1)-1;
Text = find(Rplot==0,1)-1;
C = orderedcolors("gem");
figure;
plot(Tplot,Nplot,'Color',C(1,:));
xline(Tstar,'--','LineWidth',1.5);
xline(Text,'r--','LineWidth',1.5);
xlabel('Time');
ylabel('Population density');

figure;
plot(Tplot,Rplot,'Color',C(2,:));
xline(Tstar,'--','LineWidth',1.5);
xline(Text,'r--','LineWidth',1.5);
xlabel('Time');
ylabel('Resource density');

figure;
plot(Tplot,splot,'Color',C(3,:));
yline(sstar,'k--','LineWidth',1.5,'Label','s=s^*','Interpreter','tex', ...
    'FontSize',18,'LabelVerticalAlignment','bottom');
line([Tstar Tstar], [0 sstar],'LineStyle','--','LineWidth',1.5,'Color','k');
yline(splot(Text),'r--','LineWidth',1.5,'Label',['s=',num2str(splot(Text),3)], ...
    'Interpreter','tex','FontSize',18,'LabelVerticalAlignment','top');
line([Text Text], [0 splot(Text)],'LineStyle','--','LineWidth',1.5,'Color','r');
xlabel('Time');
ylabel('Trait value');

%}