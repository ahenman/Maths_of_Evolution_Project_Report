%%Alice Henman%%
%%Chapter 2 Intro to adaptive dynamics%%

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

%% Initial evolution of population in time (simple)
r = 0.5;
d=0.2;
N0=5*rand;
tspan = [0 20];
[t,N] = ode45(@(t,N) N*(r-N*d), tspan, N0);
figure;
plot(t,N);
hold on;

plot(t,r./d*ones(size(t)),'--'); %predicted equilibrium value
title('dN/dt = N(r-Nd)')
xlabel('Time');
ylabel('Population (N)');
hold off;

%% Evolution of resident and mutant population in time (simple)
x=0.9;%rand; % resident trait
y=x+0.01; % mutant trait
N0=x/d;
tspan = [0 3000];

%ODEs for time evolution of population
function Ndot = Ndot(N,M,x,d)
    Ndot = N*(x-(N+M)*d);
end

function Mdot = Mdot(N,M,y,d)
    Mdot = M*(y-(N+M)*d);
end

function ddt = odefcn(n,x,y,d)
N = n(1);
M = n(2);
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = Ndot(N,M,x,d); % resident population
    ddt(2) = Mdot(N,M,y,d); % mutant population density
end
n0 = [N0,0.001]; %initial values, mutant population initially small
[t,n] = ode45(@(t,n) odefcn(n,x,y,d) , tspan, n0);
%{
%second mutation after resdient has equilibrated
y2=y+0.02;
N0=n(end,2);
n0 = [N0,0.001]; %initial values, mutant population initially small
[t2,n2] = ode45(@(t2,n2) odefcn(n2,x,y,d) , [500000 600000], n0);
%}
figure;
plot(t,n(:,1),'-' ,t,n(:,2),'--',LineWidth=2);
hold on;
%plot(t2,n2(:,1),'r--' ,t2,n2(:,2),'g:',LineWidth=2);

xlabel('Time');
ylabel('Population density');
%legend(sprintf('N_x, x = %.2f',x),sprintf('N_y, y = %.2f',y),'Location','northwest');
legend(sprintf('Resident population, x = %.2f',x),sprintf('Mutant Population, y = %.2f',y),'Location','east');
%legend(sprintf('Resident population, x = %.2f',x),sprintf('Mutant Population, y = %.2f',y), ...
%    sprintf('New Resident Population, x = %.2f',y),sprintf('Second Mutant Population, y = %.2f',y2),'Location','east');

%% PIP (simple)

%Invasion fitness
f = @(x,y) y-x;

y = linspace(0,1,500);
[X,Y] = meshgrid(y, y);
F = f(X,Y);

figure
hold on
% PIP
contourf(X,Y,F>0,1,'LineStyle','none')
contour(X,Y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
text(0.75,0.25,'-','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',50);
text(0.25,0.75,'+','HorizontalAlignment','left','VerticalAlignment','top','FontSize',50);
p1 = patch(NaN, NaN, [0.6 0.9 0.3]);  % green block for s>0
p2 = patch(NaN, NaN, [1 1 1]);  % gray block for s<0
hZero = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % black line for s=0

legend([p1 p2 hZero], {'f_x(y)>0 region', 'f_x(y)<0 region', 'f_x(y)=0 contour'});
%title(['PIP for σ_{α} = ',num2str(sigma)],'Interpreter','tex');

axis equal
box on

%% PIP with trait substitution sequence added

% Parameters
n_steps = 20; % number of mutations
step_size = 0.15;
initial_trait = 0.1;%x;

% Trait values and meshgrid
traits = linspace(0,1,200);
[X,Y] = meshgrid(traits, traits);

% Simulate TSS
[residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f);

% Compute PIP
F = f(X,Y);

% Plot PIP and TSS
figure;
hold on;
contourf(X,Y,F>0,1,'LineStyle','none');
contour(X,Y,F,[0 0],'k','LineWidth',2); %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0

scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k'); %plots initial value in black
scatter(residents(1:end-1), mutants, 110,'filled'); % Mutant steps
scatter(residents(2:end), residents(2:end), 110,'filled'); % Resident steps
scatter(residents(end), residents(end), 180, 'red','filled',...
    'MarkerEdgeColor','k');
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex');
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex');
%title('Trait Substitution Sequence on Pairwise Invasibility Plot');
legend('','','Starting trait value','Mutant trait','Resident trait', ...
    'final trait value','Location','northwest');
axis equal;
hold off;

% Compute the PIP
function [residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f)
    residents = zeros(n_steps+1, 1);
    mutants = zeros(n_steps, 1);
    residents(1) = initial_trait;
    for i = 1:n_steps
        res = residents(i);
        % Generate a single candidate mutant nearby
        mut = res + (2*rand-1) * step_size;
        mut = min(max(mut, traits(1)), traits(end));%makes sure mutant value is in range of the trait values
        mutants(i) = mut;
        % If mutant is successful, update resident; else, resident doesn't change
        if f(res, mut) > 0
            residents(i+1) = mut;
        else
            residents(i+1) = res;
        end
    end
end

%% Initial evolution of population in time (d varying)
a=0.5;%rand;
b=0.5;%rand;
r = 0.5;
d = @(r) a*r.^2+b;
function d=dfun(r,a,b)
    d=a*r.^2+b;
end
N0=5*rand;
tspan = [0 20];
[t,N] = ode45(@(t,N) N*(r-N*d(r)), tspan, N0);
figure;
plot(t,N);
hold on;

plot(t,r./d(r)*ones(size(t)),'--');
title('dN/dt = N(r-Nd)');
xlabel('Time');
ylabel('Population (N)');
hold off;

%% Evolution of resident and mutant population in time (d varying)
x=2*rand; % resident trait
y=x+0.2*rand-0.1; % mutant trait
N0=x/d(x);
tspan = [0 200];

%ODEs for time evolution of population
function Ndot = Ndot2(N,M,x,a,b)
    Ndot = N*(x-(N+M)*dfun(x,a,b));
end

function Mdot = Mdot2(N,M,y,a,b)
    Mdot = M*(y-(N+M)*dfun(y,a,b));
end

function ddt = odefcn2(n,x,y,a,b)
N = n(1);
M = n(2);
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = Ndot2(N,M,x,a,b); % resident population
    ddt(2) = Mdot2(N,M,y,a,b); % mutant population density
end
n0 = [N0,0.01]; %initial values, mutant population initially small
[t,n] = ode45(@(t,n) odefcn2(n,x,y,a,b) , tspan, n0);

figure;
plot(t,n(:,1),'-' ,t,n(:,2),'--',LineWidth=2);
xlabel('Time');
ylabel('Population density');
legend(sprintf('N_x, x = %.2f',x),sprintf('N_y, y = %.2f',y),'Location','northwest');

%% PIP (d varying)

%Invasion fitness
f = @(x,y) y-x.*d(y)./d(x);

z = linspace(0,3,500);
[X,Y] = meshgrid(z, z);
F = f(X,Y);

figure
hold on
% PIP
contourf(X,Y,F>0,1,'LineStyle','none')
contour(X,Y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex');
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex');
p1 = patch(NaN, NaN, [0.6 0.9 0.3]);  % green block for s>0
p2 = patch(NaN, NaN, [1 1 1]);  % white block for s<0
hZero = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % black line for s=0
%h=plot(sqrt(b/a)*ones(1,2),[0,Y(end)],'r--','LineWidth',2);

xline(sqrt(b/a),'--','f_x(y)<0','FontSize',20,'LabelOrientation','horizontal','LineWidth',2);
yline(sqrt(b/a),':','s(x)<0','FontSize',20,'LineWidth',2, ...
    'LabelHorizontalAlignment','center');
yline(sqrt(b/a),':','    s(x)>0','FontSize',20,'LineWidth',2, ...
    'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');

%legend([p1 p2 hZero h], {'f_x(y)>0 region', 'f_x(y)<0 region', 'f_x(y)=0 contour','x*'});

axis equal
box on

%% PIP with trait substitution sequence added (d varying)

% Parameters
n_steps = 25; % number of mutations
step_size = 0.15;
initial_trait = 0.3;%x;
%Colour scheme
G=[0.757 0.820 0.122
    0.596 0.788 0.169
    0.431 0.753 0.027
    0.318 0.694 0.020
    0.204 0.635 0.012
    0.224 0.580 0.027
    0.239 0.525 0.043
    0.122 0.455 0.051
    0.000 0.380 0.055];
% Trait values and meshgrid
traits = linspace(0,2,200);
[X,Y] = meshgrid(traits, traits);

% Simulate TSS
[residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f);

% Compute PIP
F = f(X,Y);

% Plot PIP and TSS
figure;
hold on;
contourf(X,Y,F>0,1,'LineStyle','none');
contour(X,Y,F,[0 0],'k','LineWidth',2); %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0

%xline(sqrt(b/a),'--','f_x(y)<0','LabelOrientation','horizontal','LineWidth',2);
%yline(sqrt(b/a),':','f_x(y)>0','LineWidth',2);

scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k'); %plots initial value in black
scatter(residents(1:end-1), mutants, 110,'filled'); % Mutant steps
scatter(residents(2:end), residents(2:end), 110,'filled'); % Resident steps
scatter(residents(end), residents(end), 180, 'red','filled',...
    'MarkerEdgeColor','k', 'LineWidth',2);
scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k'); %plots initial value in black
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex');
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex');
%title('Trait Substitution Sequence on Pairwise Invasibility Plot');
legend('','','Starting trait value','Mutant trait','Resident trait', ...
    'Final trait value','','Location','northwest');
axis equal;
hold off;