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

%% error function A and K
a=3;
b=2;
e=linspace(0,6,500);
figure;
plot(e,(erf(e-a)+1)/2,LineWidth=2);
hold on;
plot(e,(erf(e-b)+1)/2,LineWidth=2);
xlabel('efficiency (e)')
legend('A(e)','K(e)');

%% Allee Population model

r=5;
A = @(e) (erf(e-a)+1)/2;
K = @(e) (erf(e-b)+1)/2;
e=3;

Ndot = @(N,e) r*N.*(N/A(e)-1).*(1-N/K(e));

figure;
N=linspace(0,1,500);
plot(N,Ndot(N,e), LineWidth=2); %plots dN/dt by N
hold on;
xlabel('\textsf{Population Density} ($N$)',Interpreter='latex');
ylabel('$\frac{dN}{dt}$', Interpreter='latex',Rotation=0,FontSize=25);
yline(0);                       %x-axis
%adds A and K to x-axis
%xtickVals = unique([0 : 0.2 : 1, A(e),K(e)]);
%xtickLabs = compose('%.3g',xtickVals);
%xtickLabs(ismembertol(xtickVals,A(e))) = {'A'};
%xtickLabs(ismembertol(xtickVals,K(e))) = {'K'};
%set(gca,'xtick', xtickVals, 'xticklabel', xtickLabs, 'xlim', [min(xtickVals),max(xtickVals)])
A=xline(A(e),'r--','DisplayName', sprintf('A(e) = %.3f',A(e)));
K=xline(K(e),'r-.','DisplayName', sprintf('K(e) = %.3f',K(e)));
legend([A K],'Location','southwest')
%plot(0,0,'k.',K(e),0,'k.',MarkerSize=25);
%plot(A(e),0,'ko',MarkerSize=10,LineWidth=1.2);

%% Evolution of resident and mutant population in time

A = @(e) (erf(e-a)+1)/2;
K = @(e) (erf(e-b)+1)/2;

function Afun= Afun(e,a)
Afun = (erf(e-a)+1)/2;
end
function Kfun= Kfun(e,b)
Kfun = (erf(e-b)+1)/2;
end
x=5; % resident trait
y=x+0.001*rand; % mutant trait
N0=Kfun(x,b);
tspan = linspace(0,5000000,10000000);%[0 5000000];

%ODEs for time evolution of population
function Xdot = Xdot(N,M,x,r,a,b)
    Xdot = r*N.*((N+M)/Afun(x,a)-1).*(1-(N+M)/Kfun(x,b));
end

function Ydot = Ydot(N,M,y,r,a,b)
    Ydot = r*M.*((N+M)/Afun(y,a)-1).*(1-(N+M)/Kfun(y,b));
end

function ddt = odefcn(n,x,y,r,a,b)
N = n(1);
M = n(2);
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = Xdot(N,M,x,r,a,b); % resident population
    ddt(2) = Ydot(N,M,y,r,a,b); % mutant population density
end
n0 = [N0,0.01]; %initial values, mutant population initially small
[t,n] = ode45(@(t,n) odefcn(n,x,y,r,a,b) , tspan, n0);
%[t1,n1] = ode45(@(t1,n1) odefcn(n1,x,y,r,a,b) , [40001 50000], [n(end,1) n(end,2)-0.005]);

%T=[t;t1];
%N=abs([n;n1]);

figure;
plot(t,abs(n(:,1)),'-' ,t,abs(n(:,2)),'--',LineWidth=2.5);
%plot(T,abs(N(:,1)),'-' ,T,abs(N(:,2)),'--',LineWidth=2.5);
xlabel('Time');
ylabel('Population density');
legend(sprintf('N_x, x = %.4f',x),sprintf('N_y, y = %.4f',y),'Location','northwest');

%% PIP

%Invasion fitness
f = @(x,y) r*(K(x)./A(y)-1).*(1-K(x)./K(y));

y = linspace(0,6,500);
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
%text(0.75,0.25,'-','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',50);
%text(0.25,0.75,'+','HorizontalAlignment','left','VerticalAlignment','top','FontSize',50);
p1 = patch(NaN, NaN, [0.6 0.9 0.3]);  % green block for s>0
p2 = patch(NaN, NaN, [1 1 1]);  % gray block for s<0
hZero = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % black line for s=0

legend([p1 p2 hZero], {'f_x(y)>0 region', 'f_x(y)<0 region', 'f_x(y)=0 contour'});
%title(['PIP for σ_{α} = ',num2str(sigma)],'Interpreter','tex');

axis equal
box on

%% PIP with trait substitution sequence added

% Parameters
n_steps = 500; % number of mutations
step_size = 0.15;
initial_trait = 1;%x;

% Trait values and meshgrid
traits = linspace(0,6,2000);
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