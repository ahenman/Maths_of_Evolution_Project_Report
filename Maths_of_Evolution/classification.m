%%Alice Henman%%
%%Chapter 2 Intro to adaptive dynamics%%
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

%% fitness change by x
x=linspace(0,5,500);
xstar=1.7;
f=(x-xstar).^2;
figure;

plot(x,f);
hold on;
plot(xstar*ones(1,2),[0,f(end)],'--');
xlim([0,f(end)])
axis([0 x(end) 0 x(end)])
xtickVals = unique([0 : 1 : 5, xstar]);
xtickLabs = compose('%.3g',xtickVals);
xtickLabs(ismembertol(xtickVals,xstar)) = {'x*'}; 
set(gca,'xtick', xtickVals, 'xticklabel', xtickLabs, 'xlim', [min(xtickVals),max(xtickVals)])
%set(gca,xtickLabs(ismembertol(xtickVals,xstar)),'$x^*$',"TickLabelInterpreter",'latex')
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Invasion Fitness}$ ($f_x(y)$)','Interpreter','latex','FontSize',20);
axis square

%% unstable divergent strategy e.g. PIP

xbar=1.5;
%invasion fitness
f = @(x,y) (y-x).*(y+x-2*xbar);

z = linspace(0,3,500);
[X,Y] = meshgrid(z, z);
F = f(X,Y);

figure
hold on
% PIP
contourf(X,Y,F>0,1,'LineStyle','none')
contour(X,Y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);

xline(xbar,'--','f_x(y)>0','FontSize',20,'LabelOrientation','horizontal');
yline(xbar,':','f_x(y)<0','FontSize',20);

axis equal
box on

%% PIP with trait substitution sequence added

% Parameters
n_steps = 20; % number of mutations
step_size = 0.2;
initial_trait = xbar;%x;

% Trait values and meshgrid
traits = linspace(0,3,200);
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