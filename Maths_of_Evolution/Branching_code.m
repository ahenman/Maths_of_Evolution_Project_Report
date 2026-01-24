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

%% (Mexican hat) Wavelet Competition function

% z represents x-y
z=linspace(-1.5,1.5,500);

%parameters
sigma = 0.25; %0.5;%width
sigma0 = 0.5;
beta0 = 0;   %symmetrical
beta1 = 1.5; %asymmetric

function alpha=alpha(x,y,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).* ...
(1.-((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2));
end
%The beta parameter can add asymmetry to the competition
%but I'm just keeping it as beta=0 for now

%z=x-y
figure;
%Plots of competition for 2 different widths
plot(z,alpha(z,0,sigma,beta0),'DisplayName', sprintf('σ_α = %.2f',sigma), ...
    "LineWidth",1.5);
hold on;
plot(z,alpha(z,0,sigma0,beta0),'r--','DisplayName', sprintf('σ_α = %.2f',sigma0), ...
    "LineWidth",1.5);
xlabel('Difference in trait value (x-y)');
ylabel('Strength of competition α(x − y)');
grid on;
ax=gca;
ax.FontSize = 20;
legend(FontSize=20);
legend('boxoff');

%% Double Gaussian Carrying capacity K(x)
function Gauss=G(sigma,x0,x)
    Gauss = exp(-(x-x0).^2./(2*sigma^2));
end

function K=K(C1,C2,sigma1,sigma2,x0,x)
    K = C1*G(sigma1,x0,x)-C2*G(sigma2,x0,x);
end

%Parameters
C1 = 2.5;
C2 = 2;
sigma1 = 0.3;%0.5;%
sigma2 = 0.15;%0.35;%
x0 = 1.5;

x=linspace(0,3,500);

figure;
plot(x,K(C1,C2,sigma1,sigma2,x0,x),'r--',"LineWidth",1.5);
hold on;
xlabel('Trait value x');
ylabel('Carrying Capacity, K(x)');
grid on;
ax=gca;
ax.FontSize = 20;
legend("boxoff");

%% The Mathssss (The PIP)

%Parameters
r=rand+0.001; %growth rate
%stable strategies
x2=x0+sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));
x3=x0-sqrt(2*sigma1^2*sigma2^2*log(C1*sigma2^2/(C2*sigma1^2))/(sigma2^2-sigma1^2));

h1=plot(x0*ones(1,2),[0,1.5],'k-');
h2=plot(x2*ones(1,2),[0,1.5],'r-.');
h3=plot(x3*ones(1,2),[0,1.5],'b--');
legend([h1,h2,h3], {sprintf('x_0 = %.2f', x0) sprintf('x_0+x_1 = %.2f', x2) sprintf('x_0-x_1 = %.2f', x3)});

%Invasion fitness
f = @(x,y) r*(1-alpha(y,x,sigma,beta0).*K(C1,C2,sigma1,sigma2,x0,x) ...
    ./K(C1,C2,sigma1,sigma2,x0,y));

y = linspace(0.5,2.5,500);
[x,y] = meshgrid(y, y);
F = f(x,y);

figure
hold on
% PIP
contourf(x,y,F>0,1,'LineStyle','none')
contour(x,y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',17);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',17);
title(['PIP for σ_{α} = ',num2str(sigma)],'Interpreter','tex');

h1=plot(x0*ones(1,2),[0.5,y(end)],'k-');
h2=plot(x2*ones(1,2),[0.5,y(end)],'r-.');
h3=plot(x3*ones(1,2),[0.5,y(end)],'b--');

legend(Location="northwest");
legend([h1,h2,h3], {sprintf('x_0 = %.2f', x0) sprintf('x_0+x_1 = %.2f', x2) sprintf('x_0-x_1 = %.2f', x3)});
axis equal
box on

%% Plot to show that the conditions for branching are satisfied
figure;
z=linspace(0.001,1,1000); %range of sigmak1 and k2
[sk1,sk2] = meshgrid(z,z);
C=C1/C2; %ratio of C1/C2>1
h = @(x,y) 3*(y.^2-x.^2)./(2.*log(C*y.^2./x.^2));
H=h(sk1,sk2); %values of the RHS of branching condition
[h,w]=size(H);
for i = 2:h
    for j = 1:i
        if i>=j
            H(i,j)=-3;%makes sure all values with sigmak1<k2 are undefined
        end
    end
end
hlims=[-3 3]; %stops large values drowning out small variations close to 0
imagesc(z,z,H,hlims);
hold on;
cmap = [abyss;autumn];%makes difference between positive and negative obvious
colormap(cmap); 
colorbar;
%Marks all negatives as undefined, makes it obvious that values of "3"
%could be bigger
colorbar('Ticks',[-1.5,0,0.5,1,1.5,2,2.5,3],...
         'TickLabels',{'undefined','0','0.5','1','1.5','2','2.5','3+'})
xlabel("\sigma_{k1}",'Interpreter','tex');
ylabel("\sigma_{k2}",'Interpreter','tex');
set(gca,'YDir','normal'); %makes axis look normal
% title('C* colour plot on \sigma_{k1},\sigma_{k2} axis','Interpreter','tex');

%% Example PIP to show branching conditions
Y = linspace(1,1.48,500);
X = linspace(1,2*x3-1,500);
[x,y] = meshgrid(Y, Y);
F = f(x,y);

figure
hold on

contourf(x,y,F>0,1,'LineStyle','none')
contour(x,y,F,[0 0],'k','LineWidth',2) %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',17);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',17);
%title(['PIP for σ_{α} = ',num2str(sigma)],'Interpreter','tex');

h1=plot(x3*ones(1,2),[1,y(end)],'k-');
h2=plot((x3+0.05)*ones(1,2),[1,x3-0.05],'b:');
h3=plot((x3-0.05)*ones(1,2),[1,x3+0.05],'r-.');
h4=plot([1,x3-0.05],(x3+0.05)*ones(1,2),'b:');
h5=plot([1,x3+0.05],(x3-0.05)*ones(1,2),'r-.');
h6=plot(X,2*x3-X,'--',LineWidth=1.5);

plot(x3+0.05,x3-0.05,'k.','MarkerSize',18);
plot(x3-0.05,x3+0.05,'k.','MarkerSize',18);
text(x3+0.05,x3-0.05,"$f_{x'}(y')>0$",'Interpreter','latex','FontSize',18,'VerticalAlignment','bottom');
text(x3-0.05,x3+0.05,"$f_{x'}(y')>0$",'Interpreter','latex','FontSize',18,'HorizontalAlignment','right','VerticalAlignment','top');

%legend(Location="northwest");
legend([h1,h2,h3,h6], {sprintf('$x^*$') sprintf("$x'$") sprintf("$y'$") sprintf("$y=2x^*-x$")},'Interpreter','latex');
axis equal
set(gca,'XTick',[], 'YTick', [])
box on

%% PIP with trait substitution sequence added

% Parameters
n_steps = 20; % number of mutations
step_size = 0.15;
initial_trait = 1;

% Trait values and meshgrid
traits = linspace(1,2,200);
[X,Y] = meshgrid(traits, traits);

% Simulate TSS
[residents, mutants] = simulate_tss_single_mutant(initial_trait, n_steps, step_size, traits, f);

% Compute PIP
pip = f(X, Y);

% Plot PIP and TSS
figure;
hold on;
scatter(residents(1), residents(1), 250, 'k','filled',...
    'MarkerEdgeColor','k');
contour(X, Y, pip, [0 0], 'LineWidth', 2);

%contour(X, Y, pip, [-10 0 10], 'LineColor', 'none', 'FaceAlpha',0.7);
scatter(residents(1:end-1), mutants, 110,'filled'); % Mutant steps
scatter(residents, residents, 110,'filled'); % Resident steps
scatter(residents(end), residents(end), 180, 'red','filled',...
    'MarkerEdgeColor','k');
xlabel('Resident trait');
ylabel('Mutant trait');
title('Trait Substitution Sequence on Pairwise Invasibility Plot');
legend('Starting trait value','PIP (zero invasion fitness)','Mutant trait','Resident trait',...
    'final trait value','Location','northeast');
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
        mut = min(max(mut, traits(1)), traits(end));
        mutants(i) = mut;
        % If mutant is successful, update resident; else, resident doesn't change
        if f(res, mut) > 0
            residents(i+1) = mut;
        else
            residents(i+1) = res;
        end
    end
end