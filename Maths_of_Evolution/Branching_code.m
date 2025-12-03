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

%Gaussian Competition function: alpha(x-y)
function alpha = alpha(z,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).*exp(-(z+sigma.^2.*beta).^2./(2.*sigma.^2));
end
% z represents x-y
z=linspace(-3,3,500);

%parameters
sigma = 0.5; %width
sigma0 = 1;
beta0 = 0;   %symmetrical
beta1 = 1.5; %asymmetric

figure;
plot(z,alpha(z,sigma,beta0),'DisplayName', sprintf('β = %.1f',beta0), ...
    'LineWidth',1.2);
hold on;    %Plotting 2 plots on 1
plot(z,alpha(z,sigma,beta1),"--",'DisplayName', sprintf('β = %.1f',beta1), ...
    'LineWidth',1.2);  %To differentiate curves in black and white
xlabel('Difference in trait value (x-y)');
ylabel('Strength of competition α(x − y)');
ax=gca;
ax.FontSize = 18;    %sets axes labels and titles size
legend(FontSize=18); %sets legend fontsize
legend('boxoff');    %Removes box around legend



%(Mexican hat) Wavelet Competition function
function alpha1=alpha1(x,y,sigma,beta)
    alpha1 = exp(sigma.^2.*beta.^2/2).* ...
(1.-((x-y+sigma.^2.*beta)./sigma).^2).*exp(-(x-y+sigma.^2.*beta).^2./(2*sigma^2));
end
%The beta parameter can add asymmetry to the competition
%but I'm just keeping it as beta=0 for now
%z=x-y

figure;
plot(z,alpha1(z,0,sigma,beta0),'DisplayName', sprintf('σ_α = %.2f',sigma), ...
    "LineWidth",1.5);
hold on;
plot(z,alpha1(z,0,sigma0,beta0),'--','DisplayName', sprintf('σ_α = %.2f',sigma0), ...
    "LineWidth",1.5);
xlabel('Difference in trait value (x-y)');
ylabel('Strength of competition α(x − y)');
grid on;
ax=gca;
ax.FontSize = 18;
legend(FontSize=18);
legend('boxoff');

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
sigma1 = 0.3;%0.5*rand;
sigma2 = 0.15;%sigma1*sqrt(C2/C1)*rand; % condition for there to be 2 humps
x0 = 1.5;%1+5*rand;

x=linspace(0,3,500);

figure;
plot(x,K(C1,C2,sigma1,sigma2,x0,x),"LineWidth",1.5);
hold on;
xlabel('Trait value x');
ylabel('Carrying Capacity, K(x)');
grid on;
ax=gca;
ax.FontSize = 18;

%The Mathssss

%Invasion fitness
r=rand+0.001;

f = @(x,y) r*(1-alpha1(y,x,sigma,beta0)*K(C1,C2,sigma1,sigma2,x0,x) ...
    /K(C1,C2,sigma1,sigma2,x0,y));

y = linspace(0,10,500);
[x,y] = meshgrid(y, y);
F = f(x,y);

figure
hold on

contourf(x,y,F>0,1,'LineStyle','none')
% Add the zero-level curve
contour(x,y,F,[0 0],'k','LineWidth',2)
% Cosmetics
colormap([1 1 1; 0.6 0.9 0.3]);
xlabel('x'), ylabel('y');
title('PIP');

p1 = patch(NaN, NaN, [0.6 0.9 0.3]);  % green block for s>0
p2 = patch(NaN, NaN, [1 1 1]);  % gray block for s<0
hZero = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % black line for s=0

h1=plot(x0*ones(1,2),[0,y(end)],'r--');

legend([p1 p2 hZero h1], {'s>0 region', 's<0 region', 's=0 contour',sprintf('x^∗ = %.2f', x0)});
axis equal
box on