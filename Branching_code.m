%%Alice Henman%%
%%Evolutionary Branching in Models for Symmetric and Asymmetric
%%Competition%%

%Gaussian Competition function: alpha(x-y)
function alpha = alpha(z,sigma,beta)
    alpha = exp(sigma.^2.*beta.^2/2).*exp(-(z+sigma.^2.*beta).^2./(2.*sigma.^2));
end
% z represents x-y
z=linspace(-3,3,500);

%parameters
sigma = 0.65;%width
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
function alpha1=alpha1(z,sigma,beta)
    alpha1 = exp(sigma.^2.*beta.^2/2).*(2/(sqrt(3*sigma)*pi^(1/4)).* ...
        (1.-(z./sigma).^2).*exp(-(z+sigma.^2.*beta).^2./(2*sigma^2)));
end
%Parameters
beta2 = 1; %Slightly less skewed than for Gaussian

figure;
plot(z,alpha1(z,sigma,beta0),'DisplayName', sprintf('β = %.1f',beta0), ...
    "LineWidth",1.2);
hold on;
plot(z,alpha1(z,sigma,beta2),'--','DisplayName', sprintf('β = %.1f',beta2), ...
    "LineWidth",1.2);
xlabel('Difference in trait value (x-y)');
ylabel('Strength of competition α(x − y)');
grid on;
ax=gca;
ax.FontSize = 18;
legend(FontSize=18);
legend('boxoff');

%(Mexican hat) Carrying capacity K(x)

function K=K(C1,C2,sigma1,sigma2,beta1,beta2,x)
    K = C1*alpha(x,sigma1,beta1)-C2*alpha(x,sigma2,beta2);
end

%Parameters
C1 = 2;
C2 = 1;

x=linspace(-3,3,500);

figure;
plot(z,K(C1,C2,sigma,sigma,beta0,beta0,x),"LineWidth",1.2);
hold on;
xlabel('Trait value x');
ylabel('Carrying Capacity, K(x)');
grid on;
ax=gca;
ax.FontSize = 18;


