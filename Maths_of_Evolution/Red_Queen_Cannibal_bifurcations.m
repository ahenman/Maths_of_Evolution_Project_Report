%%Alice Henman%%
%%Chapter 4 Red Queen%%

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
set(0,'defaultLineLineWidth',2)

C = orderedcolors("gem12");
C(2,:)=[];
C(3,:)=[];
C(6,:)=[];

function A0= A0(A,x0,a,x)
A0 = 2*A*x.^a./((x.^2/x0).^a+(x0)^a);
end
function Aij= Aij(A,xbar,xunder,b,g,d,p,x,y)
Aij = A*2./((p*x./y).^(b)+(y./(p*x)).^(b)).*(x.^g./(xunder^g+x.^g)).*(xbar^d./(xbar^d+x.^d));
end
function hij= h(w1,w2,x)
hij = w1*x.^(-w2);
end

x0 = 0.1;
xbar = 5;
xunder = 0.5;
Ao = 1;
A=10;
a = 2;
%b=1.5;%2.5;%
g=8;
d=2;
c=1;
f=0.6;
p=0.2;
w1=0.1;
w2=0.25;
%N0=10;%500;%200;%

%Non-trivial equilibrium
function Nbar= Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
a10=A0(Ao,x0,a,x);
a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
h1=h(w1,w2,x);

B=-c-c*N0*h1.*a10+(f-1)*a11;
A1=c*h1.*a11;
C=f*N0*a10;

Nbar = 2*C./(sqrt(B.^2 +4*C.*A1)-B);%(B+ sqrt(B.^2 +4*C.*A1))./(2*A1);%
end

%Invasion fitness
function fit = fit(x,y,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x);
    a20=A0(Ao,x0,a,y);
    a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
    a21=Aij(A,xbar,xunder,b,g,d,p,y,x);
    a12=Aij(A,xbar,xunder,b,g,d,p,x,y);
    h1=h(w1,w2,x);
    h2=h(w1,w2,y);
    N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

    fit=f*(a20*N0+a21.*N)./(1+h2.*(a20*N0+a21.*N)) -a12.*N./(1+h1.*(a10*N0+a11.*N))-c*N;
end

%Selection gradient
function secgrad = secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x);
    a10dot = a*a10.*(x0^(2*a)-x.^(2*a))./(x.*(x0^(2*a)+x.^(2*a)));
    h1=h(w1,w2,x);
    h1dot = -w2*h1./x;
    N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    G =g*xunder^g./(xunder^g+x.^g);
    D =d*x.^d./(xbar^d+x.^d);
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
    a21dot = a11.*(G-D+B)./x;
    a12dot = -B.*a11./x;

    secgrad=(f*(N0*a10dot+N.*a21dot-h1dot.*(N0*a10+a11.*N).^2)-(1+h1.*(a10*N0+a11.*N)).*a12dot.*N)./(1+h1.*(a10*N0+a11.*N)).^2;
end

function fyy = fyy(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x);
    a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
    h1=h(w1,w2,x);
    h1dot = -w2*h1./x;
    h1ddot = w2*(w2+1)*h1./x.^2;
    N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    a10dot = a*a10.*(x0^(2*a)-x.^(2*a))./(x.*(x0^(2*a)+x.^(2*a)));
    a10ddot = 2*a10dot.^2./a10 - a10dot./x-a^2*a10./x.^2;
    G =g*xunder^g./(xunder^g+x.^g);
    D =d*x.^d./(xbar^d+x.^d);
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a21dot = a11.*(G-D+B)./x;
    a21ddot = a21dot.^2./a11-a21dot./x+a11.*(-b^2+B.^2-G.^2.*x.^g/xunder^g-xbar^d*D.^2./x.^d)./x.^2;
    a12dot = -B.*a11./x;
    a12ddot = 2*a12dot.^2./a11 - a12dot./x-b^2*a11./x.^2;
    u = a10*N0+a11.*N;
    H = 1+h1.*u;

    fyy=(f*(H.*(a10ddot*N0+a21ddot.*N-u.^2.*h1ddot)+2*(u.^3.*h1dot.^2- ...
        (a10dot*N0+a21dot.*N).^2.*h1-2*u.*(a10dot*N0+a21dot.*N).*h1dot))-a12ddot.*N.*H.^2)./H.^3;
end

function Nbardot= Nbardot(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
G =g*xunder^g./(xunder^g+x.^g);
D =d*x.^d./(xbar^d+x.^d);
a10=A0(Ao,x0,a,x);
a10dot = a*a10.*(x0^(2*a)-x.^(2*a))./(x.*(x0^(2*a)+x.^(2*a)));
a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
a11dot = a11.*(G-D)./x;
h1=h(w1,w2,x);
h1dot = -w2.*h1./x;
N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

B=-c-c*N0*h1.*a10+(f-1)*a11;
Bdot = -c*N0*(h1dot.*a10+h1.*a10dot)+(f-1).*a11dot;
A1=c.*h1.*a11;
A1dot=c.*(h1dot.*a11+h1.*a11dot);
C=f*N0*a10;
Cdot=f*N0.*a10dot;

Nbardot = (Bdot.*N-A1dot.*N.^2+Cdot)./sqrt(B.^2 +4*C.*A1);%sqrt(B.^2 +4*C.*A1)%(2*A1.*N-B)
end

function Nbarddot= Nbarddot(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
G =g*xunder^g./(xunder^g+x.^g);
D =d*x.^d./(xbar^d+x.^d);
a10=A0(Ao,x0,a,x);
a10dot = a*a10.*(x0^(2*a)-x.^(2*a))./(x.*(x0^(2*a)+x.^(2*a)));
a10ddot = 2*a10dot.^2./a10 - a10dot./x-a^2*a10./x.^2;
a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
a11dot = a11.*(G-D)./x;
a11ddot = a11dot.^2./a11-a11dot./x-a11.*(G.^2.*x.^g/xunder^g+xbar^d*D.^2./x.^d)./x.^2;
h1=h(w1,w2,x);
h1dot = -w2.*h1./x;
h1ddot = w2*(w2+1).*h1./x.^2;
N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
Ndot= Nbardot(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

B=-c-c*N0*h1.*a10+(f-1)*a11;
Bdot = -c*N0*(h1dot.*a10+h1.*a10dot)+(f-1).*a11dot;
Bddot = -c*N0*(h1ddot.*a10+2*h1dot.*a10dot+h1.*a10ddot)+(f-1).*a11ddot;
A1=c.*h1.*a11;
A1dot=c.*(h1dot.*a11+h1.*a11dot);
A1ddot=c.*(h1ddot.*a11+2*h1dot.*a11dot+h1.*a11ddot);
C=f*N0*a10;
Cddot=f*N0.*a10ddot;

Nbarddot = (Ndot.*(2*Bdot-4*N.*A1dot-2*A1.*Ndot)+N.*(Bddot-A1ddot.*N)+Cddot)./sqrt(B.^2 +4*C.*A1);%(2*A1.*N-B);%
end

function fxx = fxx(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x);
    a10dot = a*a10.*(x0^(2*a)-x.^(2*a))./(x.*(x0^(2*a)+x.^(2*a)));
    a10ddot = 2*a10dot.^2./a10 - a10dot./x-a^2*a10./x.^2;
    h1=h(w1,w2,x);
    h1dot = -w2.*h1./x;
    h1ddot = w2*(w2+1).*h1./x.^2;
    N =Nbar(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    Ndot = Nbardot(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    Nddot = Nbarddot(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    G =g*xunder^g./(xunder^g+x.^g);
    D =d*x.^d./(xbar^d+x.^d);
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a11=Aij(A,xbar,xunder,b,g,d,p,x,x);
    a11dot = a11.*(G-D)./x;
    a11ddot = a11dot.^2./a11-a11dot./x-a11.*(G.^2.*x.^g/xunder^g+xbar^d*D.^2./x.^d)./x.^2;
    a21d = -B.*a11./x;
    a21dd =2*a21d.^2./a11-a21d./x-b^2*a11./x.^2;
    a12d = a11.*(B+G-D)./x;
    a12dd = a12d.^2./a11-a12d./x+a11.*(-b^2+B.^2-G.^2.*x.^g/xunder^g-xbar^d*D.^2./x.^d)./x.^2;
    u = a10*N0+a11.*N;
    udot = a10dot*N0+a11dot.*N+a11.*Ndot;
    uddot = a10ddot*N0+a11ddot.*N+2*a11dot.*Ndot+a11.*Nddot;
    vdot = a21d.*N+a11.*Ndot;
    vddot = a21dd.*N+2*a21d.*Ndot+a11.*Nddot;
    wdot = a12d.*N+a11.*Ndot;
    wddot = a12dd.*N+2*a12d.*Ndot+a11.*Nddot;
    H = 1+h1.*u;

    fxx=(f*vddot.*H-2*f*h1.*vdot.^2-wddot.*H.^2+2*wdot.*(h1dot.*u+h1.*udot).*H+ ...
        a11.*N.*(h1ddot.*u+2*h1dot.*udot+h1.*uddot).*H-2*a11.*N.*(h1dot.*u+h1.*udot).^2)./H.^3-c*Nddot;
end

function ddt = odefcn1(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
N1 = N(1);
N2 = N(2);
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = dN1dt(N1,N2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A); %
    ddt(2) = dN2dt(N1,N2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A); %
end

function Nfun = Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x1);
    a20=A0(Ao,x0,a,x2);
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);

    N1=N(1);
    N2=N(2);

    u=a10*N0+a11.*N1+a12.*N2;
    v=a20*N0+a21.*N1+a22.*N2;

    N1fun=N1.*((f*u-a11.*N1)./(1+h1.*u)-a21.*N2./(1+h2.*v)-c*(N1+N2));
    N2fun=N2.*((f*v-a22.*N2)./(1+h2.*v)-a12.*N1./(1+h1.*u)-c*(N1+N2));
    Nfun = [N1fun;N2fun];
end

function dN1dt = dN1dt(N1,N2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x1);
    a20=A0(Ao,x0,a,x2);
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);

    u=a10*N0+a11.*N1+a12.*N2;
    v=a20*N0+a21.*N1+a22.*N2;

    dN1dt=N1.*((f*u-a11.*N1)./(1+h1.*u)-a21.*N2./(1+h2.*v)-c*(N1+N2));
end

function dN2dt = dN2dt(N1,N2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    a10=A0(Ao,x0,a,x1);
    a20=A0(Ao,x0,a,x2);
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);

    u=a10*N0+a11.*N1+a12.*N2;
    v=a20*N0+a21.*N1+a22.*N2;

    dN2dt=N2.*((f*v-a22.*N2)./(1+h2.*v)-a12.*N1./(1+h1.*u)-c*(N1+N2));
end

%{
fyystar=fyy(xstar,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
fxxstar=fxx(xstar,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

for i = 1:size(xstar,1)
    fprintf("x^* = %.6f\n",xstar(i));
    fprintf("f_{yy} = %.6f\n",fyystar(i));
    fprintf("f_{xx} = %.6f\n",fxxstar(i));
    fprintf("ESS:%d\n",fyystar(i)<0);
    fprintf("Converging:%d\n",fxxstar(i)>fyystar(i));
end
%}
%% Colour plot to show how many singular strategies there are
%{
figure;
betinv=linspace(1e-2,1,5000); %1/beta values %try starting at a larger value
N0ran = linspace(1e-6,1000,5000); %N0 values
H = zeros(length(betinv),length(N0ran));
[height,width]=size(H);
for i = 1:height
    fprintf("i=%d\n",i);
    for j = 1:width
        b = 1/betinv(i);
        N0 = N0ran(j);
        sx = @(x) secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
        inter = find(diff(sign(sx(linspace(1e-2,10,1000)))));%finds location of sign changes in linspace
        H(i,j) = numel(inter); %numer of sign changes
    end
end
imagesc(N0ran,betinv,H);
hold on;
colormap sky
colorbar;
xlabel("Environmental Richness (N_0)",'Interpreter','tex');
ylabel("Size range (1/β)",'Interpreter','tex');
set(gca,'YDir','normal'); %makes axis look normal
%}
%% 2 parameter bifurcation diagram for singular strategy classification
%{
betinv=linspace(1e-2,1,100); %1/beta values
N0ran = linspace(1e-6,1000,100); %N0 values

[N0grid,binv] = meshgrid(N0ran,betinv);
bgrid=1./binv;
%Xstars=cell(size(N0grid));
H = nan(size(bgrid));
FYY1 = nan(size(bgrid)); FYY2 = nan(size(bgrid)); FYY3 = nan(size(bgrid));
FXX1 = nan(size(bgrid)); FXX2 = nan(size(bgrid)); FXX3 = nan(size(bgrid));
z=linspace(1e-2,10,1000);
for k = 1:numel(N0grid)
    fprintf("k=%.4f\n",100*k/numel(N0grid));
    b=bgrid(k);
    N0=N0grid(k);
    sx = @(x) secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    inter = find(diff(sign(sx(z))));
    n=numel(inter);
    H(k)=n;
    xstar = zeros(1,n);
    for i =1:n
        xstar(i)=fzero(sx,z(inter(i)));
    end
    %xstar = roots(chebfun(sx,[1e-2 10]));
    %Xstars{k}=xstar;                %finding ss for this choice of b and N0
    fyystar=fyy(xstar,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    fxxstar=fxx(xstar,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    FYY1(k)=fyystar(1);             %entering values of fyy and fxx for first ss
    FXX1(k)=fxxstar(1);
    if size(xstar,2)>1                %if there are 3 ss, enter their fyy and fxx too
        FYY2(k)=fyystar(2);
        FXX2(k)=fxxstar(2);
        FYY3(k)=fyystar(3);
        FXX3(k)=fxxstar(3);
    end
end
ESS1 = FYY1<0;
ESS2 = FYY2<0;
ESS3 = FYY3<0; 

exists1 = ~isnan(FYY1)&~isnan(FXX1);
exists2 = ~isnan(FYY2)&~isnan(FXX2);
exists3 = ~isnan(FYY3)&~isnan(FXX3);

%I checked this
CONV1 = FYY1<FXX1; %Should always be true
CONV2 = FYY2<FXX2; %Should always be false
CONV3 = FYY3<FXX3; %Should always be true (where it exists)

BRANCH1 = ~ESS1 & exists1;    %=1 when the 1st ss is a branch point (0 when 1st ss CSS)
BRANCH3 = ~ESS3 & exists3;  %=1 when 3rd ss is branch point (0 when its CSS or doesn't exist)
CSS1 = ESS1 & exists1;      %=1 when 1st ss is CSS (0 when branch)
CSS3 = ESS3 &exists3;      %=1 when 3rd ss is CSS (0 when branch or doesn't exist)
ESS2 = ESS2 &exists2;
%branching if fxx>fyy>0 (convergent but not ESS)
%CSS if convergent and ESS

figure;
hold on
% evolutionary bifurcation diagram
Z = zeros(size(N0grid));
%I checked these individually first before plotting so I'm not overwriting
%anything here
Z(BRANCH1 > 0) = 1;
Z(CSS1 > 0)    = 3;
Z(BRANCH3 > 0) = 2;
Z(ESS2 > 0)    = 5;
Z(CSS3 > 0)    = 4;

contourf(N0grid,binv,Z,'LineStyle','none','FaceAlpha',0.5)
hold on;
%contour(N0grid, binv, Z, [0.5 1.5 2.5 3.5 4.5], 'k', 'LineWidth', 2)
contour(N0grid, binv, double(BRANCH1), [0.5 0.5], 'k', 'LineWidth', 2)
contour(N0grid, binv, double(BRANCH3), [0.5 0.5], 'k', 'LineWidth', 2)
contour(N0grid, binv, double(CSS1),    [0.5 0.5], 'k', 'LineWidth', 2)
contour(N0grid, binv, double(CSS3),    [0.5 0.5], 'k', 'LineWidth', 2)
contour(N0grid, binv, double(ESS2),    [0.5 0.5], 'k', 'LineWidth', 2)
colormap([0.7 0 0.2; 0.7 0.4  0.6; 0.7 0.8 0.9; 1 0.9 0.5; 0.9 0.4 0.3]); %red BRANCH1, purple BRANCH3, blue CSS1, yellow CSS3, orange ESS2
ylabel('$\textsf{Size range}$ ($\frac{1}{\beta}$)','Interpreter','latex','FontSize',20);
xlabel('$\textsf{Common resource density}$ ($N_0$)','Interpreter','latex','FontSize',20);

%axis equal
axis square
box on
%}

%% Create plots side-by-side
%{
%I:
b1=5; N01=500;
%II
b2=1.5; N02=750;
%III:
b3=1.5; N03=10;
%IV:
b4=1.9; N04=300;
%V
b5=2.05; N05=350;

z = linspace(0,5,1000);
[X,Y] = meshgrid(z, z);

figure;
tile = tiledlayout(2,3);

% I
ax1 = nexttile;
sx = @(x) secgrad(x,a,b1,g,d,c,f,p,w1,w2,N01,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
F = fit(X,Y,a,b1,g,d,c,f,p,w1,w2,N01,x0,xbar,xunder,Ao,A);

hold on
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap(ax1,[1 1 1; 0.7 0.8 0.9]); %blue
%xlabel(ax1,'$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel(ax1,'$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
xlim(ax1,[0 1.2]);
ylim(ax1,[0 1.2]);
axis(ax1,"square");
box on
title(ax1,'I');

% II
ax2 = nexttile;
sx = @(x) secgrad(x,a,b2,g,d,c,f,p,w1,w2,N02,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
F = fit(X,Y,a,b2,g,d,c,f,p,w1,w2,N02,x0,xbar,xunder,Ao,A);

hold on
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap(ax2,[1 1 1; 0.7 0 0.2]); %red
%xlabel(ax2,'$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
%ylabel(ax2,'$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
axis(ax2,"square");
box on
title(ax2,'II');

% III
ax3 = nexttile;
sx = @(x) secgrad(x,a,b3,g,d,c,f,p,w1,w2,N03,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
F = fit(X,Y,a,b3,g,d,c,f,p,w1,w2,N03,x0,xbar,xunder,Ao,A);

hold on
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap(ax3,[1 1 1; 0.7 0.4  0.6]); %purple
xlabel(ax3,'$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
%ylabel(ax3,'$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
xlim(ax3,[0 4]);
ylim(ax3,[0 4]);
axis(ax3,"square");
box on
title(ax3,'III');

% IV
ax4 = nexttile;
sx = @(x) secgrad(x,a,b4,g,d,c,f,p,w1,w2,N04,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
F = fit(X,Y,a,b4,g,d,c,f,p,w1,w2,N04,x0,xbar,xunder,Ao,A);

hold on
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap(ax4,[1 1 1; 0.9 0.4 0.3]); %orange
xlabel(ax4,'$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel(ax4,'$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
xlim(ax4,[0 2.6]);
ylim(ax4,[0 2.6]);
axis(ax4,"square");
box on
title(ax4,'IV');

% V
ax5 = nexttile;
sx = @(x) secgrad(x,a,b5,g,d,c,f,p,w1,w2,N05,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
F = fit(X,Y,a,b5,g,d,c,f,p,w1,w2,N05,x0,xbar,xunder,Ao,A);

hold on
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap(ax5,[1 1 1; 1 0.9 0.5]); %yellow
xlabel(ax5,'$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
%ylabel(ax5,'$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
xlim(ax5,[0 2]);
ylim(ax5,[0 2]);
axis(ax5,"square");
box on
title(ax5,'V');

tile.TileSpacing = 'tight';
% Link the axes
%{
linkaxes([ax1,ax2],'y');
linkaxes([ax3,ax4],'y');
linkaxes([ax1,ax3,ax5],'x');
linkaxes([ax2,ax4],'x');


% Move plots closer together
yticklabels(ax2,{});
yticklabels(ax4,{});
xticks(ax1,[0:5]);
xticks(ax2,[0:5]);
xticks(ax3,[0:5]);
xticks(ax4,[0:5]);
xticks(ax5,[0:5]);
%}
%}

%% Population density

b=1.9;
N0=500;
tspan = [0 5000];
opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',2000, 'MaxFunEvals',10000);
C = orderedcolors("gem12");

x1list=linspace(0.4,2.6,50);
x2=4.1;
CA=C(5,:);
CB=C(3,:);
CC=C(1,:);
for i=1:50
    x1=x1list(i);
    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
    %stable eqm
    N1star = Nstar(1);
    N2star = Nstar(2);
    %unstable
    N2=Nbar(x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%N1=0
    N1=Nbar(x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%N2=0
    
    [t,n1] = ode45(@(t,n1) odefcn1(n1,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0 0.0001]);
    [~,n2] = ode45(@(t,n2) odefcn1(n2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 0]);
    [~,n3] = ode45(@(t,n3) odefcn1(n3,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 N2]);
    [~,n4] = ode45(@(t,n4) odefcn1(n4,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 0.0001]);
    
    figure;
    if fit(x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)>0
        if fit(x2,x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)>0
            C=CA;
        else
            C=CC;
        end
    else
        C=CB;
    end
    plot(n1(:,1),n1(:,2),'LineWidth',5,'Color',C);
    hold on;
    plot(n2(:,1),n2(:,2),'LineWidth',5,'Color',C);
    plot(n3(:,1),n3(:,2),'LineWidth',5,'Color',C);
    plot(n4(:,1),n4(:,2),'LineWidth',5,'Color',C);
    %equilibria
    plot(0,0,'o',MarkerSize=15,Color=C,MarkerFaceColor='w');
    plot(0,N2,'square',MarkerSize=15,Color=C,MarkerFaceColor='w');
    plot(N1,0,'square',MarkerSize=15,Color=C,MarkerFaceColor='w');
    plot(N1star,N2star,'.',MarkerSize=52,Color=C);
    
    axis square;
    ylim([0 2.5]);
    xlim([0 4]);
end
%{
%coexistence
x1=1;
x2=7;
[~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
%stable eqm
N1star = Nstar(1);%1.3203
N2star = Nstar(2);%0.4490
%unstable
N2=Nbar(x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%N1=0
N1=Nbar(x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%N2=0

[t,n1] = ode45(@(t,n1) odefcn1(n1,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0 0.0001]);
[~,n2] = ode45(@(t,n2) odefcn1(n2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 0]);
[~,n3] = ode45(@(t,n3) odefcn1(n3,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 N2]);
[~,n4] = ode45(@(t,n4) odefcn1(n4,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 0.0001]);

figure;
plot(n1(:,1),n1(:,2),'LineWidth',5,'Color',C(5,:));
hold on;
plot(n2(:,1),n2(:,2),'LineWidth',5,'Color',C(5,:));
plot(n3(:,1),n3(:,2),'LineWidth',5,'Color',C(5,:));
plot(n4(:,1),n4(:,2),'LineWidth',5,'Color',C(5,:));
%equilibria
plot(0,0,'o',MarkerSize=15,Color=C(5,:),MarkerFaceColor='w');
plot(0,N2,'o',MarkerSize=15,Color=C(5,:),MarkerFaceColor='w');
plot(N1,0,'o',MarkerSize=15,Color=C(5,:),MarkerFaceColor='w');
plot(N1star,N2star,'.',MarkerSize=50,Color=C(5,:));

set(gca,"XTickLabel",[],"YTickLabel",[]);
axis square;

%N1=0
x1=0.1;
x2=0.5;
[t,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
%stable eqm
N1star = Nstar(1);%0
N2star = Nstar(2);%3.9317
%unstable
N1=Nbar(x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
N2=0;

[t,n1] = ode45(@(t,n1) odefcn1(n1,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0 0.0001]);
[~,n2] = ode45(@(t,n2) odefcn1(n2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 0]);
[~,n3] = ode45(@(t,n3) odefcn1(n3,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 N2+0.0001]);
[~,n4] = ode45(@(t,n4) odefcn1(n4,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1star+0.0001 N2star]);

figure;
plot(n1(:,1),n1(:,2),'LineWidth',5,'Color',C(1,:));
hold on;
plot(n2(:,1),n2(:,2),'LineWidth',5,'Color',C(1,:));
plot(n3(:,1),n3(:,2),'LineWidth',5,'Color',C(1,:));
plot(n4(:,1),n4(:,2),'LineWidth',5,'Color',C(1,:));
%equilibria
plot(0,0,'o',MarkerSize=15,Color=C(1,:),MarkerFaceColor='w');
plot(N1star,N2star,'.',MarkerSize=50,Color=C(1,:));
plot(N1,N2,'o',MarkerSize=15,Color=C(1,:),MarkerFaceColor='w');

set(gca,"XTickLabel",[],"YTickLabel",[]);
axis square;

%N2=0
x1=0.1;
x2=5;
[t,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
%stable eqm
N1star = Nstar(1);%3.3531
N2star = Nstar(2);%0
%unstable
N1=0;
N2=Nbar(x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

[t,n1] = ode45(@(t,n1) odefcn1(n1,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0 0.0001]);
[~,n2] = ode45(@(t,n2) odefcn1(n2,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.0001 0]);
[~,n3] = ode45(@(t,n3) odefcn1(n3,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1+0.0001 N2]);
[~,n4] = ode45(@(t,n4) odefcn1(n4,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1star N2star+0.0001]);

figure;
plot(n1(:,1),n1(:,2),'LineWidth',5,'Color',C(3,:));
hold on;
plot(n2(:,1),n2(:,2),'LineWidth',5,'Color',C(3,:));
plot(n3(:,1),n3(:,2),'LineWidth',5,'Color',C(3,:));
plot(n4(:,1),n4(:,2),'LineWidth',5,'Color',C(3,:));
%equilibria
plot(0,0,'o',MarkerSize=15,Color=C(3,:),MarkerFaceColor='w');
plot(N1,N2,'o',MarkerSize=15,Color=C(3,:),MarkerFaceColor='w');
plot(N1star,N2star,'.',MarkerSize=50,Color=C(3,:));

set(gca,"XTickLabel",[],"YTickLabel",[]);
axis square;
%}