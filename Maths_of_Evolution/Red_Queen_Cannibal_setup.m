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
%% Attack rate

x0 = 0.1;%3;%
Ao = 1;
a = 2;
a1=5;

function A0= A0(A,x0,a,x)
A0 = 2*A*x.^a./((x.^2/x0).^a+(x0)^a);
end
%
xl=0.5;
x=linspace(0,xl,500);
figure;
plot(x, A0(Ao,x0,a,x));
hold on;
yline(Ao,'--',"Label","A_{i0}",Interpreter="tex",LineWidth=2,FontSize=25);
xline(x0,'--',Interpreter="tex",LineWidth=2,FontSize=25);
plot(x, A0(Ao,x0,a1,x));
legend(sprintf('α=%.2f',a),'','',sprintf('α=%.2f',a1),FontSize=25);
xtickVals = unique([0 : xl/5 : xl, x0]);
xtickLabs = compose('%.3g',xtickVals);
xtickLabs(ismembertol(xtickVals,x0)) = {'x_0'};
set(gca,'xtick', xtickVals, 'xticklabel', xtickLabs, 'xlim', [min(xtickVals),max(xtickVals)],FontSize=25)
ytickVals = unique([0 : 0.2 : 1, Ao]);
ytickLabs = compose('%.3g',ytickVals);
ytickLabs(ismembertol(ytickVals,Ao)) = {'A_{i0}'};
set(gca,'ytick', ytickVals, 'yticklabel', ytickLabs, 'ylim', [min(ytickVals),max(ytickVals)],'TickLabelInterpreter','tex',FontSize=25)
xlabel("Trait value (x_i)",'Interpreter','tex',FontSize=25);
ylabel("Rate of common resource consumption (a_{i0})",'Interpreter','tex',FontSize=25);

figure;
%}
function Aij= Aij(A,xbar,xunder,b,g,d,p,x,y)
Aij = A*2./((p*x./y).^(b)+(y./(p*x)).^(b)).*(x.^g./(xunder^g+x.^g)).*(xbar^d./(xbar^d+x.^d));
end
%{
limit = 30;
[X,Y] = meshgrid(0:0.05:limit,0:0.05:limit);
xbar = 13;
xunder = 3;
A=5;
b=5;
g=2;
d=2;
p=0.5;
Z=Aij(A,xbar,xunder,b,g,d,p,X,Y);
surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.8);
hold on;
%grid
[X,Y] = meshgrid(0:0.3:limit,0:0.3:limit);
Z=Aij(A,xbar,xunder,b,g,d,p,X,Y);
surf(X,Y,Z,'EdgeAlpha',0.5,'FaceAlpha',0);
hold on;

%max value of Aij
x=linspace(0,limit);
Zmax=Aij(A,xbar,xunder,b,g,d,p,x,p.*x);
plot3(x,p.*x,Zmax,'k','LineWidth',2);
Zsame=Aij(A,xbar,xunder,b,g,d,p,x,x);
plot3(x,x,Zsame,'b--','LineWidth',2);

patch([xunder xunder xunder xunder],[0 0 limit limit], [0 max(Z(:)) max(Z(:)) 0],'red','FaceAlpha',0.3);
patch([xbar xbar xbar xbar],[0 0 limit limit], [0 max(Z(:)) max(Z(:)) 0],[1 0.5 0],'FaceAlpha',0.3);

xlabel('x_i',Interpreter='tex');
ylabel('x_j',Interpreter='tex');
zlabel('a_{ij}',Interpreter='tex')
set(gca,'YDir','reverse');

xtickVals = unique([0 : 5 : limit, xunder,xbar]);
xtickLabs = compose('%.3g',xtickVals);
xtickLabs(ismembertol(xtickVals,xunder)) = {'\underline{x}'};
xtickLabs(ismembertol(xtickVals,xbar)) = {'\bar{x}'};
set(gca,'xtick', xtickVals, 'xticklabel', xtickLabs,'TickLabelInterpreter','latex')
%}
%% handling time

function hij= h(w1,w2,x)
hij = w1*x.^(-w2);
end
%{
x=linspace(0,30,5000);

w11=0.1;
w12=0.9;
w21=0.2;
w22=0.9;

figure;
plot(x,h(w11,w21,x),"DisplayName",['w_1=',num2str(w11),', w_2=',num2str(w21)]);
hold on;
plot(x,h(w12,w21,x),"DisplayName",['\color{red}w_1=',num2str(w12),'\color{black}, w_2=',num2str(w21)]);
plot(x,h(w11,w22,x),"DisplayName",['w_1=',num2str(w11),',\color{red} w_2=',num2str(w22)]);
%plot(x,h(w12,w22,x),"DisplayName",['\color{red}w_1=',num2str(w12),'\color{red}, w_2=',num2str(w22)]);
legend('Interpreter','tex');

xlabel('Trait Value (x_i)','Interpreter','tex');
ylabel('Handling time (h_i)','Interpreter','tex');
ylim([0 2])
%}
%% Selection gradient

x0 = 0.1;
xbar = 5;
xunder = 0.5;
Ao = 1;
A=10;
a = 2;
b=1.5;%2.5;%
g=8;
d=2;
c=1;
f=0.6;
p=0.2;
w1=0.1;
w2=0.25;
N0=10;%500;%200;%

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

z = linspace(0,5,500);
sx = @(x) secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));
%

figure;
z = linspace(0,3,5000);
S=secgrad(z,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
plot(z,S);
grid on;
hold on;
xlabel("Trait value (x)");
ylabel("Selection gradient (s(x))");
legend('');
%evolutionary singular strategies
xstar=roots(chebfun(sx,[1e-6 5]));
for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end

legend;

%% Population density by trait value
figure;
N=Nbar(z,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
plot(z,N);
grid on;
hold on;
legend('');
plot(xstar(1),Nbar(xstar(1),a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A),'.',MarkerSize=25,Color=C(1,:));
plot(xstar(2),Nbar(xstar(2),a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A),'o',MarkerSize=10,Color=C(2,:));
plot(xstar(3),Nbar(xstar(3),a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A),'.',MarkerSize=25,Color=C(3,:));
legend('');
for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',2,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
xlabel("Body Size (x)");
ylabel("Equilibrium population density (N_x*)",'Interpreter','tex');
%}
%% PIP
%
[X,Y] = meshgrid(z, z);
F = fit(X,Y,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%fit(X,Y);%

figure;
hold on
% PIP
contourf(X,Y,F>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap([1 1 1; 0.6 0.9 0.3]); %green for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Mutant trait value}$ ($y$)','Interpreter','latex','FontSize',20);
legend('','');

for i = 1:size(xstar,1)
    xline(xstar(i),'LineStyle','--','LineWidth',1.5,'DisplayName',['x*=',num2str(xstar(i),3)],'Color',C(i,:)) %evolutionary singular strategy
end
legend;
axis equal
box on
%}
%% Second derivatives

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
