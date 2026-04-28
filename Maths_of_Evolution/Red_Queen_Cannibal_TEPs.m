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

sx = @(x) secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);

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

%% Coexistence
xmax=11.5;
xmin=2.5e-2;
b=1.5;
N0=750;

%Fitness and selection gradients

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

function fit12=fit12(x1,x2,y,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2)
    a10=A0(Ao,x0,a,x1);
    a20=A0(Ao,x0,a,x2);
    a30=A0(Ao,x0,a,y);
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    a13=Aij(A,xbar,xunder,b,g,d,p,x1,y);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    a23=Aij(A,xbar,xunder,b,g,d,p,x2,y);
    a31=Aij(A,xbar,xunder,b,g,d,p,y,x1);
    a32=Aij(A,xbar,xunder,b,g,d,p,y,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);
    h3=h(w1,w2,y);

    H1 = 1+h1.*(a10*N0+a11.*N1+a12.*N2);
    H2 = 1+h2.*(a20*N0+a21.*N1+a22.*N2);
    H3 = 1+h3.*(a30*N0+a31.*N1+a32.*N2);
    
    fit12 = f*(a30*N0+a31.*N1+a32.*N2)./H3-a13.*N1./H1-a23.*N2./H2-c*(N1+N2);
end

function secgrad1 = secgrad1(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2)
    
    a10=A0(Ao,x0,a,x1);
    a10dot = a*a10.*(x0^(2*a)-x1.^(2*a))./(x1.*(x0^(2*a)+x1.^(2*a)));
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    a20=A0(Ao,x0,a,x2);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);
    h1dot = -w2*h1./x1;

    G1 =g*xunder^g./(xunder^g+x1.^g);
    D1 =d*x1.^d./(xbar^d+x1.^d);
    B1 =b*(x2.^(2*b)-(p*x1).^(2*b))/((p*x1).^(2*b)+x2.^(2*b));
    B2 =b*(x1.^(2*b)-(p*x2).^(2*b))/((p*x2).^(2*b)+x1.^(2*b));
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a31dot = a11.*(B+G1-D1)./x1;
    a32dot = a12.*(B1+G1-D1)./x1;
    a13dot = -B.*a11./x1;
    a23dot = -a21.*B2./x1;

    A1 = a10*N0+a11.*N1+a12.*N2;
    H1 = 1+h1.*A1;
    H2 = 1+h2.*(a20*N0+a21.*N1+a22.*N2);

    secgrad1=(f*(a10dot*N0+a31dot.*N1+a32dot.*N2-h1dot.*A1.^2)-a13dot.*N1.*H1)./H1.^2 -a23dot.*N2./H2;
end

function secgrad2 = secgrad2(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2)
    
    a10=A0(Ao,x0,a,x1);
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    a20=A0(Ao,x0,a,x2);
    a20dot = a*a20.*(x0^(2*a)-x2.^(2*a))./(x2.*(x0^(2*a)+x2.^(2*a)));
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);
    h2d = -w2*h2./x2;

    B1 =b*(x2.^(2*b)-(p*x1).^(2*b))/((p*x1).^(2*b)+x2.^(2*b));
    B2 =b*(x1.^(2*b)-(p*x2).^(2*b))/((p*x2).^(2*b)+x1.^(2*b));
    G2 =g*xunder^g./(xunder^g+x2.^g);
    D2 =d*x2.^d./(xbar^d+x2.^d);
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a31d = a21.*(B2+G2-D2)./x2;
    a32d = a22.*(B+G2-D2)./x2;
    a13d = -B1.*a12./x2;
    a23d = -a22.*B./x2;

    A2 = a20*N0+a21.*N1+a22.*N2;
    H1 = 1+h1.*(a10*N0+a11.*N1+a12.*N2);
    H2 = 1+h2.*A2;

    secgrad2=(f*(a20dot*N0+a31d.*N1+a32d.*N2-h2d.*A2.^2)-a23d.*N2.*H2)./H2.^2 -a13d.*N1./H1;
end


z=logspace(log10(xmin),log10(xmax),200);%for log scale
%z=linspace(1e-6,xmax,100);             %for regular
z1=linspace(0.015,0.045,100);
z2=linspace(1.416,1.446,100);
[X1,X2] = meshgrid(z1, z2);
N1grid = zeros(size(X1));
N2grid = zeros(size(X2));

sx1= @(x1,x2,N1,N2) secgrad1(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2);
sx2= @(x1,x2,N1,N2) secgrad2(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2);

S1 = zeros(size(X1));
S2 = zeros(size(X2));

%using the limit of dN1/dt dN2/dt to find N to find the selection gradient
%{
tspan = [0 4000];

opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',2000, 'MaxFunEvals',10000);
%For upper half plane only
%{
k=0;
for i = 1:size(X2,1)
    for j = 1:i
        fprintf("p=%.4f\n",200*k/numel(S1));
        k=k+1;
        x1=X1(i,j);
        x2=X2(i,j);

        nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
        [t,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
        Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
    
        N1= Nstar(1);
        N2= Nstar(2);
        
        N1grid(i,j)=N1;
        N2grid(i,j)=N2;
        S1(i,j)=sx1(x1,x2,N1,N2);
        S2(i,j)=sx2(x1,x2,N1,N2);
    end
end
%}
%For full plot
%{
for k=1:numel(S1)
    fprintf("k=%.4f\n",100*k/numel(S1));
    x1=X1(k);
    x2=X2(k);
    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    [t,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);

    N1= Nstar(1);
    N2= Nstar(2);
    N1grid(k)=N1;
    N2grid(k)=N2;
    if norm(nfun(Nstar))>1e-6
        norm(nfun(Nstar))
        fprintf('x1=%.2f',x1)
        fprintf('x2=%.2f\n',x2)
        fprintf('N1=%.2f',N1)
        fprintf('N2=%.2f\n',N2)
        fprintf('N1=%.2f',n(end,1))
        fprintf('N2=%.2f',n(end,2))
        norm(nfun([n(end,1) n(end,2)]))
    end

    S1(k)=sx1(x1,x2,N1,N2);
    S2(k)=sx2(x1,x2,N1,N2);
end
%}

%% TEP
%
z=logspace(log10(xmin),log10(xmax),5000);%for log scale
%z=linspace(0,xmax,2000);               %for regular
[X,Y] = meshgrid(z1, z2);

F1 = fit(X,Y,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
F2 = fit(Y,X,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
sx = @(x) secgrad(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
xstar=roots(chebfun(sx,[1e-6 5]));

F = zeros(size(X));
F(F1>0)=1;   %fx(y)>0
F(F2>0)=2;   %fx(y)<0
F(F1>0&F2>0)=3;   %mutual invasibility
F(X>Y)=0;       %shows upper half only

%don't care about direction outside of coexistence region
%S1(fit(X1,X2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)<0|fit(X2,X1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)<0)=0;
%S2(fit(X1,X2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)<0|fit(X2,X1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)<0)=0;

%{
figure;
hold on
% PIP1
contourf(X,Y,F1>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F1,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap([1 1 1; 0.7 0.8 0.9]); %blue for f>0 white for f<0
xlabel('$\textsf{Resident trait value}$ ($x_1$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Mutant trait value}$ ($x_2$)','Interpreter','latex','FontSize',20);

axis equal
box on

figure;
hold on
% PIP2
contourf(X,Y,F2>0,1,'LineStyle','none','DisplayName','');
contour(X,Y,F2,[0 0],'k','LineWidth',2,'DisplayName',''); %black line for f=0
colormap([1 1 1; 1 0.9 0.5]); %yellow for f>0 white for f<0
ylabel('$\textsf{Resident trait value}$ ($x_2$)','Interpreter','latex','FontSize',20);
xlabel('$\textsf{Mutant trait value}$ ($x_1$)','Interpreter','latex','FontSize',20);

axis equal
box on
%}

figure;
hold on
% TEP
%contourf(log10(X), log10(Y), F, 'LineStyle','none');                   %for log scale
%contour(log10(X), log10(Y), F, [0.5 1.5 2.5], 'k', 'LineWidth', 2);    %for log scale
contourf(X,Y,F,'LineStyle','none');                                   %for regular
contour(X,Y,F,[0.5 1.5 2.5],'k','LineWidth',2); %black line for f=0   %for regular
colormap([1 1 1; 0.8 0.9 0.9; 1 0.9 0.7; 0.6 0.9 0.3]); %blue for fx(y)>0, yellow for fy(x)>0, green for both
xlabel('$\textsf{Resident trait value}$ ($x_1$)','Interpreter','latex','FontSize',20);
ylabel('$\textsf{Resident trait value}$ ($x_2$)','Interpreter','latex','FontSize',20);
axis equal
%
step = 8;%stepsize of z/25
X11= X1(1:step:end,1:step:end);
X22= X2(1:step:end,1:step:end);
S11=S1(1:step:end,1:step:end);
S22=S2(1:step:end,1:step:end);


%r = sqrt((S11./X11).^2+(S22./X22).^2);  %for log scale
r = sqrt(S11.^2+S22.^2);                 %for regular
r(r==0) = 1;

L = 0.08;
%quiver(log10(X11),log10(X22),L*(S11./X11./r),L*(S22./X22./r),0,"LineWidth",1,"Color","k",'Alignment','center');  %for log space
quiver(X11,X22,S11./r,S22./r,0.6,"LineWidth",1,"Color","k","Alignment","center");                                   %for regular
legend('','');

%contour(log10(X1),log10(X2),S1,[0 0],'r:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_1)=0');    %for log space
%contour(log10(X1),log10(X2),S2,[0 0],'b:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_2)=0');    %for log space
contour(X1,X2,S1,[0 0],'r:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_1)=0');                 %for regular
contour(X1,X2,S2,[0 0],'b:','LineWidth',2,'DisplayName','s_{x_1,x_2}(x_2)=0');                 %for regular

%xline(log10(xstar(end)), 'LineStyle','--', 'LineWidth',1.5, 'Color',C(3,:), 'DisplayName','B');       %for log space
%xline(log10(xstar(2)), 'LineStyle','--', 'LineWidth',1.5, 'Color',C(2,:), 'DisplayName','D');         %for log space
%xline(log10(xstar(1)), 'LineStyle','--', 'LineWidth',1.5, 'Color',C(1,:), 'DisplayName','CSS');       %for log space
xline(xstar(end),'LineStyle','--','LineWidth',1.5,'Color',C(3,:),'DisplayName',['B=',num2str(xstar(end),3)])                 %for regular
xline(xstar(2),'LineStyle','--','LineWidth',1.5,'Color',C(2,:),'DisplayName',['D=',num2str(xstar(2),3)])             %for regular
xline(xstar(1),'LineStyle','--','LineWidth',1.5,'Color',C(1,:),'DisplayName',['CSS=',num2str(xstar(1),3)])             %for regular

%set(gca,'XTick',[-1,0,1], 'YTick', [-1,0,1]);               %for log space
%set(gca,'XTickLabel',[0.1,1,10], 'YTickLabel', [0.1,1,10]); %for log space

legend;
axis equal
box on
%}
%% 2nd deriv

function fit1yy=fit1yy(x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    tspan = [0 4000];
    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
    opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',2000, 'MaxFunEvals',10000);
    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);

    N1= Nstar(1);
    N2= Nstar(2);
    
    a10=A0(Ao,x0,a,x1);
    a10dot = a*a10.*(x0^(2*a)-x1.^(2*a))./(x1.*(x0^(2*a)+x1.^(2*a)));
    a10ddot = 2*a10dot.^2./a10 - a10dot./x1-a^2*a10./x1.^2;
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    a20=A0(Ao,x0,a,x2);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);
    h1dot = -w2*h1./x1;
    h1ddot = w2*(w2+1)*h1./x1.^2;
    

    G1 =g*xunder^g./(xunder^g+x1.^g);
    D1 =d*x1.^d./(xbar^d+x1.^d);
    B1 =b*(x2.^(2*b)-(p*x1).^(2*b))/((p*x1).^(2*b)+x2.^(2*b));
    B2 =b*(x1.^(2*b)-(p*x2).^(2*b))/((p*x2).^(2*b)+x1.^(2*b));
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a31dot = a11.*(B+G1-D1)./x1;
    a31ddot = a31dot.^2./a11-a31dot./x1+a11.*(-b^2+B.^2-G1.^2.*x1.^g./xunder^g-xbar^d*D1.^2./x1.^d)./x1.^2;
    a32dot = a12.*(B1+G1-D1)./x1;
    a32ddot = a32dot.^2./a12-a32dot./x1+a12.*(-b^2+B1.^2-G1.^2.*x1.^g./xunder^g-xbar^d*D1.^2./x1.^d)./x1.^2;
    a13dot = -B.*a11./x1;
    a13ddot = 2*a13dot.^2./a11-a13dot./x1-b^2*a11./x1.^2;
    a23dot = -a21.*B2./x1;
    a23ddot = 2*a23dot.^2./a21-a23dot./x1-b^2*a21./x1.^2;
    
    A1=a10*N0+a11*N1+a12*N2;
    A3dot = a10dot*N0+a31dot*N1+a32dot*N2;
    A3ddot = a10ddot*N0+a31ddot*N1+a32ddot*N2;
    H1 = 1+h1.*A1;
    H2 = 1+h2.*(a20*N0+a21.*N1+a22.*N2);
    
    fit1yy = (f*(H1.*(A3ddot-A1.^2.*h1ddot)+2*(A1.^3.*h1dot.^2-A3dot.^2.*h1-2*A1.*A3dot.*h1dot))-N1*a13ddot.*H1.^2)./H1.^3-N2*a23ddot./H2;
end

function fit2yy=fit2yy(x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    tspan = [0 4000];
    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [0.2 0.2]);
    opts = optimset('Display','off', 'TolX',1e-12, 'TolFun',1e-12, 'MaxIter',2000, 'MaxFunEvals',10000);
    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);

    N1= Nstar(1);
    N2= Nstar(2);
    
    a20=A0(Ao,x0,a,x2);
    a20d = a*a20.*(x0^(2*a)-x2.^(2*a))./(x2.*(x0^(2*a)+x2.^(2*a)));
    a20dd = 2*a20d.^2./a20 - a20d./x2-a^2*a20./x2.^2;
    a11=Aij(A,xbar,xunder,b,g,d,p,x1,x1);
    a12=Aij(A,xbar,xunder,b,g,d,p,x1,x2);
    a10=A0(Ao,x0,a,x1);
    a21=Aij(A,xbar,xunder,b,g,d,p,x2,x1);
    a22=Aij(A,xbar,xunder,b,g,d,p,x2,x2);
    h1=h(w1,w2,x1);
    h2=h(w1,w2,x2);
    h2d = -w2*h2./x2;
    h2dd = w2*(w2+1)*h2./x2.^2;
    

    G2 =g*xunder^g./(xunder^g+x2.^g);
    D2 =d*x2.^d./(xbar^d+x2.^d);
    B1 =b*(x2.^(2*b)-(p*x1).^(2*b))/((p*x1).^(2*b)+x2.^(2*b));
    B2 =b*(x1.^(2*b)-(p*x2).^(2*b))/((p*x2).^(2*b)+x1.^(2*b));
    B =b*(1-p^(2*b))/(p^(2*b)+1);
    a31d = a21.*(B2+G2-D2)./x2;
    a31dd = a31d.^2./a21-a31d./x2+a21.*(-b^2+B2.^2-G2.^2.*x2.^g./xunder^g-xbar^d*D2.^2./x2.^d)./x2.^2;
    a32d = a22.*(B+G2-D2)./x2;
    a32dd = a32d.^2./a22-a32d./x2+a22.*(-b^2+B.^2-G2.^2.*x2.^g./xunder^g-xbar^d*D2.^2./x2.^d)./x2.^2;
    a13d = -B1.*a12./x2;
    a13dd = 2*a13d.^2./a12-a13d./x2-b^2*a12./x2.^2;
    a23d = -a22.*B./x2;
    a23dd = 2*a23d.^2./a22-a23d./x2-b^2*a22./x2.^2;
    
    A2=a20*N0+a21*N1+a22*N2;
    A3d = a20d*N0+a31d*N1+a32d*N2;
    A3dd = a20dd*N0+a31dd*N1+a32dd*N2;
    H1 = 1+h1.*(a10*N0+a11.*N1+a12.*N2);
    H2 = 1+h2.*A2;
    
    fit2yy = (f*(H2.*(A3dd-A2.^2.*h2dd)+2*(A2.^3.*h2d.^2-A3d.^2.*h2-2*A2.*A3d.*h2d))-N2*a23dd.*H2.^2)./H2.^3-N1*a13dd./H1;
end

function S = sx12_fun(X,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,opts,tspan)
    x1 = X(1);
    x2 = X(2);
    
    nfun = @(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A), ...
                  tspan, [0.2 0.2]);
    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);

    N1 = Nstar(1);
    N2 = Nstar(2);

    S = [ ...
        secgrad1(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2); ...
        secgrad2(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2) ...
    ];
end
%% Tests location and type of singular strategies
%{
%guesses
x1=0.02977624985;
x2=1.087314801;

%solving selection gradients=0 near these guesses
sx12 = @(X) sx12_fun(X,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,opts,tspan);
S=fsolve(sx12,[x1,x2],opts);
S(1)
S(2)

%Printing second derivatives at the singular strategy
fit1yy(S(1),S(2),a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
fit2yy(S(1),S(2),a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)
%}

%% Simulation:
%
%Stochastic simulation of trait evolution
tspan = [0 5000];
opts = optimset('Display','off');
numberofrealisations = 1; % number of sample paths
T = 4000; % end time for simulation
xstar=roots(chebfun(sx,[1e-6 5]));
xinitial = xstar(end); %branch point
N1initial = Nbar(xinitial,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)/2;
N2initial = N1initial;
% matrix to hold trait value at each timestep, for each realisation
xplot = zeros(T+1,numberofrealisations);
xbplot = zeros(T+1,numberofrealisations);
timeplot = zeros(T+1,numberofrealisations);
m=0.01;
for i = 1:numberofrealisations %Run numberofrealisations simulations
    i
    x1=xinitial;
    x2=xinitial; %resident of the branch
    N1= N1initial;
    N2= N2initial;
    time=0;
    j=1;
    xplot(j,i)=x1;       %To keep track of trait values for simulation i at "time" j
    xbplot(j,i)=x2;     %To keep track of branched trait values for simulation i at "time" j
    timeplot(j,i)=0;    %To keep track of time that trait x is obtained for simulation i
    branch = false;     %means the trait hasn't branched yet

    while time<T        %goes until specified time T
        %only runs for the trait value sufficiently far from x*
        if branch==false
            y = x1+2*m*rand-m; %mutant value is resident +-random number <= m
            if fit(x1,y,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) > 0
                if fit(y,x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)>0
                    x2=max(x1,y);%update branch value if mutally invasive s.t. x2>x1
                    x1=min(x1,y); 
                    branch=true;
                    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%check pop dens of new system
                    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 N2]);
                    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
                    N1 = Nstar(1);
                    N2 = Nstar(2);
                else
                    x1=y; %update trait value if y has positive invasion fitness
                    x2=x1;
                    N1=Nbar(x1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A)/2;
                    N2=N1;
                end
            end
        %This runs once we have got sufficiently close to "branching" point
        else
            if rand<=N1/(N1+N2) %pick which subpopulation is mutating scaled by relative abundance
                y1 = x1+2*m*rand-m;
                if fit12(x1,x2,y1,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2)>0
                    x1=y1;
                    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%check pop dens of new system
                    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 N2]);
                    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
                    N1 = Nstar(1);
                    N2 = Nstar(2);
                    if N1<1e-20%if branch 1 dies
                        N1
                        x1=x2;
                        branch=false;%we return to single resident system, with trait value x2
                        N1=N2/2;
                        N2=N1;
                    elseif N2<1e-20%if branch 2 dies
                        N2
                        x2=x1;
                        branch=false;%we return to single resident system, with trait value x1
                        N2=N1/2;
                        N1=N2;
                    end
                end
            else %if trait 2 mutates
                y2 = x2+2*m*rand-m;
                if fit12(x1,x2,y2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2)>0
                    x2=y2;
                    nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);%check pop dens of new system
                    [~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , tspan, [N1 N2]);
                    Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
                    N1 = Nstar(1);
                    N2 = Nstar(2);
                    if N1<1e-20%if branch 1 dies
                        N1
                        x1=x2;
                        branch=false;%we return to single resident system, with trait value x2
                        N1=N2/2;
                        N2=N1;
                    elseif N2<1e-20%if branch 2 dies
                        N2
                        x2=x1;
                        branch=false;%we return to single resident system, with trait value x1
                        N2=N1/2;
                        N1=N2;
                    end
                end
                if N1<0
                    N1=0;
                end
                if N2<0
                    N2=0;
                end
                %only need to check N if there's been a trait substitution
            end
        end
        time=time+1%update timestep
        %update matrices
        j=j+1;
        xplot(j,i)=x1;
        xbplot(j,i)=x2;
        timeplot(j,i)=time;
    end
end

%figure;
C = orderedcolors("gem");
hold on;
for i=1:numberofrealisations %shows trait trajectory in trait space
    k=stairs(min(xplot(:,i),xbplot(:,i)),max(xplot(:,i),xbplot(:,i)),'Color',C(mod(i,7)+1,:));
    set(k,'Linewidth',2);
end
axis square;
figure;
hold on;
for i=1:numberofrealisations  %shows branching against time
    stairs(timeplot(:,i),xplot(:,i),'Color',C(mod(i,7)+1,:),'Linewidth',2);
    stairs(timeplot(:,i),xbplot(:,i),'Color',C(mod(i,7)+1,:),'Linewidth',2);
end
%}
%Alternatively the determinsitic simulation of trait evolution using the
%canonical equation
%{
function ddt = odefcn2(X,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,opts)
x1=X(1);
x2=X(2);
nfun=@(N) Nfun(N,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A);
[~,n] = ode45(@(t,n) odefcn1(n,x1,x2,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A) , [0 1000], [0.2 0.2]);
Nstar = fsolve(nfun,[n(end,1) n(end,2)],opts);
N1 = Nstar(1);
N2 = Nstar(2);
if N1<1e-200||N2<1e-200
    fprintf("N1=%f\n",N1);
    fprintf("N2=%f\n",N2);
    fprintf("x1=%.2f\n",x1);
    fprintf("x2=%f\n",x2);
    ddt=zeros(2,1);
    return
end
ddt = zeros(2,1); % Initialize the derivative vector
    ddt(1) = N1*secgrad1(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2); %
    ddt(2) = N2*secgrad2(x1,x2,a,b,g,d,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,N1,N2); %
end

[t,x]=ode45(@(t,x) odefcn2(x,a,b,g,d,c,f,p,w1,w2,N0,x0,xbar,xunder,Ao,A,opts) , [0 500], [xinitial xinitial+0.1]);
plot(x(:,1),x(:,2))
%}