clear all
format long
set(0,'defaulttextinterpreter','Latex')

global I N St Pes kappas beta Zc M deta eta nu P1 gamma P0 t1min s1 dsdt F1 qwdot0 chi
scale =  4;
modelII = 1;
Xj = [0 0.85 1.7];
Xx = 0;
L = 0.705;
alpha = 2.8e-5; %?????????????
Tmelt = 1773;
Tw = 273+30;
Tscale = Tmelt-Tw;
%DEv = [0.025 0.05 0.075];
kaj = [0.09];
%Vwj = [ 0.2 0.22 0.24 0.26];
%Vwj = [ 2.5 5 7.5 10];
Vwj = 8;
tic
for j = 1:length(Xj)
Xx = Xj(j);
%Xx = 1.7;
M = L*Xx/(100*alpha*Tscale);
%M = Mj(j);
%M = 0;
cps = 620;
cpl = 700;
cpw = 4185.5;
E = 200e9;
%alpha =alpha;
ks = 40;
kl = 30;
kM = 380;
kw = 0.5142;
ka = 0.075;
%ka = kaj(j);
Tmelt = 1773;
Vcast = 0.025*1;
DeltaHf = 273000;
rhos = 7800;
rhow = 997;
rhoww = rhow;
nu = 0.3;
W = 0.0914;
WM = W+0.012;
%DE = DEv(j);
DE = 0.05;
%Vw = Vwj(j);
Vw = Vwj;
muw = 0.8e-3;



Tcast = Tmelt;

Re = rhow*Vw*DE/muw;
Pr = cpw*muw/kw;
Nu = 0.023*Re^(0.8)*Pr^(0.4);

hC = kw/DE*Nu;

g = 9.81;

p0 = 1.01e5;

Tscale = Tmelt-Tw;

qwdot0 = M;
Pes = rhos*cps*Vcast*L/ks*(W/L)^2;
Pel = rhos*cpl*Vcast*L/kl*(W/L)^2;
St = cps*Tscale/DeltaHf;
delta = alpha*Tscale;
kappas = kM/ks/log(WM/W);
kappal = kM/kl/log(WM/W);
kappaa = kM/ka/log(WM/W);
kappaM = hC*WM/kM;
P0 = p0/E/alpha/Tscale;
P1 = rhos*g*L/E/alpha/Tscale;
beta = kM/log(WM/W)/hC/WM;
chi = alpha*Tscale/log(WM/W);
gamma = ka/kM*log(WM/W)/alpha/Tscale;

thetaC = (Tmelt-Tcast)/Tscale;

Beta = kappal*(thetaC-1)/(1+beta);

%Zmelt = Pes*pi*thetaC^2/4/Beta^2
%Zgap = Zmelt+ (P0+P1*Zmelt)*(kappas^2/Pes/(1+beta)^2*(St)-nu*P1-qwdot0)^-1
Zgap = 2*P1*(St^2*kappas^3/Pes^2/(1+beta)^3-St*kappas*(M+P1*(nu-1))/Pes/((1+beta)))^(-1)
%Zgap = 5e-3;
%Zgap = Zmelt+(P0+P1*Zmelt)*(kappas^2*(St+2*nu)/(1-nu)/Pes/(1+beta)^2-nu*P1-qwdot0)^-1
%Zgap-Zmelt;


etamin = 0;
etamax = 1;

I = 50*scale;

deta = (etamax-etamin)/I;
eta = etamin:deta:etamax;

N = 2500*scale;
NN = 2500*scale;

t1min = 1e-13;
Zgapguess = Zgap*2;
Zgapguess = Zgap;
Zgap = fzero('cylinder_model_zgap_fn',Zgapguess)
%Zgap = 0.1;
dt = (1-t1min)/N;
t = t1min:dt:1;

cfl = dt/deta;
be = dt/deta^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialize matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1old = zeros(I + 1,1);
F1new = zeros(I + 1,1);
F1 = zeros(I+1,N+1);

dsdt(1) = -(St*Zgap/Pes)*kappas/(1+beta);
s1(1) = 1-1e-13;

for n = 1:N+1
    %th0(n) = 1 - phi*exp(-t1(n)/Zc);
    rhow1(n) = M*t(n)*Zgap;
end 


for i = 1:I+1
    F1(i,1) = kappas/(1+beta)*(1-eta(i)); %initial condition
    F1old(i) = F1(i,1);
end


rhow1;
for n = 1:N
    
    Tri_diagL = sparse(I + 1, I + 1);     
    Tri_diagR = sparse(I + 1, I + 1); 
      
    % discretise bc at eta = 0
    Tri_diagL(1,1) = 1;  Tri_diagL(1,2) = -4 - 2*kappas*deta*(1-s1(n))/(1+beta+chi*rhow1(n)); Tri_diagL(1,3) = 3;
    Tri_diagR(1,1) = 0; Tri_diagR(1,2) = 0; 
    
    % create tridiagonal part 
    for m = 2:I
        etam = (1/2)*(eta(m) + eta(m-1));
        etap = (1/2)*(eta(m) + eta(m+1));
        Tri_diagL(m,m-1) = (be*Zgap/Pes)*(1 - (1-s1(n))*etam) + (cfl/2)*eta(m)*dsdt(n)*(1 - (1-s1(n))*eta(m))*(1-s1(n)); 
        Tri_diagL(m,m) = - ( (be*Zgap/Pes)*(1 - (1-s1(n))*etap) + (be*Zgap/Pes)*(1 - (1-s1(n))*etam) - (dt/1)*(1 - (1-s1(n))*eta(m))*(1-s1(n))*dsdt(n) + (1-(1-s1(n))*eta(m))*(1-s1(n))^2); 
        Tri_diagL(m,m+1) = (be*Zgap/Pes)*(1 - (1-s1(n))*etap) - (cfl/2)*eta(m)*dsdt(n)*(1 - (1-s1(n))*eta(m))*(1-s1(n)); 
        %Tri_diagR(m,m) = - (dt1/2)*(1 - (1-s1(n))*xi(m))*(1-s1(n))*dsdt1(n) - (1-(1-s1(n))*xi(m))*(1-s1(n))^2;
        Tri_diagR(m,m) = - (1-(1-s1(n))*eta(m))*(1-s1(n))^2;
    end
    
    % discretise bc at xi = 1
    Tri_diagL(I+1,I+1) = 1;  Tri_diagL(I+1,I) = 0;
    Tri_diagR(I+1,I+1) = 0;  Tri_diagR(I+1,I) = 0;
    
    X = zeros(I+1,1);
    X(1) = -2*kappas*deta/(1+beta+chi*rhow1(n)); %BC
    
    % solve for new temperatures
    F1new = Tri_diagL\(Tri_diagR*F1old + X);
     
    dsdt(n+1) = (St*Zgap/Pes)*(3*F1new(I+1) - 4*F1new(I) + F1new(I-1))/2/deta;
    
    F1old = F1new;
    F1(:,n+1) = F1new;
        
    s1(n+1) = s1(n) + dt*dsdt(n+1);
        
end


for n = 1:N
zgapint(n) = -4*nu*s1(n)^2*(-(kappas*((1-s1(n))*F1(1,n)-1))/(1+beta+chi*rhow1(n))/Pes+s1(n)*dsdt(n)/Zgap/St)...
    /((1+nu)*s1(n)^2+1-nu)/(1-s1(n)^2)/(1-nu) - (2*nu*s1(n)^2*P1+(s1(n)^2+1)*qwdot0)/...
    (((1+nu)*s1(n)^2+1-nu))+ (1-s1(n))*(F1(1,n)-F1(1,n+1))/dt/Zgap-dsdt(n)*F1(1,n)/Zgap;
v1(n) = -(kappas*((1-s1(n))*F1(1,n)-1))/(1+beta+chi*rhow1(n));
SigmaRint(n) = ((-0.2e1 * St * (P1 - qwdot0 / 0.2e1) * s1(n) ^ 2 + 0.2e1 * s1(n) * dsdt(n)/Zgap - St * qwdot0) * Pes + 0.2e1 * v1(n) * St) / St / ((1 + nu) * s1(n) ^ 2 + 0.1e1 - nu) / Pes;
end

G = P0 -  trapz(t(1:N),Zgap*zgapint(1:N))

H = trapz(t(1:N),SigmaRint(1:N))
SigmaR = -P0+H

Zgap;

figure(1)
plot(Zgap*t*L*1000,Tmelt-Tscale*(1-s1).*F1(1,:),'-','Linewidth',1.2)
hold on
xlabel('$z$ [mm]')
ylabel('$T$ [K]')

figure(2)
plot(Zgap*t*L*1000,(1-s1)*W,'k-','Linewidth',1.2)
hold on
xlabel('$z$ [mm]')
ylabel('$r$ [mm]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zgap < Z < 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zmin = Zgap;


dZ = (1-Zmin)/NN;
Z = Zmin:dZ:1;

cfl = dZ/deta;
be = dt/deta^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialize matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fold = zeros(I + 1,1);
Fnew = zeros(I + 1,1);
F = zeros(I+1,NN+1);

%thIC = 1 - phi*exp(-tmin/Zc);

s(1) = s1(end);
dsdZ(1) = dsdt(end)/Zgap;

for i = 1:I+1
    F(i,1) = F1(i,end);
    Fold(i) = F(i,1);
end

for n = 1:NN+1
    %theta0(n) = 1 - phi*exp(-t(n)/Zc);
    rhow(n) = M*Z(n);
end
Thd = -(F(2,1) - F(1,1))/deta;
%Thd = -(3*F(3,1) - 4*F(2,1) + F(1,1))/2/deta;
Thdt = Thd;
Thd = -dsdZ(1)*F(1,1) + (1-s(1))*(3*F(1,3)-4*F(1,2)+F(1,1))/2/dZ;
if modelII == 1
    Ud(1) = 2 * (s(1) ^ 2 * P1 - (-(3*F(3,1)-4*F(2,1)+F(1,1))/2/deta) / Pes - s(1) * dsdZ(1) / St) / (-s(1) ^ 2 + 1);
else
    Ud(1) = 4*nu*(s(1)^2)*(-(3*F(3,1)-4*F(2,1)+F(1,1))/2/deta + Pes*s(1)*dsdZ(1)/St)/Pes/(1-nu)/(1-s(1)^4) ...
    + 2*nu*(s(1)^2)*P1/(1+s(1)^2) - (((1+nu)*s(1)^2 + 1 - nu)/(1-nu)/(1+s(1)^2))*Thd;
end

U(1) = -rhow(1)+1e-9 ;

%    Ud(1) = -10;
%    U(1) = -0.00004;

rhoa(1) = -(1)*U(1);
a(1) = rhoa(1) - rhow(1);
%SigmaRdot(1) = ((((-Thd + (-2 * nu + 2) * P1) * s(1) ^ 2 + Thd) * St - 2 * s(1) * dsdZ(1)) * Pes - 2 * -(3*F(3,1)-4*F(2,1)+F(1,1))/2/deta * St) / (-1 + nu) / (s(1) ^ 2 + 1) / Pes / St;

for n = 1:NN
    
    Tri_diagL = sparse(I + 1, I + 1);     
    Tri_diagR = sparse(I + 1, I + 1); 
      
    % discretise bc at eta = 0
    Tri_diagL(1,1) = gamma*(1 +beta +chi*rhow(n)) - (1)*U(n) + rhow(1) - 1*rhow(n);  
    Tri_diagL(1,2) = -4*(gamma*(1 +beta +chi*rhow(n)) - (1)*U(n) +rhow(1) - 1*rhow(n))...
        - 2*gamma*kappas*deta*(1-s(n)); 
    Tri_diagL(1,3) = 3*(gamma*(1 +beta +chi*rhow(n)) - (1)*U(n) +rhow(1) - 1*rhow(n));
    
    
    for m = 2:I      
        etam = (1/2)*(eta(m) + eta(m-1));
        etap = (1/2)*(eta(m) + eta(m+1));
        Tri_diagL(m,m-1) = (be/Pes)*(1 - (1-s(n))*etam) + (cfl/2)*eta(m)*dsdZ(n)*(1 - (1-s(n))*eta(m))*(1-s(n)); 
        Tri_diagL(m,m) = - ( (be/Pes)*(1 - (1-s(n))*etap) + (be/Pes)*(1 - (1-s(n))*etam) - (dZ/1)*(1 - (1-s(n))*eta(m))*(1-s(n))*dsdZ(n) + (1-(1-s(n))*eta(m))*(1-s(n))^2); 
        Tri_diagL(m,m+1) = (be/Pes)*(1 - (1-s(n))*etap) - (cfl/2)*eta(m)*dsdZ(n)*(1 - (1-s(n))*eta(m))*(1-s(n));
        Tri_diagR(m,m) = - (1-(1-s(n))*eta(m))*(1-s(n))^2;
    end
    
    
    Tri_diagL(I+1,I+1) = 1;
    
    X = zeros(I+1,1);
    X(1) = -2*gamma*kappas*deta;
    
    Fnew = Tri_diagL\(Tri_diagR*Fold + X);
     
    dsdZ(n+1) = (St/Pes)*(1*Fnew(I+1) - 1*Fnew(I) )/1/deta;
    
    Fold = Fnew;
    F(:,n+1) = Fnew;
    s(n+1) = s(n) + dZ*dsdZ(n+1);
    
    
    th1 = -(3*F(3,n+1) - 4*F(2,n+1) + F(1,n+1))/2/deta;
     
    Thd = -dsdZ(n+1)*F(1,n+1) + (1-s(n+1))*(F(1,n+1)-F(1,n))/dZ;
    if modelII ==1
        Ud(n+1) = 2 * (s(n+1) ^ 2 * P1 - th1 / Pes - s(n+1) * dsdZ(n+1) / St) / (-s(n+1) ^ 2 + 1);
    else
        Ud(n+1) = 4*nu*(s(n+1)^2)*(th1 + Pes*s(n+1)*dsdZ(n+1)/St)/Pes/(1-nu)/(1-s(n+1)^4) + 2*nu*(s(n+1)^2)*P1/(1+s(n+1)^2) -...
            (((1+nu)*s(n+1)^2 + 1 - nu)/(1-nu)/(1+s(n+1)^2))*Thd;
    end
    U(n+1) = trapz(Z(1:n+1),Ud(1:n+1)) - 1*rhow(1);
    
    %SigmaRdot(n+1) = ((((-Thd + (-2 * nu + 2) * P1) * s(n) ^ 2 + Thd) * St - 2 * s(n) * dsdZ(n)) * Pes - 2 * th1 * St) / (-1 + nu) / (s(n) ^ 2 + 1) / Pes / St;
    SigmaRdot(n+1) = ((((-Thd + (-2 * nu + 2) * P1) * s(n) ^ 2 + Thd) * St - 2 * dsdZ(n) * s(n)) * Pes - 2 * th1 * St) / (-1 + nu) / Pes / St / (s(n) ^ 2 + 1);                   
%       Ud(n+1) = -10;
%       U(n+1) = -00.4;
    % U(2) = -0.1
    rhoa(n+1) = -(1/1)*U(n+1);
    %rhow(n+1) = rhoa(n+1);  %% include for ideal taper, comment out for linear taper
    a(n+1) = rhoa(n+1) - rhow(n+1);
    
    
end

SigmaRBig = trapz(Z(1:NN+1),SigmaRdot(1:NN+1));

rm = s;
ra = 1 - delta*rhoa;
rw = 1 - delta*rhow;
rw1 = 1-delta*rhow1;
r = 1- (1-s(end))*eta;
theta = (1-s(end))*F(:,end);

rwD1 = W*rw1;
rmD = W*rm;
raD = W*ra;
rwD = W*rw;
aD = W*delta*a;
rD = W*r;
tD = L*t;
TT = Tmelt-Tscale*(1-s).*F(1,:);
figure(3)
plot(Z*L*1000,Tmelt-Tscale*(1-s).*F(1,:),'k-','Linewidth',1.2)
hold on
%plot(Z*L*1000,Tmelt-Tscale*(1-s).*F(1,:),'r*')
%hold on
xlabel('$z$ [mm]')
ylabel('$T_{r = r_a}$ [K]')


figure(4)
plot(Z*L*1000,(1-s)*W*1000,'k-','Linewidth',1.2)
ylim([0 W*1000/4])
hold on
xlabel('$z$ [mm]')
ylabel('$W-r_m$ [mm]')
ylim([0 30])

figure(5)
plot(Z*L*1000,(1-ra)*W*1000,'k-','Linewidth',1.2)
ylim([0 3])
hold on
xlabel('$z$ [mm]')
ylabel('$W-r_a$ [mm]')

Zt = [t*Zgap Z(2:end)];
Ft = [F1(1,:) F(1,2:end)];
st = [s1 s(2:end)];
rwt = [rwD1 rwD(2:end)];

figure(6)
plot(Zt*L*1000,Tmelt-Tscale*(1-st).*Ft-273,'k-','Linewidth',1.2)
hold on
xlabel('$z$ [mm]')
ylabel('$T$ [K]')


figure(7)
hold on
plot(Z*L*1000,(W-rwD)*1000,'k--','Linewidth',1.2)
plot(Z*L*1000,(W-raD)*1000,'k-','Linewidth',1.2)
%ylim([0 W])
hold on
xlabel('$z$ [mm]')
ylabel('$r$ [mm]')


figure(8)
plot(Z*L*1000,aD*1000,'k--','Linewidth',1.2)
hold on
xlabel('$z$ [mm]')
ylabel('$a_g$ [mm]')
TTs_ra = (Tmelt-Tscale*(1-s).*F(1,:));

for n = 1:N+1
   hh1(n) = -kM/(rwD1(n)*(log(rwD1(n)/WM)-kM/(hC*WM)));
%    Trw(n) = Zgap*t(n)*L*1000,Tmelt-Tscale*(1-s1(n)).*F1(1,n);
%    AM1 (n) = (Trw(n) - Tw)/(log(rwD1(n)/WM)-kM/(hC*WM));
end

for n = 1:NN+1
%     AAM(n) = ka*(TTs_ra(n) - TT0(n))/(ka*log(W/WM) - kM*log(W/rraD(n)));
%     QQ(n) = kM*AAM(n)/W;
    
    hh(n) = ka*kM/(raD(n)*(kM*log(rwD(n)/raD(n))-ka*(log(rwD(n)/WM)-kM/(WM*hC))));
    AM(n) = ka*(TTs_ra(n)-Tw)/(ka*(log(rwD(n)/WM)-kM/WM/hC)-kM*log(rwD(n)/raD(n)));
    Aa(n) = AM(n)*kM/ka;
    TTm_WM(n) = AM(n)*(log(((WM)/WM))-kM/WM/hC)+Tw;
    QQ(n) = kM*ka*(TTs_ra(n)-Tw)/(raD(n)*(ka*(log(rwD(n)/WM)-kM/WM/hC)-kM*log(rwD(n)/raD(n))));
    QWM (n) = kM*AM(n)/WM; 
    TTm1(n) = AM(n)*(log(((rwD(n)+0.00275)/WM))-kM/WM/hC)+Tw;
    TTm2(n) = AM(n)*(log(((rwD(n)+0.00875)/WM))-kM/WM/hC)+Tw;
    TTm3(n) = AM(n)*(log(((WM)/WM))-kM/WM/hC)+Tw;
end

hht = [hh1 hh(2:end)];

figure(9)
plot(Z*L,hh,'k-','Linewidth',1.2)
hold on
xlabel('$z$ [mm]')
ylabel('$h$ [Wm$^{-2}$K$^{-1}$]')
xlim([-L/10 L])

III = trapz(Z*L,TTm_WM-Tw);
DDT = 2*WM*hC*III/cpw/rhoww/Vw/DE/(DE+2*WM)
figure(10)
hold on
plot(Z*L,TTm_WM,'k-','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('$T_{r = W_{M}}$ [K]')
xlim([-L/10 L])

DEt = linspace(1/100,1,100);
Vwt = muw*1e4/rhoww./DEt;
% figure(11)
% plot(DEt,Vwt,'k-')

figure(12)
hold on
plot((W-raD)*1000,Z*L*1000,'k--','Linewidth',1.2)
plot((W-rmD)*1000,Z*L*1000,'k:','Linewidth',1.2)
plot(zeros(1,length(Z)),Z*L*1000,'k-','Linewidth',1.2)
plot((W-rwD)*1000,Z*L*1000,'r','Linewidth',1.2)
plot(W*1000*ones(1,length(Z)),Z*L*1000,'k-','Linewidth',1.2)
% plot(-1000*(WM-W)*ones(1,length(Z)),Z*L*1000,'k-')
% plot(-1000*((WM+DE)-W)*ones(1,length(Z)),Z*L*1000,'k-')
%plot(-(WM+DE)*ones(1,length(Z)),Z*L*1000,'k-')
set(gca,'Ydir','reverse')
xlim([-((W/10)*1000)-10 W*1000])
xlabel('$r$ [mm]')
ylabel('$z$ [mm]')

figure(13)
semilogy(Z*L*1000,hh,'k-','Linewidth',1.2)
hold off
hold on
xlabel('$z$ [mm]')
ylabel('$h$ [Wm$^{-2}$K$^{-1}$]')
xlim([-L*1000/10 L*1000])

figure(14)
semilogy(Z*L*1000,QQ,'k-','Linewidth',1.2)
hold off
hold on
xlabel('$z$ [mm]')
ylabel('$Q$ [Wm$^{-2}$]')
xlim([-L/10*1000 L*1000])

figure(15)
semilogy(Zt*L*1000,hht,'k-','Linewidth',1.2)
hold off
hold on
xlabel('$z$ [mm]')
ylabel('$h$ [Wm$^{-2}$K$^{-1}$]')
xlim([-L*1000/10 L*1000])

figure(16)
hold on
plot(Z*L,TTm1-273,'k-','Linewidth',1.2)
plot(Z*L,TTm2-273,'k--','Linewidth',1.2)
plot(Z*L,TTm3-273,'k:','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('$T_{r = r_m(z)+\epsilon}$ [K]')
xlim([0.01 L])
%ylim([273 573])

figure(25)
hold on
plot(Z*L*1000,TTm1,'k-','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('Temperature [K]')
xlim([0 L*1000])
ylim([0+273 500+273])

figure(26)
hold on
plot(Z*L*1000,TTm2,'k-','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('Temperature [K]')
xlim([0 L*1000])
ylim([0+273 500+273])

figure(27)
semilogy(Z*L*1000,QWM,'k-','Linewidth',1.2)
hold off
hold on
xlabel('$z$ [mm]')
ylabel('$Q$ [Wm$^{-2}$]')
xlim([-L/10*1000 L*1000])

figure(28)
hold on
plot(L*Z*1000,1000*(W-W*(1-delta*rhow)),'k','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('$W-r_w$ [mm]')

figure(29)
hold on
plot(Z,rhow,'k','Linewidth',1.2)
% figure(28)
% plot(L*Z,W*(1-delta*rhow(Z)))
xlabel('$Z$ ')
ylabel('$\varrho_w$')
figure(30)
hold on
plot((W-rwD)*1000,Z*L*1000,'k','Linewidth',1.2)
set(gca,'Ydir','reverse')
xlabel('$W-r_w$ [mm]')
ylabel('$z$ [mm]')

figure(31)
hold on
plot(L*Z*1000,1000*(W-W*(1-delta*rhow)),'k','Linewidth',1.2)
plot(L*Z*1000,1000*(W-W*(1-delta*rhoa)),'--r','Linewidth',1.2)
xlabel('$z$ [mm]')
ylabel('$W-r_w$ [mm]')

end
figure(3)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,1200,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

figure(4)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,15,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

figure(5)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,1,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

figure(6)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,1200,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

figure(8)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,0.8,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

figure(10)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(0.5,1200,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])


figure(13)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,1200,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])
xlim([0 100])

figure(14)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,-1e5,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])
ylim([-5e6 -0.5e4])

figure(15)
str1 = '$ V_w = 1, 2, 4, 12$ ms$^{-1}$';
text(500,100,str1,'Interpreter','latex','Fontsize',10)
annotation('arrow','X',[0.51,0.35],'Y',[0.51,0.35])

rhowater = 997;
Vwt = linspace(0.8,20,100);
Ret = rhowater.*Vwt.*DE./muw;
Prt = cpw.*muw./kw;
Nut1 = 0.023.*Ret.^(0.8).*Prt.^(0.4);
c1t = 0.88-0.24./(4+Prt);
c2t = 1/3+0.5.*exp(-0.6.*Prt);
Nut2 = 5+0.015.*Ret.^c1t.*Prt.^c2t;
hCt1 = kw./DE.*Nut1;
hCt2 = kw./DE.*Nut2;

rhowater = 997;
DED = linspace(0.01,0.05,100);
RetD = rhowater.*Vw.*DED./muw;
PrtD = cpw.*muw./kw;
Nut3 = 0.023.*RetD.^(0.8).*PrtD.^(0.4);
c1t = 0.88-0.24./(4+PrtD);
c2t = 1/3+0.5.*exp(-0.6.*PrtD);
Nut4 = 5+0.015.*RetD.^c1t.*PrtD.^c2t;
hCt3 = kw./DE.*Nut3;
hCt4 = kw./DE.*Nut4;

Ret = linspace(3e4,5e5,1000);
Nu1 = zeros(1,length(Ret));
Nu2 = zeros(1,length(Ret));
for i = 1:length(Ret)
    Nu1(i) = 0.023*Ret(i)^(0.8)*Pr^(0.4);
     c1 = 0.88-0.24/(4+Pr);
     c2 = 1/3+0.5*exp(-0.6*Pr);
    Nu2(i) = 5+0.015*Ret(i)^c1*Pr^c2;
end
Dvec = linspace(0.001,0.1,100);
Rev = zeros(1,length(Dvec));
for i = 1:length(Dvec)
    Rev(i) = 997*Vw*Dvec(i)/muw;
end
figure(21)
hold on
plot(Ret,Nu1,'-','Linewidth',1.3)
plot(Ret,Nu2,'--','Linewidth',1.3)
xlabel('$ Re $')
ylabel(' $ Nu$')
bleg = legend('Dittus-Boelter','Sleicher \& Rouse'); 
set(bleg,'Interpreter','latex','Fontsize',11,'Location','northwest');


MM = csvread('KellyData3.csv');

M8 =  csvread('KellyTemp8_75.csv');
M2 =  csvread('KellyTemp2_75.csv');
M8e = csvread('KellyTempExp8_75.csv');
M2e = csvread('KellyTempExp2_75.csv');

figure(8)
hold on
plot(MM(:,1)-51,MM(:,2),'-','color',[0.047 0.341 0.847])
M8 =  csvread('KellyTemp8_75.csv');
M2 =  csvread('KellyTemp2_75.csv');
xlim([-10 790])

figure(16)
hold on
plot(M8(:,1)/1000-0.0556,M8(:,2),'b','Linewidth',1.2)
plot(M2(:,1)/1000-0.0556,M2(:,2),'r','Linewidth',1.2)

figure(25)
hold on
plot((M2(:,1)-51)/1,M2(:,2)+273,'-','color',[0.980 0.282 0.235],'Linewidth',1.2)
plot((M2e(:,1)-51)/1,M2e(:,2)+273,'x','color',[0.980 0.282 0.235],'Linewidth',1.2)
ylim([0+273 400+273])
xlim([-10 790])
figure(26)
hold on
plot((M8(:,1)-51)/1,M8(:,2)+273,'-','color',[0.047 0.341 0.847],'Linewidth',1.2)
plot((M8e(:,1)-51)/1,M8e(:,2)+273,'+','color',[0.047 0.341 0.847],'Linewidth',1.2)
ylim([0+273 400+273])
xlim([-10 790])
T0min = Tw;
T00 = 781.5;
dT00 = -3e5;

zc = -(T00-T0min)/dT00;
phi = (T00-T0min)/(Tmelt-T0min);
Ttest = T0min + (T00-T0min).*exp(-Z.*L./zc);

figure(10)
hold on
plot(Z*L,Ttest,'--r')

toc