function H = cylinder_model_zgap_fn(x)

global I N St Pes kappas beta Zc M deta eta nu P1 gamma P0 t1min s1 dsdt F1 qwdot0 chi

Zgap = x;
tmin = 1e-13;

dt = (1-tmin)/N;
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
    F1(i,1) = kappas/(1+beta+chi*rhow1(1))*(1-eta(i)); %initial condition
    F1old(i) = F1(i,1);
end



for n = 1:N
    
    Tri_diagL = sparse(I + 1, I + 1);     
    Tri_diagR = sparse(I + 1, I + 1); 
      
    % discretise bc at xi = 0
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

H = trapz(t(1:N),SigmaRint(1:N))
G = P0 -  trapz(t(1:N),Zgap*zgapint(1:N));
