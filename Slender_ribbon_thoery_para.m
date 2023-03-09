function [R,F,T,a,Mat,Imnf,ILmnf,tmnf,tLmnf,Tmnf,NLa,NLb] = Slender_ribbon_thoery_para(r1,r2,r3,t1,t2,t3,Th1,Th2,Th3, rho,bl, Ud1, Ud2)
%SLENDER_RIBBON_THOERY Single program to do all the slender ribbon thoery
%compuations at a single time step. parallel process formualtion.
%IMPUTS
% ri  - the ith component of the ribbon vector as function of s
% ti  - the ith component of tangent vector to the centreline of the ribbon as function of s
% Thi  - the ith component of direction of ribbon plane vector as function of s
% bl - ribbon width/ ribbon length
%rho - cross section of ribbon
% Ud1 - defromation velocity as function of s
% Ud2 - deformation velocity as function of s multiplied by s2
% OUTPUTS
% F - total force on ribbon
% T -total torque on ribbon
% a - distribution of the force along centreline (legendre polynomials)

N = 30; %number of legendre polynomials to use
Imnf=zeros(3*N); %Matrix containing Identy terms
ILmnf=zeros(3*N); %Matrix containing Identy times log terms
tmnf=zeros(3*N); %Matrix containing tt terms
tLmnf=zeros(3*N);%Matrix containing tt time log terms
Tmnf=zeros(3*N); %Matrix containing TT terms
NLa =zeros(3*N); %Matrix containing non-local eigenfunction terms
NLb =zeros(3*N); %Matrix containing non-local intergal terms
R=zeros(6); % resistance matrix


% Integreals to evaluate with m and n
% U integral
vILmn = zeros(N*(N+1)/2);
vtmn11 = zeros(N*(N+1)/2);
vtmn12 = zeros(N*(N+1)/2);
vtmn13 = zeros(N*(N+1)/2);
vtmn22 = zeros(N*(N+1)/2);
vtmn23 = zeros(N*(N+1)/2);
vtmn33 = zeros(N*(N+1)/2);
vtLmn11 = zeros(N*(N+1)/2);
vtLmn12 = zeros(N*(N+1)/2);
vtLmn13 = zeros(N*(N+1)/2);
vtLmn22 = zeros(N*(N+1)/2);
vtLmn23 = zeros(N*(N+1)/2);
vtLmn33 = zeros(N*(N+1)/2);
vTmn11 = zeros(N*(N+1)/2);
vTmn12 = zeros(N*(N+1)/2);
vTmn13 = zeros(N*(N+1)/2);
vTmn22 = zeros(N*(N+1)/2);
vTmn23 = zeros(N*(N+1)/2);
vTmn33 = zeros(N*(N+1)/2);

%local components

parfor q=1:(N*(N+1)/2) %evaluate local integrals
    m = floor((sqrt(8*(q-1)+1)-1)/2);
    n=(q-1)-m*(m+1)/2;
    
    %Identity integral
    vILmn(q)= integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x), -1,1);
    
    % tt integrals
    vtmn11(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t1(x).*t1(x),-1,1);
    vtmn12(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t1(x).*t2(x),-1,1);
    vtmn13(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t1(x).*t3(x),-1,1);
    vtmn22(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t2(x).*t2(x),-1,1);
    vtmn23(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t2(x).*t3(x),-1,1);
    vtmn33(q)=integral(@(x) LegendreP(n,x).*LegendreP(m,x).*t3(x).*t3(x),-1,1);
    
    vtLmn11(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t1(x),-1,1);
    vtLmn12(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t2(x),-1,1);
    vtLmn13(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t3(x),-1,1);
    vtLmn22(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t2(x).*t2(x),-1,1);
    vtLmn23(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t2(x).*t3(x),-1,1);
    vtLmn33(q)=integral(@(x) log( (16*(1-x.^2))./((bl^2)*(rho(x).^2)) ).*LegendreP(n,x).*LegendreP(m,x).*t3(x).*t3(x),-1,1);
    
    % TT integrals
    vTmn11(q)=integral(@(x) Th1(x).*Th1(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    vTmn12(q)=integral(@(x) Th1(x).*Th2(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    vTmn13(q)=integral(@(x) Th1(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    vTmn22(q)=integral(@(x) Th2(x).*Th2(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    vTmn23(q)=integral(@(x) Th2(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    vTmn33(q)=integral(@(x) Th3(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
    
end



for q=1:(N*(N+1)/2) %assign values into matrices
    m = floor((sqrt(8*(q-1)+1)-1)/2);
    n=(q-1)-m*(m+1)/2;
    
    
    
    
    ILmn = vILmn(q)*eye(3);
    
    
    tmn=zeros(3);
    tmn(1,1)=vtmn11(q);
    tmn(2,2)=vtmn22(q);
    tmn(3,3)=vtmn33(q);
    tmn(1,2)=vtmn12(q);
    tmn(2,1)=vtmn12(q);
    tmn(1,3)=vtmn13(q);
    tmn(3,1)=vtmn13(q);
    tmn(2,3)=vtmn23(q);
    tmn(3,2)=vtmn23(q);
    
    tLmn=zeros(3);
    tLmn(1,1)=vtLmn11(q);
    tLmn(2,2)=vtLmn22(q);
    tLmn(3,3)=vtLmn33(q);
    tLmn(1,2)=vtLmn12(q);
    tLmn(2,1)=vtLmn12(q);
    tLmn(1,3)=vtLmn13(q);
    tLmn(3,1)=vtLmn13(q);
    tLmn(2,3)=vtLmn23(q);
    tLmn(3,2)=vtLmn23(q);
    
    Tmn=zeros(3);
    Tmn(1,1)=vTmn11(q);
    Tmn(2,2)=vTmn22(q);
    Tmn(3,3)=vTmn33(q);
    Tmn(1,2)=vTmn12(q);
    Tmn(2,1)=vTmn12(q);
    Tmn(1,3)=vTmn13(q);
    Tmn(3,1)=vTmn13(q);
    Tmn(2,3)=vTmn23(q);
    Tmn(3,2)=vTmn23(q);
    
    if m==n
        Imnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= 2/(2*m+1)*eye(3);
        
        ILmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= ILmn;
        
        tmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= tmn;
        tLmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= tLmn;
        Tmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= Tmn;
        
    else
        ILmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= ILmn;
        ILmnf((1+3*n):(3+3*n),(1+3*m):(3+3*m))= ILmn;
        
        tmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= tmn;
        tmnf((1+3*n):(3+3*n),(1+3*m):(3+3*m))= tmn;
        
        tLmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= tLmn;
        tLmnf((1+3*n):(3+3*n),(1+3*m):(3+3*m))= tLmn;
        
        Tmnf((1+3*m):(3+3*m),(1+3*n):(3+3*n))= Tmn;
        Tmnf((1+3*n):(3+3*n),(1+3*m):(3+3*m))= Tmn;
    end
end

%Non-local terms

for m=0:(N-1) %determin non-local eigenfunction matrix
    
    if m==0
        Lm=0;
    else
        Lm=2*sum( 1./(1:m));
    end
    
    NLa(:,(1+3*m):(3+3*m))=-Lm*(Imnf(:,(1+3*m):(3+3*m))+tmnf(:,(1+3*m):(3+3*m)));

end


vkmni =zeros(N*N);
vkmn11 =zeros(N*N);
vkmn12 =zeros(N*N);
vkmn13 =zeros(N*N);
vkmn22 =zeros(N*N);
vkmn23 =zeros(N*N);
vkmn33 =zeros(N*N);

parfor q=1:(N*N) %determine non-local integral terms
    n=floor((q-1)/N)
    m=rem(q-1,N)
    
    vkmnia=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*(1./sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)-1./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmnib=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*(1./sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)-1./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmni(q)=vkmnia+vkmnib;
    
    vkmn11a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r1(x)-r1(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t1(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn11b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r1(x)-r1(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t1(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn11(q)=vkmn11a+vkmn11b;
    
    vkmn12a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t2(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn12b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t2(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn12(q)=vkmn12a+vkmn12b;
    
    vkmn13a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn13b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn13(q)=vkmn13a+vkmn13b;
    
    vkmn22a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t2(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn22b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t2(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn22(q)=vkmn22a+vkmn22b;
    
    vkmn23a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn23b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn23(q)=vkmn23a+vkmn23b;
    
    vkmn33a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r3(x)-r3(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t3(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
    vkmn33b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r3(x)-r3(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t3(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
    vkmn33(q)=vkmn33a+vkmn33b;
end

for q=1:N*N
    n=floor((q-1)/N);
    m=rem(q-1,N);
    
    Kmn=zeros(3);
    Kmn(1,1)=vkmni(q) + vkmn11(q);
    Kmn(2,2)=vkmni(q) + vkmn22(q);
    Kmn(3,3)=vkmni(q) + vkmn33(q);
    Kmn(1,2)=vkmn12(q);
    Kmn(2,1)=vkmn12(q);
    Kmn(1,3)=vkmn13(q);
    Kmn(3,1)=vkmn13(q);
    Kmn(2,3)=vkmn23(q);
    Kmn(3,2)=vkmn23(q);
    
    NLb((1+3*m):(3+3*m),(1+3*n):(3+3*n))=Kmn;
end

Mat=NLa+NLb+ILmnf+tLmnf-2*tmnf+2*Tmnf;

%the velocity evalatuions
v= zeros(3*N,1); %vector deformation terms
v1= zeros(3*N,1); %vector for rigid body motion in 1
v2= zeros(3*N,1); %vector for rigid body motion in 2
v3= zeros(3*N,1); %vector for rigid body motion in 3
v4= zeros(3*N,1); %vector for rigid body rotation in 1
v5= zeros(3*N,1); %vector for rigid body rotation in 2
v6= zeros(3*N,1); %vector for rigid body rotation in 1

Um1=@(x) [1,0,0] + cross([0,0,0],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2
Um2=@(x) [0,1,0] + cross([0,0,0],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2
Um3=@(x) [0,0,1] + cross([0,0,0],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2
Um4=@(x) [0,0,0] + cross([1,0,0],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2
Um5=@(x) [0,0,0] + cross([0,1,0],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2
Um6=@(x) [0,0,0] + cross([0,0,1],[r1(x),r2(x),r3(x)]); %velcity not depdnent on s2

Ur1=@(x) bl.*cross([0,0,0],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2
Ur2=@(x) bl.*cross([0,0,0],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2
Ur3=@(x) bl.*cross([0,0,0],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2
Ur4=@(x) bl.*cross([1,0,0],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2
Ur5=@(x) bl.*cross([0,1,0],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2
Ur6=@(x) bl.*cross([0,0,1],rho(x).*[Th1(x),Th2(x),Th3(x)]); %veclotiy dependent on s2

for m=0:(N-1)
     v((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Ud1(x),-1,1,'ArrayValued',true); %complete LHS
     
    v1((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um1(x),-1,1,'ArrayValued',true); %complete LHS
    
    v2((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um2(x),-1,1,'ArrayValued',true); %complete LHS
    
    v3((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um3(x),-1,1,'ArrayValued',true); %complete LHS
    
    v4((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um4(x),-1,1,'ArrayValued',true); %complete LHS
    
    v5((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um5(x),-1,1,'ArrayValued',true); %complete LHS
    
    v6((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um6(x),-1,1,'ArrayValued',true); %complete LHS
end


%Solve for the coefficents of the force distubtion
a=Mat\v;

%total force
F=2*pi*a(1:3); 

%total torque calculations
Tint= pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ud2(x)')'),-1,1,'ArrayValued',true)./2;

Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end

T=Tc+Tint;  

%resistance matrix calcuations

a1=Mat\v1;
R(1:3,1)=2*pi*a1(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur1(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a1((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,1)=Tc+Tint;

a2=Mat\v2;
R(1:3,2)=2*pi*a2(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur2(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a2((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,2)=Tc+Tint;

a3=Mat\v3;
R(1:3,3)=2*pi*a3(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur3(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a3((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,3)=Tc+Tint;


a4=Mat\v4;
R(1:3,4)=2*pi*a4(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur4(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a4((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,4)=Tc+Tint;

a5=Mat\v5;
R(1:3,5)=2*pi*a5(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur5(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a5((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,5)=Tc+Tint;

a6=Mat\v6;
R(1:3,6)=2*pi*a6(1:3); 
Tint =pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur6(x)')'),-1,1,'ArrayValued',true)./2;
Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a6((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end
R(4:6,6)=Tc+Tint;

end

