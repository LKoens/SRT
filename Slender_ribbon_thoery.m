function [F,T,a] = Slender_ribbon_thoery(r1,r2,r3,t1,t2,t3,Th1,Th2,Th3, rho, bl, Ud1, Ud2,UR,OR)
%SLENDER_RIBBON_THOERY Single program to do all the slender ribbon thoery
%compuations at a single time step.
%IMPUTS
% ri  - the ith component of the ribbon vector as function of s
% ti  - the ith component of tangent vector to the centreline of the ribbon as function of s
% Thi  - the ith component of direction of ribbon plane vector as function of s
% rho - cross-section of ribbon shape as function of s
% bl - ribbon width/ ribbon length
% Um - centreline velocity as function of s
% Ur - rotation velocity as function of s
% OUTPUTS
% F - total force on ribbon
% T -total torque on ribbon
% a - distribution of the force along centreline (legendre polynomials)

N = 20; %number of legendre polynomials to use
Mat =zeros(3*N); %Matrix containing all the terms
v= zeros(3*N,1); %vector for all the velocity terms

Um=@(x) UR + Ud1(x)+ cross([r1(x),r2(x),r3(x)],OR); %velcity not depdnent on s2
Ur=@(x) Ud2(x) +bl.*cross(rho(x).*[Th1(x),Th2(x),Th3(x)],OR); %veclotiy dependent on s2

% U integral
for m=0:(N-1)
    v((1+3*m):(3+3*m))=integral( @(x) 8*LegendreP(m,x).*Um(x),-1,1,'ArrayValued',true); %complete LHS
    
    if m==0
        Lm=0;
    else
        Lm=sum( 1./(1:m));
    end
    
    for n=0:m
        %Identity integral
        vImn= integral(@(x) log((16*(1-x.^2))./((bl^2)*(rho(x).^2))).*LegendreP(n,x).*LegendreP(m,x), -1,1);
        
        Imn = vImn*eye(3);
        
        % tt integrals
        vtmn11=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t1(x),-1,1);
        vtmn12=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t2(x),-1,1);
        vtmn13=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t1(x).*t3(x),-1,1);
        vtmn22=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t2(x).*t2(x),-1,1);
        vtmn23=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t2(x).*t3(x),-1,1);
        vtmn33=integral(@(x) (log((16*(1-x.^2))./((bl^2)*(rho(x).^2)))-2).*LegendreP(n,x).*LegendreP(m,x).*t3(x).*t3(x),-1,1);
        
        tmn=zeros(3);
        tmn(1,1)=vtmn11;
        tmn(2,2)=vtmn22;
        tmn(3,3)=vtmn33;
        tmn(1,2)=vtmn12;
        tmn(2,1)=vtmn12;
        tmn(1,3)=vtmn13;
        tmn(3,1)=vtmn13;
        tmn(2,3)=vtmn23;
        tmn(3,2)=vtmn23;
        
        % TT integrals
        vTmn11=integral(@(x) 2*Th1(x).*Th1(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        vTmn12=integral(@(x) 2*Th1(x).*Th2(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        vTmn13=integral(@(x) 2*Th1(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        vTmn22=integral(@(x) 2*Th2(x).*Th2(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        vTmn23=integral(@(x) 2*Th2(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        vTmn33=integral(@(x) 2*Th3(x).*Th3(x).*LegendreP(n,x).*LegendreP(m,x),-1,1);
        
        Tmn=zeros(3);
        Tmn(1,1)=vTmn11;
        Tmn(2,2)=vTmn22;
        Tmn(3,3)=vTmn33;
        Tmn(1,2)=vTmn12;
        Tmn(2,1)=vTmn12;
        Tmn(1,3)=vTmn13;
        Tmn(3,1)=vTmn13;
        Tmn(2,3)=vTmn23;
        Tmn(3,2)=vTmn23;
        
        %Put matrix components together of parts that are symetric to the
        %interchange of m and n
        if m==n
            Mat((1+3*m):(3+3*m),(1+3*n):(3+3*n))= Imn +  tmn +Tmn - 2*Lm*eye(3)/(2*m+1);
        else
            Mat((1+3*m):(3+3*m),(1+3*n):(3+3*n))= Imn +  tmn +Tmn;
            Mat((1+3*n):(3+3*n),(1+3*m):(3+3*m))= Imn +  tmn +Tmn;
        end
        
    end
end

for m=0:(N-1)
    for n=0:(N-1)
        % Non-local integrals
        vkmnia=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*(1./sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)-1./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmnib=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*(1./sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)-1./sqrt((x-y).^2)),-1,1,@(x) x,1);
       
        vkmn11a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r1(x)-r1(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t1(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn11b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r1(x)-r1(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t1(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        vkmn12a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t2(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn12b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t2(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        vkmn13a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn13b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r1(x)-r1(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t1(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        vkmn22a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t2(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn22b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r2(x)-r2(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t2(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        vkmn23a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn23b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r2(x)-r2(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t2(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        vkmn33a=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r3(x)-r3(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t3(x).*t3(x)./sqrt((x-y).^2)),-1,1,-1,@(x) x);
        vkmn33b=quad2d(@(y,x) LegendreP(n,y).*LegendreP(m,x).*((r3(x)-r3(y)).*(r3(x)-r3(y))./((sqrt((r1(x)-r1(y)).^2+(r2(x)-r2(y)).^2 +(r3(x)-r3(y)).^2)).^3)-t3(x).*t3(x)./sqrt((x-y).^2)),-1,1,@(x) x,1);
        
        Kmn=zeros(3);
        Kmn(1,1)=vkmnia+vkmnib + vkmn11a+ vkmn11b;
        Kmn(2,2)=vkmnia+vkmnib + vkmn22a+ vkmn22b;
        Kmn(3,3)=vkmnia+vkmnib + vkmn33a+ vkmn33b;
        Kmn(1,2)=vkmn12a+ vkmn12b;
        Kmn(2,1)=vkmn12a+ vkmn12b;
        Kmn(1,3)=vkmn13a+ vkmn13b;
        Kmn(3,1)=vkmn13a+ vkmn13b;
        Kmn(2,3)=vkmn23a+ vkmn23b;
        Kmn(3,2)=vkmn23a+ vkmn23b;
        
        %Matrix compoent together (add last bit)
        Mat((1+3*m):(3+3*m),(1+3*n):(3+3*n))=Mat((1+3*m):(3+3*m),(1+3*n):(3+3*n))+Kmn;
        
    end
end


a=Mat\v;%Solve for the coefficents

F=2*pi*a(1:3); %total force



Tint= pi.*bl.*integral(@(x)cross(rho(x).*[Th1(x),Th2(x),Th3(x)],(4*(eye(3) - kron([t1(x), t2(x), t3(x)],[t1(x),t2(x),t3(x)]')/2)*Ur(x)')'),-1,1,'ArrayValued',true)./2;

Tc=0;
for m=0:(N-1)
    Tc=Tc+pi*integral(@(x) cross([r1(x),r2(x),r3(x)],a((1+3*m):(3+3*m)).*LegendreP(m,x)),-1,1,'ArrayValued',true);
end

T=Tc+Tint;  %total torque
end

