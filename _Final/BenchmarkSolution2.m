classdef BenchmarkSolution2

   properties
       alpha,beta,w,V,L,dL,q_l,int_q_l,q_0,dwdt,Omega,k,M,a,w_0,...
           A,C,D,powers,sym_powers,L_inv,q_l_star,q_l_norm;
       
%               alpha,beta,A,Y,D,C,C_2,dVdx,V1,A_1,U_0,...
%            w,V,L,L_inv,dL,q_l,q_l_star,q_l_norm,q_l_b,q_l_a,q_l_c,int_q_l,q_0,...
%            dOmegadt,dwdx,d2wdx,dwdt,Omega,k,M,a,...
%            powers,sym_powers;
   end
    
   methods
       
       function obj = BenchmarkSolution2(k,M,a)
           
           obj.k=k;
           obj.M=M;
           obj.a=a;
           
           y=1/5;
           ww=1.2;
%            ww=1;

           
%            AA=1;
%            BB=.1;
           
           AA=-1/8/exp(1)*(1/3-2*y/(3*y+1));
           BB=0.05;
           
           obj.alpha=1/3;
           obj.beta=4/3;
           
           obj.w_0=a^y.*ww;
           obj.A=[1/3 4/3];
           
           h=@(x)(1-x).^(1/3)+AA*(1-x).^(4/3)+BB*(1-x).^(7/3);
           HH=@(x)-(3*(1 - x).^(4/3).*(14*BB*(x - 1).^2 - 20*AA*(x - 1) + 35))/140;
           dh=@(x)-1/3*(1-x).^(-2/3)-4/3*AA*(1-x).^(4/3)-7/3*BB*(1-x).^(4/3);
           
           
           
           dh3=@(x)3*(AA*(1 - x).^(4/3) + BB*(1 - x).^(7/3)...
               + (1 - x).^(1/3)).^2.*((4*AA*(1 - x).^(1/3))/3 ...
           + (7*BB*(1 - x).^(4/3))/3 + 1./(3*(1 - x).^(2/3))).^2 ...
           + (AA*(1 - x).^(4/3) + BB*(1 - x).^(7/3) ...
           + (1 - x).^(1/3)).^3.*((4*AA)./(9*(1 - x).^(2/3)) ...
           + (28*BB*(1 - x).^(1/3))./9 - 2./(9*(1 - x).^(5/3)));

           %h=@(x)(1-x).^(1/3).*(1+AA*(1-x)+BB*(1-x).^2);
           %dh=@(x)- (4*AA*(1 - x).^(1/3))/3 - (7*BB*(1 - x).^(4/3))/3 - 1/(3*(1 - x).^(2/3));
           %dh3=@(x) h(x).^2.*dh(x).*(1/3*(1-x).^(-2/3)+16/3*AA*(1-x).^(4/3)+49/3*BB*(1-x).^(4/3));
           %HH=@(x)-3/4*(1-x).^(4/3)-3/7*AA*(1-x).^(7/3)-3/10*BB*(1-x).^(10/3);
           
           B=-sqrt((3*k*ww^3*(3*y+1))/(2*M));
           C=(3*y-1)/2;
           D=sqrt((2*k*ww^3)/(3*M*(3*y+1)));
           E=(3*y+1)/2;
           F=ww*(3*y+1)/2;
           G=-sqrt((3*k*ww*(3*y+1))/(2*M))*ww^2;
           
           H=(5*y-1)/2;
           
           %I=h(0)^3*dh(0); 
           I=-(1/3+4/3*AA+7/3*BB)*(1+AA+BB)^3; 
           
           obj.w=@(t,x)((a+t).^y).*(ww*h(x));
           obj.dwdt=@(t,x)y*((a+t).^(y-1)).*(ww*h(x));
           obj.V=@(t,x)B*(a+t).^C.*(h(x).^2).*dh(x);
           obj.L=@(t)D*(a+t).^E;  
           obj.dL=@(t)D*E*(a+t).^(E-1);
           
           obj.q_l=@(t,x)F*(a+t).^(y-1).*(x.*dh(x)+3*dh3(x))-y*ww*h(x).*(a+t).^(y-1);          
           obj.int_q_l=@(t,x)F*(a+t).^(y-1).*(x.*h(x)-HH(x)+3*h(x).^3.*dh(x))...
               -y*ww*HH(x).*(a+t).^(y-1);    
           obj.q_0=@(t)G*((a+t).^H)*I;
                   
           obj.Omega=@(t,x)((a+t).^y).*ww.*(3/4*(1-x).^(4/3)+...
               AA*3/7*(1-x).^(7/3)+BB*3/10*(1-x).^(10/3));
           
           
       end
   end
    
end



