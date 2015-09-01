classdef Benchmark_Custom_2
    %Benchmark_Carter zestaw wszystkich funkcj benchmarka i parametórw
    %   Detailed explanation goes here
    
   properties
       alpha,beta,A,Y,D,C,C_2,dVdx,V1,A_1,U_0,...
           w,V,L,L_inv,dL,q_l,q_l_star,q_l_norm,q_l_b,q_l_a,q_l_c,int_q_l,q_0,...
           dOmegadt,dwdx,d2wdx,dwdt,Omega,k,M,a,...
           powers,sym_powers;
   end
    
   methods
       
       function obj = Benchmark_Custom_2(k,M,a)
           
           obj.k=k;
           obj.M=M;
           obj.a=a;
           
           
%            1/3+2/3=1
%            1/3+1/2=5/6
%            1/3+1/6=1/2

           
           y=1/5;
           obj.Y=y;
           obj.U_0=1;
           u_0=obj.U_0;    
           obj.C=1;
           
           A_j=sym([1 1 1]);
           alpha_j=sym([1/3 1/2 4/3]);
%            alpha_j=sym([1/3 5/6 4/3]);
%            alpha_j=sym([1/3 1 5/3]);
           
           obj.powers=double(alpha_j);
           obj.A=double(A_j);
           obj.sym_powers=alpha_j;
           

           sym_w_0=strcat(char(A_j(1)),'*(1-x)^(',char(alpha_j(1)),')');
           for i=2:length(A_j)
              sym_w_0=strcat(sym_w_0,'+',char(A_j(i)),'*(1-x)^(',char(alpha_j(i)),')');
           end
           
           obj.alpha=double(alpha_j(1));
           obj.beta=double(alpha_j(2));
           
           
           sym_w_0=sym(sym_w_0);
           
           w_0=matlabFunction(sym_w_0)
           dw_0=matlabFunction(diff(sym_w_0,'x'));
           ddw_0=matlabFunction(diff(diff(sym_w_0,'x'),'x'));

%            w_0=@(x)(1-x).^(1/3)+a_1*(1-x).^(1/2);
%            dw_0=@(x)-1/3*(1-x).^(-2/3)-1/2*a_1*(1-x).^(-1/2);
%            ddw_0=@(x)-2/9*(1-x).^(-5/3)-1/4*a_1*(1-x).^(-3/2);
           

           obj.w=@(t,x)u_0*(a+t).^y.*w_0(x);
           obj.dwdt=@(t,x)u_0*y*(a+t).^(y-1).*w_0(x);
           obj.dwdx=@(t,x)u_0*y*(a+t).^y.*dw_0(x);
           obj.d2wdx=@(t,x)u_0*y*(a+t).^y.*ddw_0(x);
           
           obj.L=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y+1)/2);
           obj.dL=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y-1)/2)*((3*y+1)/2);
           
           obj.L_inv=@(x)(x.*sqrt((3*M*(3*y+1))/(2*k*u_0^3))).^(2/(3*y+1))-a;

           
%            obj.C=@(t)(a+t)^(1/2);
           %obj.C=@(t)A_j(1)*(3*y+1)^(1/2)*7/8*sqrt(2)*(a+t)^(y-1/2);
           obj.D=@(t)sqrt((3*y+1)/(2*(a+t)));
           
%            obj.q_0=@(t)-u_0^(5/2)*(a+t).^((5*y-1)/2)/(2*M/3/k/(3*y+1))^.5*...
%                w_0(0)^3*dw_0(0);

           obj.q_0=@(t)-(3*k*u_0^5*(3*y+1)/2/M)^.5*...
               (a+t).^((5*y-1)/2)*...
               w_0(0)^3*dw_0(0);
           
           obj.q_l=@(t,x)u_0*(a+t).^(y-1).*(...
                   -y*w_0(x)+3/2*(3*y+1)*...
                   (1/3*x.*dw_0(x)+3*w_0(x).^2.*dw_0(x).^2+w_0(x).^3.*ddw_0(x)));
               
           obj.q_l_norm=@(t,x)obj.C./sqrt(t-obj.L_inv(obj.L(t).*x));
           obj.q_l_star=@(t,x)obj.q_l(t,x)-obj.q_l_norm(t,x);

               
%            obj.q_l=@(t,x)u_0*(a+t).^(y-1).*((3*y+1)*...
%                (7/8*a_1*(1-x).^(-1/2)+5/2*a_1^2*(1-x).^(-1/3)...
%                +55/24*a_1^3*(1-x).^(-1/6)+3/4*a_1^4*(1-x).^(0))+...
%                1/6*(1-3*y)*(1-x).^(1/3)+1/4*(1-y)*a_1*(1-x).^(1/2));

%            obj.q_l=@(t,x)u_0*(a+t).^(y-1).*(-y*w_0(x)+...
%                3/2*(3*y+1)*(1/3*x.*dw_0(x)+3*w_0(x).^2.*dw_0(x).^2+...
%                w_0(x).^3.*ddw_0(x)));
           
%            obj.q_l_c=@(t,x)u_0*(a+t).^(y-1).*(...
%                -y*(1-x).^(1/3)-y*a_1*(1-x).^(1/2)+3/2*(3*y+1)*...
%                (1/9*(1-x).^(1/3)+1/6*a_1*(1-x).^(1/2)+...
%                7/12*a_1*(1-x).^(-1/2)+5/3*a_1^2*(1-x).^(-1/3)+...
%                55/36*a_1^3*(1-x).^(-1/6)+1/2*a_1^4*(1-x).^(0));
           
           
           obj.V=@(t,x)-k/M*u_0^3./obj.L(t).*(a+t).^(3*y).*w_0(x).^2.*dw_0(x);
           obj.dVdx=@(t,x)-k/M*u_0^3./obj.L(t).*(a+t).^(3*y).*...
               (2*w_0(x).*dw_0(x).^2+w_0(x).^2.*ddw_0(x));
           obj.V1=@(t)k/3/M./obj.L(t).*u_0^3.*(a+t).^(3*y);
           
           w_0_int=int(sym_w_0,'x');
           fff=matlabFunction(w_0_int);
           obj.Omega=@(t,x)u_0*(a+t).^y.*fff(x);
           obj.dOmegadt=@(t,x)u_0*y*(a+t).^(y-1).*fff(x);
           
           w_0=sym_w_0;
           dw_0=diff(sym_w_0,'x');
           ddw_0=diff(diff(sym_w_0,'x'),'x');
                                   
           syms('t');
           syms('x');
           sq_l=u_0*(a+t).^(y-1).*(-y*w_0+3/2*(3*y+1)*...
               (1/3*x.*dw_0+3*w_0.^2.*dw_0.^2+w_0.^3.*ddw_0));
           sq_l=int(sq_l,'x');
           sq_l=matlabFunction(sq_l,'vars',{'t','x'});
           obj.int_q_l=@(t,x)sq_l(t,1)-sq_l(t,x);
           
% % %            %obj.Omega=@(t,x)u_0*(a+t).^y.*(3/4*(1-x).^(4/3)+2/3*a_1*(1-x).^(3/2));
% % %            obj.dOmegadt=@(t,x)u_0*y*(a+t).^(y-1).*(3/4*(1-x).^(4/3)+2/3*a_1*(1-x).^(3/2));
% % %            
% % % %            obj.int_q_l=@(t,x)u_0*(a+t).^(y-1).*((3*y+1)*...
% % % %                (7/4*a_1*(1-x).^(1/2)+15/4*a_1^2*(1-x).^(2/3)+...
% % % %                11/4*a_1^3*(1-x).^(5/6)+3/4*a_1^4*(1-x).^(1))+...
% % % %                +1/8*(1-3*y)*(1-x).^(4/3)+1/6*(1-y)*a_1*(1-x).^(3/2));
% % %            integral=@(t,x)-(u_0*(a + t).^(y - 1).*(42*a_1*(1 - x).^(1/2) + 4*a_1*(1 - x).^(3/2) - 9*y*(1 - x).^(4/3) + 3*(1 - x).^(4/3) - 18*a_1^4*x + 90*a_1^2*(1 - x).^(2/3) + 66*a_1^3*(1 - x).^(5/6) - 54*a_1^4*x*y + 270*a_1^2*y*(1 - x).^(2/3) + 198*a_1^3*y*(1 - x).^(5/6) + 126*a_1*y*(1 - x).^(1/2) - 4*a_1*y*(1 - x).^(3/2)))/24;
% % %            obj.int_q_l=@(t,x)integral(t,1)-integral(t,x);
           
       end
   end
    
end



