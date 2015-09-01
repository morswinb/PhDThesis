classdef Benchmark_CC_2

   properties
       alpha,beta,w,V,L,dL,q_l,int_q_l,q_0,dwdt,Omega,k,M,a,w_0,...
           A,C,D,powers,sym_powers,L_inv,q_l_star,q_l_norm,dVdt,dwdx,dVdx,dwdxt
       ql_star;
       
%               alpha,beta,A,Y,D,C,C_2,dVdx,V1,A_1,U_0,...
%            w,V,L,L_inv,dL,q_l,q_l_star,q_l_norm,q_l_b,q_l_a,q_l_c,int_q_l,q_0,...
%            dOmegadt,dwdx,d2wdx,dwdt,Omega,k,M,a,...
%            powers,sym_powers;
   end
    
   methods
       
       function obj = Benchmark_CC_2(k,M,a,ver,ratio)
           
           version=ver;
           %1 h=@(x)(1-x).^(1/3)+c1*(1-x).^(1/2)+c2*(1-x).^(4/3);
           %2 h=@(x)(1-x).^(1/3)+c1*(1-x).^(4/3)+c2*(1-x).^(7/3);
           %3 h=@(x)(1-x).^(1/3)+c1*(1-x).^(5/6)+c2*(1-x).^(4/3);
           
           obj.k=k;
           obj.M=M;
           obj.a=a;
           
           y=1/3;
%            u_0=1.2;
           u_0=(3/2*(3*y+1))^(1/3);
           
           if(ratio==1)
               % 1:1 ratio (0.992922)
%                C1=.1;
%                C2=.5;
% 
              C1=0.19;
              C2=0.41;

%               C1=0.74;
%               C2=-0.13;

%               C1=0.02;
%               C2=0.1;

%               C1=0.03;
%               C2=0.08;

%               C1=0.15;
%               C2=-0.02;


           elseif(ratio==2)
               C1=.021;
               C2=.1;
           elseif(ratio==3)
               C1=1;
               C2=1;
           end
           
%             C1=.0125;
%             C2=.9;
           

%            C1=1/8/exp(1)*(1/3-2*y/(3*y+1));
%            C2=-0.05;
           
           obj.alpha=1/3;
           switch version
               case 1
                   obj.beta=1/2;
               case 3
                   obj.beta=4/3;
               case 2
                   obj.beta=5/6;
           end           

           syms('x');
           w_0=(1-x).^(1/3)+C1*(1-x).^(1/2)+C2*(1-x).^(4/3);
           dw_0=diff(w_0,'x');
           ddw_0=diff(diff(w_0,'x'),'x');
           syms('t');
           syms('x');
           sq_l=u_0*(a+t).^(y-1).*(-y*w_0+3/2*(3*y+1)*...
               (1/3*x.*dw_0+3*w_0.^2.*dw_0.^2+w_0.^3.*ddw_0));
           sq_l=int(sq_l,'x');
           sq_l=matlabFunction(sq_l,'vars',{'t','x'});
           obj.int_q_l=@(t,x)sq_l(t,x)-sq_l(t,0);
           
           ww_0=matlabFunction(w_0);
           dw_0=matlabFunction(diff(w_0,'x'));
           obj.q_0=@(t)-(3*k*u_0^5*(3*y+1)/2/M)^.5*...
               (a+t).^((5*y-1)/2)*ww_0(0)^3*dw_0(0);
                    
           if(version==1)
              
               sym_w_0=sym(w_0);
               w_0=matlabFunction(sym_w_0);
               dw_0=matlabFunction(diff(sym_w_0,'x'));
               ddw_0=matlabFunction(diff(diff(sym_w_0,'x'),'x'));

               obj.w=@(t,x)u_0*(a+t).^y.*w_0(x);
               obj.dwdt=@(t,x)u_0*y*(a+t).^(y-1).*w_0(x);
               obj.dwdx=@(t,x)u_0*(a+t).^y.*dw_0(x);
               obj.dwdxt=@(t,x)u_0*y*(a+t).^(y-1).*dw_0(x);
%                obj.d2wdx=@(t,x)u_0*y*(a+t).^y.*ddw_0(x);

               obj.L=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y+1)/2);
               obj.dL=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y-1)/2)*((3*y+1)/2);

               obj.L_inv=@(x)(x.*sqrt((3*M*(3*y+1))/(2*k*u_0^3))).^(2/(3*y+1))-a;

               obj.D=@(t)sqrt((3*y+1)/(2*(a+t)));

               obj.q_0=@(t)-(3*k*u_0^5*(3*y+1)/2/M)^.5*...
                   (a+t).^((5*y-1)/2)*...
                   w_0(0)^3*dw_0(0);

               obj.q_l=@(t,x)u_0*(a+t).^(y-1).*(...
                       -y*w_0(x)+3/2*(3*y+1)*...
                       (1/3*x.*dw_0(x)+3*w_0(x).^2.*dw_0(x).^2+w_0(x).^3.*ddw_0(x)));
                   
               obj.ql_star=@(t,x)obj.q_l(t,x)-1./sqrt(t-obj.L_inv(x*obj.L(t)));
                   
                   
                obj.V=@(t,x)-k/M./obj.L(t).*obj.w(t,x).^2.*obj.dwdx(t,x);

%                obj.q_l_norm=@(t,x)obj.C./sqrt(t-obj.L_inv(obj.L(t).*x));
%                obj.q_l_star=@(t,x)obj.q_l(t,x)-obj.q_l_norm(t,x);
% 
%                obj.V=@(t,x)-k/M*u_0^3./obj.L(t).*(a+t).^(3*y).*w_0(x).^2.*dw_0(x);
%                
%  			   obj.dVdt=@(t,x)-k/M*u_0^3.*w_0(x).^2.*dw_0(x).*...
% 				   (3*y*(a+t).^(3*y-1)./obj.L(t)-2*(a+t).^(3*y)/obj.L(t).^2*obj.dL(t));
			   
               obj.dVdx=@(t,x)-k/M*u_0^3./obj.L(t).*(a+t).^(3*y).*...
                   (2*w_0(x).*dw_0(x).^2+w_0(x).^2.*ddw_0(x));
%                obj.V1=@(t)k/3/M./obj.L(t).*u_0^3.*(a+t).^(3*y);

               w_0_int=int(sym_w_0,'x');
               fff=matlabFunction(w_0_int);
               obj.Omega=@(t,x)u_0*(a+t).^y.*fff(x);
%                obj.dOmegadt=@(t,x)u_0*y*(a+t).^(y-1).*fff(x);

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
                             
           else
               q_0s=obj.q_0(0);
               Q_ls=obj.int_q_l(0,1);
                              
               syms x c1 c2;
               if(version==3)
                   sym_w_0=(1-x).^(1/3)+c1*(1-x).^(4/3)+c2*(1-x).^(7/3);
               else
                   sym_w_0=(1-x).^(1/3)+c1*(1-x).^(5/6)+c2*(1-x).^(4/3);
               end              
               
               w_0=matlabFunction(sym_w_0,'vars',{'x','c1','c2'});
               dw_0=matlabFunction(diff(sym_w_0,'x'),'vars',{'x','c1','c2'});
%                ddw_0=matlabFunction(diff(diff(sym_w_0,'x'),'x'),'vars',{'x','c1','c2'});
               
               q_0_cc=@(t,c1,c2)-(3*k*u_0^5*(3*y+1)/2/M)^.5*...
                   (a+t).^((5*y-1)/2)*w_0(0,c1,c2)^3*dw_0(0,c1,c2);               
               syms t;
               w_0=sym_w_0;
               dw_0=diff(sym_w_0,'x');
               ddw_0=diff(diff(sym_w_0,'x'),'x');
               sq_l=u_0*(a+t).^(y-1).*(-y*w_0+3/2*(3*y+1)*...
                   (1/3*x.*dw_0+3*w_0.^2.*dw_0.^2+w_0.^3.*ddw_0));
               sq_l=int(sq_l,'x');
               sq_l=matlabFunction(sq_l,'vars',{'t','x','c1','c2'});
               int_q_l_cc=@(t,x,c1,c2)sq_l(t,x,c1,c2)-sq_l(t,0,c1,c2);
               
%                disp('starting least square')
               [c,resnorm] = ...
                   lsqnonlin(@(cc)Benchmark_CC_2.cc_min_fun(q_0_cc,int_q_l_cc,q_0s,Q_ls,cc),[1 1]);
%                fprintf('found C1=%f C2=%f with squared 2-norm %f',c(1),c(2),resnorm);

               w_0=matlabFunction(sym_w_0,'vars',{'x','c1','c2'});
               dw_0=matlabFunction(diff(sym_w_0,'x'),'vars',{'x','c1','c2'});
               ddw_0=matlabFunction(diff(diff(sym_w_0,'x'),'x'),'vars',{'x','c1','c2'});
               w_0=@(x)w_0(x,c(1),c(2));
               dw_0=@(x)dw_0(x,c(1),c(2));
               ddw_0=@(x)ddw_0(x,c(1),c(2));
               
               obj.w=@(t,x)u_0*(a+t).^y.*w_0(x);
               obj.dwdt=@(t,x)u_0*y*(a+t).^(y-1).*w_0(x);
%                obj.dwdx=@(t,x)u_0*y*(a+t).^y.*dw_0(x);
%                obj.d2wdx=@(t,x)u_0*y*(a+t).^y.*ddw_0(x);
               obj.L=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y+1)/2);
               obj.dL=@(t)(2*k*u_0^3/3/M/(3*y+1))^.5*(a+t).^((3*y-1)/2)*((3*y+1)/2);
               obj.L_inv=@(x)(x.*sqrt((3*M*(3*y+1))/(2*k*u_0^3))).^(2/(3*y+1))-a;
               obj.D=@(t)sqrt((3*y+1)/(2*(a+t)));
               
               w_0_int=int(sym_w_0,'x');
               fff=matlabFunction(w_0_int,'vars',{'x','c1','c2'});
               fff=@(x)fff(x,c(1),c(2));
               obj.Omega=@(t,x)u_0*(a+t).^y.*fff(x);

%                obj.dOmegadt=@(t,x)u_0*y*(a+t).^(y-1).*fff(x);
               
               obj.q_0=@(t)q_0_cc(t,c(1),c(2));
               obj.q_l=@(t,x)u_0*(a+t).^(y-1).*(-y*w_0(x)+3/2*(3*y+1)*...
                   (1/3*x.*dw_0(x)+3*w_0(x).^2.*dw_0(x).^2+w_0(x).^3.*ddw_0(x)));
               
               obj.V=@(t,x)-k/M*u_0^3./obj.L(t).*(a+t).^(3*y).*w_0(x).^2.*dw_0(x);


               w_0=sym_w_0;
               dw_0=diff(sym_w_0,'x');
               ddw_0=diff(diff(sym_w_0,'x'),'x');
               syms('t');
               syms('x');
               sq_l=u_0*(a+t).^(y-1).*(-y*w_0+3/2*(3*y+1)*...
                   (1/3*x.*dw_0+3*w_0.^2.*dw_0.^2+w_0.^3.*ddw_0));
               sq_l=int(sq_l,'x');
               sq_l=matlabFunction(sq_l,'vars',{'t','x','c1','c2'});
               obj.int_q_l=@(t,x)sq_l(t,1,c(1),c(2))-sq_l(t,x,c(1),c(2));
               
%                obj.q_l_norm=@(t,x)obj.C./sqrt(t-obj.L_inv(obj.L(t).*x));
%                obj.q_l_star=@(t,x)obj.q_l(t,x)-obj.q_l_norm(t,x);
               
%                q_0s
%                Q_ls
%                q_0_cc
%                int_q_l_cc
           end        
           
           obj.w_0=@(t)u_0*(a+t).^y;
       end
   end
   
   methods(Static)
       
       function f=cc_min_fun(q_0_cc,int_q_l_cc,q_0s,Q_ls,cc)
           f(1)=q_0_cc(0,cc(1),cc(2))-q_0s;
           f(2)=int_q_l_cc(0,1,cc(1),cc(2))-Q_ls;
       end
       
   end
    
end



