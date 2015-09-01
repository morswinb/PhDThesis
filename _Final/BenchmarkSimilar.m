classdef BenchmarkSimilar
   properties
       alpha,beta,w,V,L,dL,q_l,int_q_l,q_0,dwdt,Omega,k,M,a,w_0,...
           A,C,D,powers,sym_powers,L_inv,q_l_star,q_l_norm,dwdx,...
           Yb,YY,ql_star;
   end
    
   methods
       
       function obj = BenchmarkSimilar(a)
           
           obj.a=a;
				 		   
		   kl=4;
		   ke=1;
		   q0=1;
		   qn=1;
           tn=1;
           
		   xn=((kl*ke)/4)^(1/5)*qn^(3/5)*tn^(4/5);
		   wn=qn*tn/xn;
                      
		   b0=1;
		   b1=-1/16;
		   b2=-15/224*b1;
		   b3=-3/80*b2;
		   b4=-11/5824*b3;
		   		   
		   xi_s=1.3208446*(q0/qn)^(0.6);

           obj.Yb=@(x)0.6*xi_s^2*((1-x).^(1)...
               +b1/2*(1-x).^(2)+b2/3*(1-x).^(3)...
               +b3/4*(1-x).^(4)+b4/5*(1-x).^(5));
           

 		   Y=@(x)0.6*xi_s^2*(b0/1*(1-x).^(1)+b1/2*(1-x).^(2)+...
 			   +b2/3*(1-x).^(3)+b3/4*(1-x).^(4)+b4/5*(1-x).^(5));                        
           obj.YY=Y;

           
           dxYb=@(x)0.6*xi_s^2*(...
               -1 ...
               +2*1/16/2*(1-x).^(1)+...
			   +3*15/224*-1/16/3*(1-x).^(2)...
               +4*3/80*-15/224*-1/16/4*(1-x).^(3)...
               +5*11/5824*-3/80*-15/224*-1/16/5*(1-x).^(4));
                                
		   obj.k=4;	   
		   obj.M=1;
		   obj.w=@(t,x)wn*(Y(x)).^(1/3).*(a+t).^(1/5);
           obj.dwdt=@(t,x)1/5*wn*(Y(x)).^(1/3).*(a+t).^(-4/5);
           obj.dwdx=@(t,x)...
               wn*1/3*(Y(x)).^(-2/3).*dxYb(x).*(a+t).^(1/5);
		   obj.q_0=@(t)q0;	
           obj.ql_star=@(t,x) t*0+x*0;
 		   obj.L=@(t)xi_s*xn*(a+t).^(4.0/5.0);
           obj.V=@(t,x)-obj.k/obj.M./obj.L(t).*...
               obj.w(t,x).^2.*obj.dwdx(t,x);
		   obj.L_inv=@(x)(x./xi_s./xn).^(5.0/4.0)-a;  
		   obj.alpha=1.0/3.0;  
		   obj.beta=4.0/3.0;
           obj.Omega=@(t,x)quad(@(x)obj.w(t,x),1,x,10^-12);  
		   obj.q_l=@(t,x)0+0*t+0*x;                     
           
       end
   end
    
end



