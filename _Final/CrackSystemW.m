classdef CrackSystemW<handle
    
    properties,
        function_calls,xi,N,k,M,...
            left_BC,right_BC,ic,diffs;
    end
    
    methods
        %constructor
        %ic set of initial conditions
        %left_BC used for BC at x=0
        %right_BC used for BC at x=1
        %diffs used for dw/dx and d2w/dx2 approximation
        function this=CrackSystemW(ic,left_BC,right_BC,diffs)
            this.N=ic.grid.N;
            this.ic=ic;
            this.left_BC=left_BC;
            this.right_BC=right_BC;
            this.diffs=diffs;      
            this.xi=ic.grid.xi;
            this.k=ic.k;
            this.M=ic.M;
        end
        
         
        function [sol]=solve(this,t_spam)
            
            %reset function calls counter
            this.function_calls=0;
            
            n=this.N;
            A=diag(ones(1,n+1),0)+... %tridiagonal part
                diag(ones(1,n),1)+...
                diag(ones(1,n),-1);
            A(:,end-2:end)=1;   %three extra columns
            A(1,3)=1;  %attributed to left BC
            
            options = odeset('RelTol',1e-8,'AbsTol',1e-8,...
                'Jpattern',A); 

            sol=ode15s(@(t,w)ODE(this,t,w),t_spam,...
                [this.ic.w;this.ic.L_0^2],options); 	
			
            
        end
        
        function dw=ODE(this,t,w)
                      
            this.function_calls=this.function_calls+1;
            %allows to escape form endless computation
            if(this.function_calls>100000)
                throw(MException('Id:id','iterations exeeded'));
            end
            
            % assign local variables
            n=this.N;
            L_t=w(end)^.5;
            x=this.xi;
            
            %allocated output of operators A and B
            %{A_1,A_2,...,A_N,B}
            dw=zeros(n+1,1);        
            
            %computing derevative approximations
            this.diffs.preset(w);
            dwdx=this.diffs.calcFirstDer();
            d2wdx2=this.diffs.calcSecondDer();  
                      
            %getting boundary derevative values, and asymptotics term w_0
            [dwdx(1),d2wdx2(1)]=...
                this.left_BC.get_left(t,L_t,w,dwdx,d2wdx2);
            [w_0,dwdx(n),d2wdx2(n)]=...
                this.right_BC.get_right(t,L_t,w,dwdx,d2wdx2);   
            
            %computes operator A
            dw(1:n)=this.k/this.M/L_t^2.*dwdx(1:n).*w(1:n).^3.*...
                (1/3*w_0^3*x(1:n)./w(1:n).^3+3*dwdx(1:n)./w(1:n)...
                +d2wdx2(1:n)./dwdx(1:n))-this.ic.q_l(t,x(1:n));
            
            %computes operator B
            dw(end)=2*this.k/this.M/3*w_0^3;
      
        end  
    end
    
end

