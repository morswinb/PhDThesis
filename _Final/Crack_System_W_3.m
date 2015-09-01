classdef Crack_System_W_3<handle
    
    properties,
        function_calls,q_l,xi,N,k,M,L_handle,left_handle,right_handle,Initial_Condition,diffs,Leak_handle,L0;
        open_array,oN,iterations_exeded,leak_type,CL,da;
    end
    
    methods
        % constructor, przypisanie zmienych
        % use_preconditoning
%         use_preopened
%         leak_type 0 no leak 1 carter 2 p carter
        function this=Crack_System_W_3(Initial_Condition,diffs,...
                use_preconditoning,use_preopened,leak_type,L_factor,L_extra)
            
            this.leak_type=leak_type;
            this.N=Initial_Condition.grid.N;
            this.Initial_Condition=Initial_Condition;
            this.Initial_Condition.alpha=1/3;
            this.Initial_Condition.alpha=1/2;  
            
            if(this.leak_type==0);
                this.Initial_Condition.alpha=4/3;  
            elseif(this.leak_type==1);
                this.Initial_Condition.alpha=1/2;  
            elseif(this.leak_type==2);
                this.Initial_Condition.alpha=5/6;  
            end
           

            this.L_handle=dL_2_handle(Initial_Condition);
            this.left_handle=left_handle_1(Initial_Condition);
            this.right_handle=right_handle_1(Initial_Condition);
            this.diffs=diffs;
            this.function_calls=0;
            this.xi=Initial_Condition.grid.xi;
            this.k=Initial_Condition.k;
            this.M=Initial_Condition.M;    
            this.oN=this.N;
            this.open_array=ones(1,this.N);
            
            this.CL=1;
            
            
            this.da=AsymFD(Initial_Condition,Initial_Condition.alpha,Initial_Condition.beta);
            
%             this.Leak_handle=Leak_handle_Carter_fast4a(Initial_Condition);
%             this.Leak_handle=Leak_handle_Carter_fast_BACK3(Initial_Condition);
%             this.Leak_handle=Leak_handle_Carter_fast_BACK4(Initial_Condition,use_preconditoning,use_preopened);

            this.Leak_handle=LeakOff.LeakLoader.loadBackLeak();
%             this.Leak_handle.setUp(this.xi,Initial_Condition.a,use_preconditoning,use_preopened,L_factor,L_extra);
%             Initial_Condition.a
%             pause
            
        end
        
		
        function sol=solve(this,t_spam)
    
            w=this.Initial_Condition.w;
            L_t=this.Initial_Condition.L_0();
            this.L0=L_t;
			
            %szukanie V_0
            t=0;
            dwdx=this.diffs.calcFirstDer(w);
            d2wdx2=this.diffs.calcSecondDer(w);                 
            [dwdx(1),d2wdx2(1)]=this.left_handle.get_left(t,L_t,w,dwdx,d2wdx2);
            [w_0,~,~]=this.right_handle.get_right(t,L_t,w,dwdx,d2wdx2);
            V_0=this.k/3/this.M/L_t*w_0^3;
            
            this.Leak_handle.add_event_L(V_0,L_t,0);

            n=this.N;
            A=diag(ones(1,n+1),0)+diag(ones(1,n),1)+diag(ones(1,n),-1);
            A(:,end)=1;
            A(:,end-1)=1;
            A(:,end-2)=1;
            A(1,3)=1;
            
            A(end+1,end+1)=1;
            A(end,end-1)=1;
            A(end-1,end)=1;
            
            options = odeset('RelTol',1e-10,'AbsTol',1e-10,...
                'Stats','off','Jconstant','off','BDF','on',...
                'Jpattern',A,'Events',@(t,y)this.event1(t,y)); %%%,'NonNegative',ones(n+1,1));    
            
%             figure
%             this.Leak_handle.printLeak(0,this.xi);
%             hold on;
            
            try

                i_sol=ode15s(@(t,w)ODE(this,t,w),t_spam,...
                    [this.Initial_Condition.w;this.L_handle.get_L_0();0],options);

                sol=i_sol;
                sol.y=sol.y(:,:);
%                 sol.y=sol.y(n+1:n+2,:);
                while(sol.x(end)~=t_spam(end))
                    
                    tt=[sol.x(end) t_spam(end)];
                                        
                    w=this.remesh(sol.x(end),i_sol.y(1:n,end));

                    
                    L_t=this.xi(end)*sol.y(n+1,end).^.5;
                    tttt=sol.x(end);
                    dwdx=this.diffs.calcFirstDer(w);
                    d2wdx2=this.diffs.calcSecondDer(w);                 
                    [dwdx(1),d2wdx2(1)]=this.left_handle.get_left(tttt,L_t,w,dwdx,d2wdx2);
                    [w_0,~,~]=this.right_handle.get_right(tttt,L_t,w,dwdx,d2wdx2);
                    V_0=this.k/3/this.M/L_t*w_0^3;
                    
                    this.Leak_handle.add_event_L(V_0,L_t,sol.x(end));
%                     this.Leak_handle.printLeak(sol.x(end),this.xi);
                    
%                     this.Leak_handle.add_Initial_L(V_0,L_t,0)
                    
                    i_sol=ode15s(@(t,w)ODE(this,t,w),tt,[w;L_t*L_t;i_sol.y(n+2,end)],options);
                               
                    if(this.iterations_exeded)
                        break; 
                    end
                    if(abs(i_sol.x(end)-sol.x(end))/i_sol.x(end)<10^-8)
                        break;
                    end

                    
                    sol.x=[sol.x i_sol.x];
                    sol.y=[sol.y i_sol.y];

%                     sol.y=[sol.y i_sol.y(n+1:n+2,:)];

                end
                
            catch e
                e.message
                if(isa(e, 'matlab.exception.JavaException'))
                    ex = e.ExceptionObject;
                    assert(isJava(ex));
                    ex.printStackTrace;
                end
            end
            
        end

        function dw=ODE(this,t,w)
            
%            t
           
           
           this.function_calls=this.function_calls+1;
           if(this.function_calls>100000)
               this.iterations_exeded=1;
               throw(MException('Id:id','iterations exeeded'));
           end
           
           % przepisanie wartosci do lokalnych zmienych
           n=this.N;
           ww=w(1:n+1);
           w=0;
           w=ww;
           L_t=this.L_handle.get_L(w);
           x=this.xi;

           
%            fprintf('L0=%f L=%f t %f \n',this.L0,L_t,t); 
           
           % przygotowanie miejsca na odpowiedz
           dw=zeros(n+1,1);        
           
           %policznie pochodnych po x na przedziale
%            dwdx=this.diffs.calcFirstDer(w);
%            d2wdx2=this.diffs.calcSecondDer(w);
            this.da.preset(w(1:n),t,L_t);
            dwdx=this.da.calcFirstDer();
            d2wdx2=this.da.calcSecondDer();  

           %policzenie pochodnych na granicach przedzia?u, i w_0
           [dwdx(1),d2wdx2(1)]=this.left_handle.get_left(t,L_t,w,dwdx,d2wdx2);
           [w_0,dwdx(n),d2wdx2(n)]=this.right_handle.get_right(t,L_t,w,dwdx,d2wdx2);   
                   
           V_0=this.k/3/this.M/L_t*w_0^3;
           this.Leak_handle.add_L(V_0,L_t,t);
           
           %polidzenie dL/dt
           dw(n+1)=this.L_handle.get_dL(w_0,w);
           
           %zastosowanie rówania na dw/dt (w,x) na ca?ym przedziale
           dw(1:n)=this.k/this.M/L_t^2*(1/3*w_0^3*x(1:n).*dwdx(1:n)...
               +3*w(1:n).^2.*(dwdx(1:n).^2)+w(1:n).^3.*d2wdx2(1:n));
            
           
           if(this.leak_type==0)
                leak=zeros(n,1);
           elseif(this.leak_type==1)
%                leak=this.Leak_handle.get_Leak(L_t,t,x);
                leak=this.Leak_handle.get_Leak(L_t,t)*this.CL;
           elseif(this.leak_type==2)
                leak=this.Leak_handle.get_Leak(L_t,t);
%                 leak=this.Leak_handle.get_Leak(L_t,t,x);
                leak=leak.*this.k.*w(1:n)*this.CL;
           elseif(this.leak_type==3)
               leak=this.Initial_Condition.q_l(t,x(1:n));
           end
%            elseif(this.leak_type==4)
%                 leak=this.Leak_handle.get_Leak(L_t,t);
% %                leak=this.Leak_handle.get_Leak(L_t,t)+this.Initial_Condition.bench.ql_star(t,x);
%            end
           
%            leak
%            pause
                     
%             leak=this.Leak_handle.get_Leak(t);

%             if(t>10^-8)
%                 figure
%                 plot(x,leak);
%                 
%                 pause
%             end


           dw(1:n)=dw(1:n)-leak;
           dw(n+2)=trapz(x,leak)*L_t;
                      
%            dw(1:n)=dw(1:n).*this.open_array';
%             fprintf('t %f\n',t);

        end
        
        function w_new=remesh(this,t,w)
            
            n=this.N;
            L_t=this.L_handle.get_L(w);
            dwdx=this.diffs.calcFirstDer(w);
            d2wdx2=this.diffs.calcSecondDer(w);                 
            [dwdx(1),d2wdx2(1)]=this.left_handle.get_left(t,L_t,w,dwdx,d2wdx2);
            [w_0,dwdx(n),d2wdx2(n)]=this.right_handle.get_right(t,L_t,w,dwdx,d2wdx2);           
                      
            a=this.Initial_Condition.alpha;
            b=this.Initial_Condition.beta;
            [A B A1 A2 B1 B2]=this.right_handle.get_AB(w);
            V_0=this.k/3/this.M/L_t*w_0^3;
            
            f=@(x) A*(1-x).^(a)+B*(1-x).^(b);
            df=@(x) -a*A*(1-x).^(a-1)-b*B*(1-x).^(b-1);
            dl=dwdx(1)-df(0);
            dr=dwdx(n)-df(this.xi(n));
            wr=w(1:this.N)-f(this.xi);
            Sp=spline(this.xi,[dl;wr;dr]);
            old_w=@(x) f(x)+fnval(Sp,x);
            w(1:n)=old_w(this.xi(1:n)*this.xi(end));            
                     

            dwdx=this.diffs.calcFirstDer(w);
            d2wdx2=this.diffs.calcSecondDer(w);                           
            [dwdx(1),d2wdx2(1)]=this.left_handle.get_left(t,L_t*this.xi(end),w,dwdx,d2wdx2);
            [w_0,dwdx(n),d2wdx2(n)]=this.right_handle.get_right(t,L_t*this.xi(end),w,dwdx,d2wdx2);   
            [k kk A1 A2 B1 B2]=this.right_handle.get_AB(w);
            V_0=this.k/3/this.M/L_t*this.xi(end)*w_0^3;  

            w_new=real(w);
        end
        
        function [value,isterminal,direction] = event1(this,t,y) %->
            value(1)=y(this.N);
            isterminal=[1];
            direction=[-1];
        end
        
    end
    
end

