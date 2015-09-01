classdef Nowa_Siatka_2
    
    properties
        N,xi,dxi,epsilon,xi_n_plus_one,dxi_n;
    end
    
    methods
        function this=Nowa_Siatka_2(N,epsilon,alpha,beta)

            this.N=N;
            this.epsilon=epsilon;

            delta=epsilon;
            
% %             alpha=3;
%             beta=1.5;
            eps=delta;
            M=10*N+1;
            x=0:1/(M-1):1;
            xx=0:1/(N-1):1;
            xx=(1-eps^(1/alpha))*xx;

            mm=floor(M/2)+1;

            AAA=(1-(x(mm))^alpha)/(x(mm))^beta;

            for ii=1:M
                if ii<M/2
                    w(ii)=AAA*x(ii)^beta;
                else
                    w(ii)=1-(1-x(ii))^alpha;
                end
            end

            for ii=1:M
                if ii<M/2
                    w1(ii)=beta*AAA*x(ii)^(beta-1);
                else
                    w1(ii)=alpha*(1-x(ii))^(alpha-1);
                end
            end

            r1=round(M/6);
            r2=round(M/5);

            r3=round(2*M/3);
            r4=round(7*M/8);
            in=0;
            for ii=1:M
                if x(ii)<x(r2)
                    in=in+1;
                    xe(in)=x(ii);
                    we(in)=AAA*x(ii)^beta;
                elseif x(ii)>x(r3)
                    in=in+1;
                    xe(in)=x(ii);
                    we(in)=1-(1-x(ii))^alpha;

                end
            end

            y=spaps(xe,we,0,3); 

            for ii=1:length(xx)
                if xx(ii)==0
                    wg(ii)=0;
                else
                    wg(ii)=fnval(xx(ii),y);
                end
            end

            ww=wg;
            
            this.xi=ww';
            this.dxi=diff(this.xi);
%             this.dxi_n=this.xi_n_plus_one-this.xi(end);
        end
    end
end

