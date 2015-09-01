classdef MixedGrid

    properties
        N,xi,dxi,epsilon,xi_n_plus_one,dxi_n;       
    end
       
    methods    
        function this=MixedGrid(N,epsilon)
            M=N/5;
            
            N=N-1;
            this.N=N;
            this.epsilon=epsilon;
            
            MN=M-N;
            b1=1.0e-3;
            b2=0.999999999999;

            while abs(b2-b1)>eps
                f1=exp(1/MN*log((1-b1)/epsilon))-1+b1/M/(1-b1);
                %f2=exp(1/MN*log((1-b2)/epsilon))-1+b2/M/(1-b2);
                b0=(b1+b2)/2;
                f0=exp(1/MN*log((1-b0)/epsilon))-1+b0/M/(1-b0);
                if  f1*f0>0
                    b1=b0;
                else
                    b2=b0;
                end
            end
            
            b=(b1+b2)/2;
            LL=exp(1/MN*log((1-b)/epsilon));
            D=(1-b)*LL;


            for ii=1:M+1
                x(ii)=b*(ii-1)/M;
            end

            for ii=M+2:N+1
                x(ii)=1-D*exp(log(LL)*(ii-M-2));
            end
            
            this.xi=x';
            this.dxi=diff(x)';
            
            this.xi_n_plus_one=1-D*exp(log(LL)*(N+2-M-2));
            this.dxi_n=this.xi_n_plus_one-this.xi(end);
            this.N=N+1;
        end
    end
end
