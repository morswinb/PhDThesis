classdef ArcTanhGrid

    properties
        N,xi,dxi,epsilon,xi_n_plus_one,dxi_n;
    end
    
    methods
        
        function this=ArcTanhGrid(N,epsilon)
            this.N=N;
            this.epsilon=epsilon;

            c=atanh(1-epsilon)/N;
            for i=0:N-1
                this.xi(i+1)=tanh(c*i);
            end
            this.xi_n_plus_one=tanh(c*(N+1));
            this.xi=this.xi';
            this.dxi=diff(this.xi);
            this.dxi_n=this.xi_n_plus_one-this.xi(end);
        end
    end
end



