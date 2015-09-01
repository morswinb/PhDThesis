classdef QuadraticGrid

    properties
        N,xi,dxi,epsilon,xi_n_plus_one,dxi_n;
    end
    
    methods
        
        function this=QuadraticGrid(N,epsilon)

            this.N=N;
            this.epsilon=epsilon;

            y=2;
            
            this.xi(1)=0;
            for i=1:(N-1)
                this.xi(i+1)=1-(1-(1-epsilon^(1/y))*i/(N-1))^y;
            end
            
            this.xi_n_plus_one=(1+this.xi(end))/2;
            this.xi=this.xi';
            this.dxi=diff(this.xi);
            this.dxi_n=this.xi_n_plus_one-this.xi(end);

        end
    end
end



