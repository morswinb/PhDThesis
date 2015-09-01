classdef QuadraticGrid_2

    properties
        N,xi,dxi,epsilon;
    end
    
    methods
        
        function this=QuadraticGrid_2(N,epsilon,y)
            this.N=N;
            this.epsilon=epsilon;
            for i=1:(N-1)
                this.xi(i+1)=1-(1-(1-epsilon^(1/y))*i/(N-1))^y;
            end
            this.xi=this.xi';
            this.dxi=diff(this.xi);
        end
    end
end



