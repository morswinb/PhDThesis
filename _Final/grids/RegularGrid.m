classdef RegularGrid
    %uniform grid
    properties
        N,xi,dxi,epsilon;       
    end
       
    methods    
        function this=RegularGrid(N,epsilon)
            this.N=N;
            this.epsilon=epsilon;
            this.xi=linspace(0,1-epsilon,N)';        
            this.dxi=diff(this.xi);
        end
    end
end
