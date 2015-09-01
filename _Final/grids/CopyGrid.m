classdef CopyGrid

    properties
        N,xi,dxi,epsilon,xi_n_plus_one,dxi_n;
    end
    
    methods
        
        function this=CopyGrid(oldGrid,N)

            this.N=N;
            this.epsilon=oldGrid.xi(N+1)-oldGrid.xi(N);
            
            this.xi=oldGrid.xi(1:N);
            
            a=(1-this.epsilon)/this.xi(this.N);
            this.xi=this.xi*a;
            
            
%             this.epsilon=1-this.xi(this.N);
            this.dxi=diff(this.xi);

        end
    end
end


