classdef SpectralGrid

    properties
        N,xi,dxi,;
    end
    
    methods
        
        function this=SpectralGrid(N,y)

            this.N=N;

            
            this.xi(1)=0;
            for i=1:N+1
                this.xi(i)=1-(1-(i-1)/N)^y;
            end
            
            this.xi=this.xi';
            this.dxi=diff(this.xi);
        end
    end
end



