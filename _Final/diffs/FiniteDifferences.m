classdef FiniteDifferences<handle
    %handle for computing Crack derevative
    %central FD scheme with variable interval
    %length is used
    
    properties
        D,F,G,N,grid,dy;
    end
    
    methods
        function this=FiniteDifferences(grid)
            this.N=grid.N;
            n=this.N-1;
            dxi=grid.dxi;
            
            %precompute constant parameters
            this.D=dxi(1:n-1);
            this.F=dxi(2:n);    
            this.G=1./dxi(1:n-1)+1./dxi(2:n);
        end
        
        function preset(this,y)
            n=this.N;
            this.dy=diff(y(1:n));  
        end
        
        function fstder=calcFirstDer(this)
            n=this.N;
            fstder=zeros(n,1);
            fstder(2:n-1)=1/2*(this.dy(2:n-1)...
                ./this.F+this.dy(1:n-2)./this.D);
        end
        
        function secder=calcSecondDer(this)
            n=this.N;
            secder=zeros(n,1);
            secder(2:n-1)=1/2*this.G.*(this.dy(2:n-1)...
                ./this.F-this.dy(1:n-2)./this.D);
        end
    end    
end

