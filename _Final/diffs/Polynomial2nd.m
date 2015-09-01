classdef Polynomial2nd<handle
    %handle for computing Crack derevative
    %approximation where derevative at x_i is given by 2ax+b and 2a 
    %as approximated by ax^2+bx+c on points x_i-1,x_i,x_i+1
    properties
        N,grid,A,B,C,E,F,dy;
    end
    
    methods
        function this=Polynomial2nd(grid)       
            this.grid=grid;
            this.N=grid.N;
            
            %precompute constant parameters
            n=this.N;
            dxi=grid.dxi;
            this.A=1./(dxi(2:n-1)+dxi(1:n-2));
            this.B=dxi(1:n-2)./dxi(2:n-1);
            this.C=dxi(2:n-1)./dxi(1:n-2);
            this.E=1./dxi(1:n-2);
            this.F=1./dxi(2:n-1); 
        end
        
        function preset(this,y)
            n=this.N;
            this.dy=diff(y(1:n));  
        end
        
        function fstder=calcFirstDer(this)
            n=this.N;
            fstder=zeros(n,1);
            fstder(2:n-1)=this.A.*(this.dy(2:n-1)...
                .*this.B+this.dy(1:n-2).*this.C);          
        end
        
        function secder=calcSecondDer(this)
            n=this.N;
            secder=zeros(n,1);
            secder(2:n-1)=2*this.A.*(this.dy(2:n-1)...
                .*this.F-this.dy(1:n-2).*this.E);
        end
    end  
end

