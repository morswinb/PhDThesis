classdef DiffSpline<handle
    %handle for computing Crack derevative
    %spline is constructed by spline() function
    %and subsequently its derevative given by fnder
    %is used
    
    properties
        N,x,pp;
    end
    
    methods
        function this=DiffSpline(grid)       
            this.x=grid.xi;
            this.N=grid.N;
        end
        
        function preset(this,y)
            n=this.N;
            %spline is constructed
            this.pp=spline(this.x,y(1:n));
        end
        
        function fstder=calcFirstDer(this)
            n=this.N;
            fstder=zeros(n,1);
            %function for first derevative
            f=@(x)ppval(fnder(this.pp,1),x);
            %evaluation over grid x
            fstder(2:n-1)=f(this.x(2:n-1));               
        end
        
        function secder=calcSecondDer(this)
            n=this.N;
            secder=zeros(n,1);
            %function for second derevative
            f=@(x)ppval(fnder(this.pp,2),x);
            %evaluation over grid x
            secder(2:n-1)=f(this.x(2:n-1));     
        end
    end
    
end

