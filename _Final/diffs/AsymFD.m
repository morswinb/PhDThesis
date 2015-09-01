classdef AsymFD<handle
    %handle for computing Crack derevative 
    %using 2 terms asymptotics.
    %Asymptotic approximation is used to calculate most 
    %of the value, thous better accuracy is expected. 
    %The reminder is treated with central FD
    
    properties
        grid,N,A1,A2,B1,B2,xn1,xn2,alpha,beta,dy,f,fx,fxx,...
            A,B,C,D,F,G,xa,xb,xaa,xbb,xaaa,xbbb;
    end
    
    methods
        %constructor grid, first asym power, second asym power
        function this=AsymFD(grid,ic)  
            %assign grid, N
            this.N=grid.N;
            this.grid=grid;
            a=ic.alpha; % first asym term power (usually 1/3)
            b=ic.beta; % second asym term power
            this.alpha=a;
            this.beta=b;
            %working tip parameters for finding asym terms
            x1=1-grid.xi(end);
            x2=1-grid.xi(end-1);
            this.A1=(x2^(a)-x2^(b)/x1^(b-a))^(-1);  
            this.A2=-(x2/x1)^(b)*this.A1;            
            this.B1=-this.A1/x1^(b-a);
            this.B2=1/x1^(b)-this.A2/x1^(b-a);
            this.xn1=x1;
            this.xn2=x2;
            
            %precomputing constant asym parameters
            x=grid.xi;
            this.xa=(1-x).^a;
            this.xb=(1-x).^b;
            this.xaa=-a*(1-x).^(a-1);
            this.xbb=-b*(1-x).^(b-1);
            this.xaaa=a*(a-1)*(1-x).^(a-2);
            this.xbbb=b*(b-1)*(1-x).^(b-2);
            
            %precomputing constant FD parameters
            dxi=grid.dxi;
            n=this.N-1;
            this.A=1./(dxi(2:n)+dxi(1:n-1));
            this.B=dxi(1:n-1)./dxi(2:n);
            this.C=dxi(2:n)./dxi(1:n-1);
            this.D=dxi(1:n-1);
            this.F=dxi(2:n);    
            this.G=1./dxi(1:n-1)+1./dxi(2:n);
        end
        
        function preset(this,y)
            n=this.N;
            %finds asym terms on 2 last data points
            a=this.A1*y(n-1)+this.A2*y(n);
            b=this.B1*y(n-1)+this.B2*y(n); 
            
            %calculates asym based approximation for
            %F dFdx and d2Fdx2
            this.f=a*this.xa+b*this.xb;
            this.fx=a*this.xaa+b*this.xbb;              
            this.fxx=a*this.xaaa+b*this.xbbb;
            
            %difference of asym approximation and acctual value
            this.dy=diff(y(1:n)-this.f);  
        end
        
        function fstder=calcFirstDer(this)
            n=this.N;
            fstder=zeros(n,1);
            fstder(2:n-1)=this.fx(2:n-1)... %asym derivative term
                +1/2*(this.dy(2:n-1)./this.F...%FD derivative term
                +this.dy(1:n-2)./this.D);
        end
         
        function secder=calcSecondDer(this)
             n=this.N;
             secder=zeros(n,1); 
             secder(2:n-1)=this.fxx(2:n-1)...%asym derivative term
                 +1/2*this.G.*...            %FD derivative term
                 (this.dy(2:n-1)./this.F-this.dy(1:n-2)./this.D);
        end
    end   
end

