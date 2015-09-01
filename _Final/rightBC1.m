classdef rightBC1
    %rigthBC1 deals with BC at x=1
    
    properties
        A1,A2,B1,B2,xn1,xn2,N,alpha,beta;
    end
    
    methods

        function this=rightBC1(ic)
                                  
            x1=1-ic.grid.xi(end);
            x2=1-ic.grid.xi(end-1);
                       
            a=ic.alpha;
            b=ic.beta;

            this.alpha=a;
            this.beta=b;           
            this.xn1=x1;
            this.xn2=x2;

            
            this.A1=(x2^(a)-x2^(b)/x1^(b-a))^(-1);  
            this.A2=-(x2/x1)^(b)*this.A1;            
            this.B1=-this.A1/x1^(b-a);
            this.B2=1/x1^(b)-this.A2/x1^(b-a);
            this.N=ic.grid.N;
            
        end
        
        function [w_0,dwdx_n,d2wdx2_n]=get_right(this,~,~,w,~,~)
                        
            n=this.N;
            
            %finds w_0 and w_1
            A=this.A1*w(n-1)+this.A2*w(n);
            B=this.B1*w(n-1)+this.B2*w(n);
            
            a=this.alpha;
            b=this.beta;
            w_0=A;
            
            %direct result of differentiating asymptotics
            dwdx_n=-a*A*this.xn1^(a-1)+...
                -b*B*this.xn1^(b-1);
            
            %direct result of differentiating asymptotics
            d2wdx2_n=a*(a-1)*A*this.xn1^(a-2)+...
                b*(b-1)*B*this.xn1^(b-2);
        end
    end
end