classdef leftBC1
    %leftBC1 deals with BC at x=0
    
    properties
        q_0,k,M,x;
    end
    
    methods
        function this=leftBC1(ic)
            this.q_0=ic.q_0;
            this.k=ic.k;
            this.M=ic.M;
            %the first 3 points are extracted
            this.x=ic.grid.xi(1:3);
        end
        
        function [dwdx_0,d2wdx2_0]=get_left(this,t,L_t,w,~,~)
            
            %pumping rate q_o is proportional to first derevative
            dwdx_0=-this.M/this.k*L_t/w(1)^3*this.q_0(t);
            
            x1=this.x(1);
            x2=this.x(2);
            x3=this.x(3);
            da=-1/2./(x2-x1);
            db=1/2*(1./(x2-x1)-1/(x3-x2));
            dc=1/2./(x3-x2);
            %aproximation of a special polynomial
            d2wdx2_0=(-2/x2*da+4/x2*this.q_0(t)*...
                this.M/this.k*L_t*w(1)^(-4)-6*x2^-2)*w(1)+...
                (-2/x2*db+6*x2^(-2))*w(2)+...
                (-2/x2*dc)*w(3);
        end
    end
    
end

