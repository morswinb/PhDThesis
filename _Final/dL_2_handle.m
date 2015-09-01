classdef dL_2_handle
    %dL_2_handle liczenie L po jako L^2 i pochodna z L^2 
    
    properties
        L_0,k,M;
    end
    
    methods
        function this=dL_2_handle(initial_condition)
            this.L_0=initial_condition.L_0;
            this.k=initial_condition.k;
            this.M=initial_condition.M;
        end
        
        function L=get_L(this,w)
            L=w(end)^.5;
        end
        
        function dL=get_dL(this,w_0,w)
            dL=2*this.k/this.M/3*w_0^3;
        end
        
        function L_0=get_L_0(this)
            L_0=this.L_0^2;
        end
    end
    
end

