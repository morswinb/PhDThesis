classdef dL_handle
    %dL_handle liczenie L i pochodnej po L
    properties
        L_0,k,M;
    end
    
    methods
        function this=dL_handle(initial_condition)
            this.L_0=initial_condition.L_0;
            this.k=initial_condition.k;
            this.M=initial_condition.M;
        end
        
        function L=get_L(this,w)
            L=w(end);
        end
        
        function dL=get_dL(this,w_0,w)
            dL=this.k/this.M/3/w(end)*w_0^3;
        end
        
        function L_0=get_L_0(this)
            L_0=this.L_0;
        end
    end
    
end

