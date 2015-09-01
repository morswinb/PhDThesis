classdef InitialCondition
    %InitialCondition encapsulation of simulation parameters
    %   includes various parameters that might be used
    %   a proxy to benchmark in this particular implementation
    
    properties
        a,w,w_0,V,Omega,L_0,L_inv,k,M,q_0,...
            q_l,q_l_star,q_l_norm,int_q_l,C,...
            grid,bench,alpha,beta;
    end
    
    methods
        function this=InitialCondition(grid,bench)
            
            this.a=bench.a;
            this.L_0=bench.L(0);
            this.grid=grid;
            this.alpha=bench.alpha;
            this.beta=bench.beta;
            this.w=bench.w(0,grid.xi);
            this.V=bench.V(0,grid.xi);
            this.k=bench.k;
            this.M=bench.M;
            this.q_0=bench.q_0;
            this.q_l=bench.q_l;
            this.int_q_l=bench.int_q_l;
            this.q_l_star=bench.q_l_star;
            this.q_l_norm=bench.q_l_norm;
            this.L_inv=bench.L_inv;
            this.C=bench.C;
        end
    end
    
end

