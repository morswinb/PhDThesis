classdef RegularGrid_2
    %RegularGrid_2 regularna siatka z punktem N+1 ustawionym na 1-eps 
    %(do uzycia dla "normalnych" warunków)


    properties
        N,xi,dxi,epsilon,xi_n_plus_one;       
    end
       
    methods    
        function this=RegularGrid_2(N,epsilon)
            this.N=N;
            this.epsilon=epsilon;
            this.xi=linspace(0,1-epsilon,N+1)';
            this.xi_n_plus_one=this.xi(N+1);
            this.xi=this.xi(1:N);
            this.dxi=diff(this.xi);
        end
    end
end
