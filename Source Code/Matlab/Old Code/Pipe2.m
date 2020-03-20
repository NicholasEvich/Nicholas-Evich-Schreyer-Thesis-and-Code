classdef Pipe2
    % The idea of this interpretation of the pipe idea is that there is
    % only one object, and it moves  through the pipe in an incremental
    % fashion (as opposed to one bulky pipe object that contains a large
    % ammount of smaller incremental segment objects which in turn each
    % contain a fluid proprties object... see the problem?)
    
    properties
        N
        L
        d0
        phi
        theta
        w
        dz
        angle
        Segments
        SA
    end
    
    methods
        % How can I make the declarations below simpler (like JavaScript)?
        % And, how can I declare constant properties in the constructor?
        % I'd like to make the angles constant, but I still need the
        % freedom to actually initialize them here
        function obj = Pipe2(L,d0,N,phi,theta,w)
            obj.L = L;
            obj.d0 = d0;
            obj.N = N;
            obj.phi = phi;
            obj.theta = theta;
            obj.w = w;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

