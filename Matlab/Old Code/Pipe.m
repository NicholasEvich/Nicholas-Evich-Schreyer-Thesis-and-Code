classdef Pipe
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N
        L
        d0
        dz
        angle
        Segments
        SA
    end
    
    methods
        function obj = Pipe(N, L, d0, angle)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.N = N;
            obj.L = L;
            obj.d0 = d0;
            obj.dz = L/N;
            obj.angle = angle;
            % Refer to N or obj.N for the rest of this function?
            Segments(1,250) = PipeSegment(250);
            obj.Segments = Segments;
            obj.SA = pi*d0*L + 0.5*pi*L*L*tand(angle);
            obj.Segments = obj.initSegments(N);
        end
        
        % Two ways to do this: the way it is, or make the constructor of
        % the pipe segment class do literally nothing: then call the
        % diameter and surface area calculaton methods in this function
        function Segments = initSegments(obj, N)            
            for i = 1:N/2
               z = obj.L*i/N;
               d = obj.d0 + 2*z*tand(obj.angle);
               obj.Segments(i) = PipeSegment(obj.dz, d);
            end
            
            for i = N/2:N
                z = obj.L*i/N;
                d = obj.d0 + 2*obj.L*tand(obj.angle) - 2*z*tand(obj.angle);
                obj.Segments(i) = PipeSegment(obj.dz, d);
            end
            
            Segments = obj.Segments;
        end
        
    end
end

