classdef PipeSegment
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dz
        diam
        area
        x_in
        x_out
        alpha_in
        alpha_out
        P_abs_in
        P_abs_out
        deltaP
        v_g_in
        v_g_out
        v_f_in
        v_f_out
        h_g_in
        h_g_out
        h_f_in
        h_f_out
    end
    
    methods
        function obj = PipeSegment(dz, diam)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 1
                obj.dz = dz;
                obj.diam = diam;
                obj.area = 0.7854*diam*diam; % This is very inefficient
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

