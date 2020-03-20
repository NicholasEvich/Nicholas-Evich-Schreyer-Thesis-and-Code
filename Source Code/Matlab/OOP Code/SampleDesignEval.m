% Nicholas Evich
% 2/11/2020
% Sample script for creating and evaluating a channel design

clear;clc

L = 1; % channel length, in meters
N = 250; % number of segments to break the channel into
w = 0.5; % mass flow rate in kg/s
Q = 500000; % total heat input (W)

% r = @(z) 0.1 + 0.01*sin(0.5*pi*z); % add verification to make sure the
% surfaces do not overlap

% Note that for this test case, all three functions below are
% constant-valued
r = @(z) 0.05; % function for channel radius vs axial distance
centerline = @(z) 0; % function for centerline height vs axial distance
qflux = @(z) Q/(L*pi*0.1); % function for heat flux vs axial distance

A = Pipe(L,N,w,r,centerline,qflux); % init pipe with no incline

% A is an object/instance of class "Pipe"
T0 = 373.15; P0 = 101300; x0 = 0; % entrance temp, pressure, flow quality
% initFluidProps is a member function of the 
% Pipe class. This function can be used as an operation on any Pipe object,
% and it initializes a data member of class Water. The Water class inherits
% from both the HEM class and the HeatTransfer classes, the two classes
% that determine channel pressure drops and heat transfer performance
A.initFluidProps(T0,P0,x0); 

[deltaP, h_overall, xout] = A.evalDesign();