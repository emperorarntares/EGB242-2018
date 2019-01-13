%% Assignment 1 - Part A: GenerateDataAssignment1A.m
%  Clearing workspace
clear; close all; clc;
%  This file generates the assignment data needed to complete Assignment 1 - Part A. 
%  Data is generated based on your student numbers.
%  This script only needs to be executed ONCE.
%
%  Make sure that the following files are in the same directory.
%  The directory path should not have spaces nor special characters.
%  1) data
%  2) initialise1A.p
%  3) GenerateDataAssignment1A.m
%  4) preparation.m
%  5) mission.m
%
%  Enter your student number into lines 22-24.
%  Student numbers have the format "n01234567".
%  Omit the leading 'n' and leading '0', and enter it as 1234567.
%  It should only be a 7 or 8 digit number.
%  8 digit numbers will be truncated.
%  Groups of less than 3 members set the other student numbers to 0

% Student numbers
sid1 = 8797617; 
sid2 = 9934731;
sid3 = 9889574;

% Check student id length
if (length(num2str(sid1)) == 7)||(length(num2str(sid1)) == 8)
    %  Generate test signals and initialise the assignment variables.
    % Generate the noisy speech signal and 
    % the scale factor and DC shift of the noise signal before affecting the original speech signal
    [A, B, C, noiseSound] = initialise1A(sid1,sid2,sid3); % do not modify this line.
    fprintf('\nScript ran successfully.\n')
else
    fprintf('Invalid student number.\nPlease enter a valid student number.\n');
end

%  GenerateDataAssignment1A.m file will load and work with the Data1A file directly.
%  You may close this script once Data1A has been saved.
%  This script does not need to be executed again.
clear all