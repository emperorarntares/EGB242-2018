%% Assignment 1 - Part B: GenerateDataAssignment1B.m
%  Clearing workspace
clear all; close all; clc;
%  This file generates the assignment data needed to complete Assignment 1 - Part B. 
%  Data is generated based on individual student numbers.
%  This script only needs to be executed ONCE.
%
%  Make sure that the following files are in the same directory.
%  The directory path should not have spaces nor special characters.
%  1) initialise1B.p
%  2) A1BData.bin
%  3) A1BLPF.p
%  4) A1BTextdecode.p
%  5) training.m
%  6) mission.m
%
%  Enter your student numbers into line 21.
%  Student numbers have the format "n01234567".
%  Omit the leading 'n' and leading '0', and enter it as 1234567.
%  It should only be a 7 or 8 digit number.
st = 9934731; % insert student number here

%  Initialise the assignment variables and save them into Data1B.mat.
initialise1B(st); % do not modify this line.

%  The initialise.p function will automatically save the variables.
%  A message will be displayed on the Command Window upon successful 
%  execution of the code.
%
%  The training.m file will load and work with the Data1B.mat file directly.
%  You may close this script once Data1B.mat has been saved.

%  This script only needs to be executed ONCE!!