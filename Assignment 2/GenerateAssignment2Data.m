% EGB242 - Assignment 2
% Generate Assignment 2 Data
%
% This file generates the assignment data needed to complete assignment 2. 
% Data is generated based on the student numbers of group members.
% This script only needs to be executed ONCE.
%
% Make sure that the following files are in the same directory.
% 1) genP1Data.p
% 2) genP2Data.p
% 3) genP3Data.p
% 4) genData.dat
%
% Student numbers have the format "n01234567".
% Omit the leading 'n' and leading '0', and enter it as 1234567.
% It should only be a 7/8 digit number.
%
% If you do not have a group of 3, enter in 0 from the missing member(s).
%
% Enter student numbers of group members below:
student1 = 8797617;
student2 = 9934731;
student3 = 9889574;


% Do not modify below this line
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
genP1Data([student1, student2, student3]);
genP2Data([student1, student2, student3]);
genP3Data([student1, student2, student3]);
