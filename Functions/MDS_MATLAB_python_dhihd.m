% MDS_MATLAB_python.m

% Author:  Michael Ross
% Description:  This code calls a Python script to extract data from
% MDSplus and save that data as a CSV file.  Then, this script reads in the
% data from the CSV file for use in MATLAB.  This script enables the user
% to avoid the memory leakage issues MATLAB experiences when extracting
% data directly from MDSplus data trees.

function [signal] = MDS_MATLAB_python_dhihd(shotnum,signame)

% Call Python script in Windows DOS to extract data from MDSplus:
status = dos(['python M:\Users\Michael_Ross\Coding\Functions\MDS_extract_dhihd.py ',num2str(shotnum),' ',signame]);

% Read in data from CSV file stored by Python script:
data_mat = csvread('M:\Users\Michael_Ross\Coding\Functions\MDS_temp_dhihd.csv');

% Store data and export out of this MATLAB function:
signal = data_mat;
% signal = data_mat(:,2);
