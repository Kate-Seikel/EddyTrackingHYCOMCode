%% StepsMain.m (version 1.0)
%Author: Paul Ernst
%Date Created: 5/10/2021
%Date of Last Update: 5/10/2021
%What was last update?
%Created.
%--------------------------------------
%Purpose: Main function for the steps of the Eddy Tracking code. 
%Inputs: None.
%Outputs: We'll figure that one out together.
%% Inputs
basepath = '/Volumes/NEMO3/GOFFISH/Emily/Programs/Eddy_Tracking_HYCOM_NEW/'; %this is the path where all of the step files sit: the directory you're currently running this on
pathtodata = '/Volumes/LACIE-GOM/GOM-HYCOM-1993-2020-Daliy3D-Ebenezer/daily/'; %this is the path to all of the data you'll be crunching initially
% We're assuming that the paths to functions are in a subfolder ("FUNCTIONS")
% We're also assuming that you have a file tree of EXTRACTION -> INDIVIDUAL
% FILES / EDDY PROPERTIES / EDDY TRAJECTORIES set up under basepath
% The above pathnames MUST BE USING SINGLE QUOTES and end with a /
% We're assuming that under your pathtodata folder, there are subfolders
% labelled with the years you list in the below variable (i.e. "2020")
years = ["2019","2020","2021"];%enter all years
% % years = ["2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"]; %enter all years formatted like ["2019", "2020"]-- DOUBLE QUOTES ESSENTIAL
latlonbounds = [31, 18, -80, -98]; % [N, S, E, W] lat long boundaries
yearmetadata = [20190101, 20211231]; %[YYYYMMDD, YYYYMMDD] start and end dates of the data
fullextract = 1; %input('Basic (1: EKE) or Full Extraction (2: velocity, vorticity, OW components, etc.);
slaoradt = 3; %0 = SLA, 1 = ADT, 2 = MADT, 3 = surf_el set for your preferred extracted variable
minlife = 0; %minimum life of eddy to be processed
minamp = 0; %minimum amplitude of eddy to be processed
minrad = 0; %minimum radius of eddy to be processed
cut=(1:365); %These are the days you will use in your Step 5 analysis
               %(152:304) covers June-September (so think days from Jan. 1)
addpath(basepath);
%% Run Associated Steps.
%% Step 0.
%This is a reminder to check if you have a C Compiler, and if you've
    %recreated the inpoly.mexmaci64 file for this computer.
    %If you haven't recreated it for the most recent version of matlab, run
    %mex inpoly.c and place the resulting file in FUNCTIONS. 
    %Theoretically, the following lines should do it for you:
funcpath = strcat('/Users/Kate/Downloads/inpoly.c');
mex -setup;
mex(funcpath);
movepath = strcat(basepath, 'FUNCTIONS');
[status,message,messageId] = movefile('inpoly.mexmaci64', movepath);
%% Step 1.
Step1_HYCOM(basepath, pathtodata, years, latlonbounds, yearmetadata, fullextract, slaoradt);
%% Step 2.
Step2_Identification_eddies_HYCOM(basepath, minamp, minrad, fullextract);
%% Step 3.
Step3_Merge_eddymat_HYCOM(basepath, fullextract);
%% Step 4.
Step4_Stat_trajectories_HYCOM(basepath);
%% Step 5.
Step5_Trajectories_HYCOM(basepath,cut,fullextract);
%% Step 6.
Step6_Filtered_Trajectories_HYCOM(minlife,minamp,minrad,basepath);
