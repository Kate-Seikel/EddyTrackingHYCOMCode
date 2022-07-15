%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step1_extraction_SSH; Extract SLA data on an Area Of Interest
%Compute Parameters; Created Jan 27, 2019 by CT
%Updated with main function union by Paul Ernst (May 2021)
% Make a structure storing all filenames
function Step1_HYCOM(basepath, pathtodata, pathtoOcolor, years, latlonbounds, yearmetadata, fullextract, slaoradt)
%% Initialize
close all
addpath(strcat(basepath, 'FUNCTIONS'));
clc
%% Set loop for desired number of years
for i=1:length(years)
    direct = pathtodata;
    if exist(direct)==7 % 7 means name is a folder
        input_str = (pathtodata + years(i) + '/');%uigetdir(direct,'DIRECTORY OF AVISO SLA OR MADT FILES');
        input_dir = convertStringsToChars(input_str);
    else
        input_dir = uigetdir('','DIRECTORY OF AVISO SLA OR MADT FILES');
    end
    uvdir = input_dir; %uvdir='/Volumes/Lacie-SAN/SAN2/CMEMS/SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047/2019/';
    katepath = '/Volumes/Kate-Research/Data/Eddy_Extraction/';
    output_dir= strcat(katepath);
    [status,message,messageid] = mkdir(output_dir);
    NL = latlonbounds(1);
    SL = latlonbounds(2);
    WL = latlonbounds(4);
    EL = latlonbounds(3);
    Begin = yearmetadata(1);
    End = yearmetadata(2);
    extrac = fullextract;
    yearnum = str2double(years(i));
    yrb = fix(Begin/10000); datebeg = datenum(num2str(Begin),'yyyymmdd');
    yre = fix(End/10000); dateend = datenum(num2str(End),'yyyymmdd');
    if extrac ==1
        Extraction_Type = 'Basic Extraction';
    elseif extrac==2
        Extraction_Type = 'Full Extraction';
    end
    Directory_Of_NetCDF_Files = input_dir; Directory_Of_Extracted_Files = output_dir;
    Northern_Limit = NL; Southern_Limit = SL; Eastern_Limit = EL; Western_Limit = WL;
    Beginning_Date = Begin; Final_Date = End;
    Date_Of_Extraction = datestr(now);
    savedest = strcat(basepath, 'EXTRACTION/ConfigFile_fulltime.mat');
    save  (savedest, 'Date_Of_Extraction', 'Extraction_Type', 'Directory_Of_NetCDF_Files', 'Directory_Of_Extracted_Files', 'Northern_Limit', 'Southern_Limit', 'Eastern_Limit', 'Western_Limit', 'Beginning_Date', 'Final_Date');
    %% Loop through an individual year at this point
%     for mon=1
%         if mon < 10
%             monstr=['0' int2str(mon)];
%         else
%             monstr=int2str(mon);
%         end
%         yrstr=monstr;
        %yrstr
        list_SSH = [];list_UV = []; filenamesO = [];
        list_SSH = [list_SSH;dir([input_dir '/*.nc'])]; % make sure this directory exists
        list_UV = [list_UV;dir([ input_dir '/*.nc'])]; % this one too
        IND=[]; %Index

        for i=1:length(list_SSH) %5800
            name = list_SSH(i).name;
            ind = strfind(name,['4_' num2str(yrb)]);
            if isempty(ind)==0
                if datenum(name(ind+2:ind+9),'yyyymmdd')<datebeg
                    IND = [IND i];
                end
            end
            ind = strfind(name,['4_' num2str(yre)]);
            if isempty(ind)==0
                if datenum(name(ind+2:ind+9),'yyyymmdd')>dateend
                    IND = [IND i];
                end
            end
        end
        list_SSH(IND) = [];
        list_UV(IND) = [];
        clear list_SSH1 list_UV1
        for i=1:length(list_SSH) %5800
            name =  list_SSH(i).name;
            list_SSH1(i,:) = [input_dir  name]; % not calling proper directory
            name =  list_UV(i).name;
            list_UV1(i,:) = [uvdir  name]; %% DT TO NRT!!!
        end
        renvar list_SSH1 list_SSH
        renvar list_UV1 list_UV
        %Case adjustment
        if slaoradt ==0
            SSH_VAR = 'adt';
        elseif slaoradt == 1
            SSH_VAR = 'sla';
        elseif slaoradt == 2
            SSH_VAR = 'MADT';
        elseif slaoradt == 3
            SSH_VAR = 'surf_el';
        end


        %Create a waiting bar
        h = waitbar(0,['Extraction of ' SSH_VAR]);

        switch extrac
            case 1
                %Make on loop on AVISO files
                [M,N] = size(list_SSH);
                for i = 1:M %5800
                    waitbar(i/M,h)
                    filename_SSH = list_SSH(i,:);
                    filename_UV = list_UV(i,:);
                    disp(['Extracting:    ' filename_SSH ''])
                    SSH=double(ncread(filename_SSH,SSH_VAR));
                    U=double(ncread(filename_UV,'u_velocity'));
                    V=double(ncread(filename_UV,'v_velocity'));
                    % Compute on the first time parameters that do not depend of AVISO date
                    if i<365 %5800
                        % Read information needed in AVISO netcdf file
                        NbLongitudes=double(ncread(filename_SSH,'longitude'));
                        NbLatitudes=double(ncread(filename_SSH,'latitude'));
                        % Change [0,360] coordinates in [-180,180] if it is more convenient
                        % for the are of study
                        % Change as longitude are in 0 360 in AVISO FILE
    %                         if WL<0; WL=WL+360; end
    %                         if EL<0; EL=EL+360; end
    %                         if WL==EL; WL=0; EL=360;end
                        if WL<EL
                            indlon = find(NbLongitudes>=WL & NbLongitudes<=EL);
                            indlat = find(NbLatitudes>=SL & NbLatitudes<=NL);
                            X = NbLongitudes(indlon);
                            Y = NbLatitudes(indlat);
                        else
                            indlon = [find(NbLongitudes>=WL);find(NbLongitudes<=EL)];
                            indlat = find(NbLatitudes>=SL & NbLatitudes<=NL);
                            X = NbLongitudes(indlon);
                            Y = NbLatitudes(indlat);
                            X(X>180) = X(X>180)-360;
                        end
                        % Define a metadatavariable to store parameters of extraction
                        Metadata=struct('DATA','AVISO','Northern_Limit',NL,'Southern_Limit',SL,'Western_Limit',WL,'Eastern_Limit',EL,'Type_of_SSH',SSH_VAR);

                        if WL==0 && EL==360 || WL==360 && EL==0
                            ghost_point = 4;
                            ghp=ceil(ghost_point/(X(2)-X(1)));
                            dx= X(2)-X(1);
                            X2=[X(1)+(X(1)-ghp*dx:dx:(X(1)-dx)) X(1:end) X(end)+dx:dx:(X(end)+ghp*dx)];
                            renvar X2 X

                            % Write number of days used for stat
                            Metadata().IMPORTANT1 = 'MAP 0-360 in long thus periodicity satified by ghost points, use UNIQUE or center position to discretize during result'   ;
                            Metadata().ghost_point = ['Number of extra degree in each side =',num2str(ghost_point)]   ;

                        end

                        % Create matrix of coriolis parameters and define gravitationnel
                        % acceleration
                        f = repmat(coriolis(Y)',length(X),1);
                        g = 9.81;

                        % Create Matrices of distance between 2 points in x and y direction
                        delta_x=NaN(length(X),length(Y));
                        delta_y=NaN(length(X),length(Y));
                        delta_x(:,1) = ac_distance(Y(1),X(3),Y(1),X(1))*1000;
                        delta_x(:,end) =  ac_distance(Y(end),X(3),Y(end),X(1))*1000;
                        dy = Y(2)-Y(1);
                        delta_y(:,1) =  ac_distance(Y(1)+dy,X(1),Y(1)-dy,X(1))*1000;
                        dy = Y(end)-Y(end-1);
                        delta_y(:,1) =  ac_distance(Y(end)+dy,X(1),Y(end)-dy,X(1))*1000;
                        for j=2:length(Y)-1
                            delta_x(:,j) = ac_distance(Y(j), X(3),Y(j), X(1))*1000;
                            delta_y(:,j) =  ac_distance(Y(j+1),X(1),Y(j-1),X(1))*1000;
                        end
                    end

                    % Reduce SSH to area of interest
                    SSH = SSH(indlon,indlat);
                    % compute geostrophic velocity
                    U = U(indlon,indlat);
                    V = V(indlon,indlat);
                    if WL==0 && EL==360 || WL==360 && EL==0
                        SSH =cat(1,SSH(end-ghp+1:end,:),SSH,SSH(1:ghp,:));
                        U =cat(1,U(end-ghp+1:end,:),U,U(1:ghp,:));
                        V =cat(1,V(end-ghp+1:end,:),V,V(1:ghp,:));
                    end

                    %Norm of speed and EKE
                    EKE=(U.^2+V.^2)./2;

                    % Take dtime values and store it in a datenum format
                    if yearnum < 2021
                        date_num=datenum(filename_SSH(end-19:end-12),'yyyymmdd');
                        % Store variables in a .mat file
                        filename_out = [output_dir SSH_VAR '_uv_' filename_SSH(78:85)];
                    else
                        date_num=datenum(filename_SSH(end-10:end-3),'yyyymmdd');
                        % Store variables in a .mat file
                        filename_out = [output_dir SSH_VAR '_uv_' filename_SSH(end-10:end-3)];
                    end


                    U = single(U);
                    V = single(V);
                    EKE = single(EKE);
                    SSH = single(SSH);

                    if strcmp(SSH_VAR,'adt')  %'adt')
                        renvar SSH ADT % SLA
                        save(filename_out,'Metadata','X','Y','date_num','U','V','EKE','ADT');
                    elseif strcmp(SSH_VAR,'MADT')
                        renvar SSH MADT
                        save(filename_out,'Metadata','X','Y','date_num','U','V','EKE','MADT');
                    elseif strcmp(SSH_VAR, 'sla')
                        renvar SSH SLA
                        save(filename_out,'Metadata','X','Y','date_num','U','V','EKE','SLA');
                    elseif strcmp(SSH_VAR, 'surf_el')
                        renvar SSH surf_el
                        save(filename_out,'Metadata','X','Y','date_num','U','V','EKE','surf_el');
                    end
                end
            close(h);

        case 2
                    %Make on loop on AVISO files
                    [M,N] = size(list_SSH);
                    for i = 1:M
                        waitbar(i/M,h)
                        % Display file treated name and create associated filename
                        filename_SSH = list_SSH(i,:);
                        filename_UV = list_UV(i,:);
                        disp(['Extracting:    ' filename_SSH ''])
                        %gunzip(filename_SSH)
                        %filename_SSH = filename_SSH(1:end-3);
                        %gunzip(filename_UV)
                        %filename_UV = filename_UV(1:end-3);
                        SSH=double(ncread(filename_SSH,SSH_VAR));
                        U=double(ncread(filename_UV,'u'));
                        V=double(ncread(filename_UV,'v'));
                        % Compute on the first time parameters that do not depend of AVISO date
                        if i==1
                            % Read information needed in AVISO netcdf file
                            NbLongitudes=double(ncread(filename_SSH,'lon'));
                            NbLatitudes=double(ncread(filename_SSH,'lat'));
                            % Change [0,360] coordinates in [-180,180] if it is more convenient
                            % for the are of study
                            % Change as longitude are in 0 360 in AVISO FILE
%                                 if WL<0; WL=WL+360; end
%                                 if EL<0; EL=EL+360; end
%                                 if WL==EL; WL=0; EL=360;end
                            if WL<EL
                                indlon = find(NbLongitudes>=WL & NbLongitudes<=EL);
                                indlat = find(NbLatitudes>=SL & NbLatitudes<=NL);
                                X = NbLongitudes(indlon);
                                Y = NbLatitudes(indlat);
                            else
                                indlon = [find(NbLongitudes>=WL);find(NbLongitudes<=EL)];
                                indlat = find(NbLatitudes>=SL & NbLatitudes<=NL);
                                X = NbLongitudes(indlon);
                                Y = NbLatitudes(indlat);
                                X(X>180) = X(X>180)-360;
                            end

                            % Define a metadatavariable to store parameters of extraction
                            Metadata=struct('DATA','AVISO','Northern_Limit',NL,'Southern_Limit',SL,'Western_Limit',WL,'Eastern_Limit',EL,'Type_of_SSH',SSH_VAR,'Ghost_points',4,'IMPORTANT1','MAP 0-360 in long thus periodicity satified by ghost points, use UNIQUE or center position to discretize during result'  );

                            if WL==0 && EL==360 || WL==360 && EL==0
                                ghost_point = 4;
                                ghp=ceil(ghost_point/(X(2)-X(1)));
                                dx= X(2)-X(1);
                                X2=[X(1)+(X(1)-ghp*dx:dx:(X(1)-dx)) X(1:end)' X(end)+dx:dx:(X(end)+ghp*dx)];
                                renvar X2 X

                                % Write number of days used for stat
                                %Metadata().IMPORTANT1 = 'MAP 0-360 in long thus periodicity satified by ghost points, use UNIQUE or center position to discretize during result'   ;
                                %Metadata().ghost_point = ['Number of extra degree in each side =',num2str(ghost_point)]   ;

                            end

                            % Create matrix of coriolis parameters and define gravitationnel
                            % acceleration
                            f = repmat(coriolis(Y)',length(X),1);
                            g = 9.81;

                            % Create Matrices of distance between 2 points in x and y direction
                            delta_x=NaN(length(X),length(Y));
                            delta_y=NaN(length(X),length(Y));
                            delta_x(:,1) = ac_distance(Y(1),X(3),Y(1),X(1))*1000;
                            delta_x(:,end) =  ac_distance(Y(end),X(3),Y(end),X(1))*1000;
                            dy = Y(2)-Y(1);
                            delta_y(:,1) =  ac_distance(Y(1)+dy,X(1),Y(1)-dy,X(1))*1000;
                            dy = Y(end)-Y(end-1);
                            delta_y(:,1) =  ac_distance(Y(end)+dy,X(1),Y(end)-dy,X(1))*1000;
                            for j=2:length(Y)-1
                                delta_x(:,j) = ac_distance(Y(j), X(3),Y(j), X(1))*1000;
                                delta_y(:,j) =  ac_distance(Y(j+1),X(1),Y(j-1),X(1))*1000;
                            end
                        end


                        % Reduce SSH to area of interest
                        SSH = SSH(indlon,indlat);
                        % compute geostrophic velocity
                        U = U(indlon,indlat);
                        V = V(indlon,indlat);
                        if WL==0 && EL==360 || WL==360 && EL==0
                            SSH =cat(1,SSH(end-ghp+1:end,:),SSH,SSH(1:ghp,:));
                            U =cat(1,U(end-ghp+1:end,:),U,U(1:ghp,:));
                            V =cat(1,V(end-ghp+1:end,:),V,V(1:ghp,:));
                        end

                        % compute Derivative of veloticy

                        % calcul de dU/dx
                        dU = NaN*U;
                        dU(2:end-1,:) = (U(3:end,:)-U(1:end-2,:));
                        dX = repmat(delta_x(1,:),[size(U,1) 1]);
                        Ux = dU./dX;

                        %calcul de dV/dx
                        dV = NaN*V;
                        dV(2:end-1,:) = (V(3:end,:)-V(1:end-2,:));
                        Vx = dV./dX;

                        % calcul de dU/dy
                        dU = NaN*U;
                        dU(:,2:end-1) = (U(:,3:end)-U(:,1:end-2));
                        dy = (Y(3)-Y(1))*ac_distance(0,0,1,0)*1000;
                        dY = repmat(dy,[size(U,1) size(U,2)]);
                        Uy = dU./dY;

                        % calcul de dV/dy
                        dV = NaN*V;
                        dV(:,2:end-1) = (V(:,3:end)-V(:,1:end-2));
                        Vy = dV./dY;

                        % Vorticity
                        Vorticity=Vx-Uy;

                        % Shear, normal componant of strain and Okubo-Weiss parameter
                        Ss=Vx+Uy;
                        Sn=Ux-Vy;
                        OW=Sn.^2+Ss.^2-Vorticity.^2;

                        %Norm of speed and EKE
                        Speed=sqrt(U.^2+V.^2);
                        EKE=(U.^2+V.^2)./2;

                        % Take dtime values and store it in a datenum format
                        date_num=datenum(filename_SSH(end-19:end-12),'yyyymmdd');


                        % Store variables in a .mat file
                        filename_out = [output_dir 'EXTRACTION/INDIVIDUAL_FILES/' SSH_VAR '_uv_' filename_SSH(78:85)];

                        U = single(U);
                        V = single(V);
                        EKE = single(EKE);
                        SSH = single(SSH);
                        Vorticity = single(Vorticity);
                        Ss = single(Ss);
                        Sn = single(Sn);
                        OW = single(OW);
                        Speed = single(Speed);
                        if strcmp(SSH_VAR,'sla')
                            renvar SSH SLA
                            save(filename_out,'Metadata','X','Y','date_num','U','V','Vorticity','Ss','Sn','OW','Speed','EKE','SLA');
                        elseif strcmp(SSH_VAR,'MADT')
                            renvar SSH MADT
                            save(filename_out,'Metadata','X','Y','date_num','U','V','Vorticity','Ss','Sn','OW','Speed','EKE','MADT');
                        elseif strcmp(SSH_VAR,'surf_el')
                            renvar SSH surf_el
                            save(filename_out,'Metadata','X','Y','date_num','U','V','Vorticity','Ss','Sn','OW','Speed','EKE','surf_el');
                        end
                        %delete(filename_SSH);
                        %delete(filename_UV);
                    end
                close(h);
            end
end
