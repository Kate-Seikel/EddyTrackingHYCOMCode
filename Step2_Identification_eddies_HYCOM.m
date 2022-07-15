%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       CETTE FONCTION PERMET D' IDENTIFIER LES TOURBILLONS           %%%
%%%       EN UTILISANT LA TECHNIQUE DES CONTOURS FERME DE SSH           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function Step2_Identification_eddies(input_dir,min_amp,minDist_deg)
%Identify each eddies using SSH contour from daily matrices of SSH where
%several parameters have been computed (ex EKE, Vorticity, etc)
%Last modification 06/05/14
function Step2_Identification_eddies_BoB(basepath, minamp, minrad, fullextract)
addpath(strcat(basepath, 'FUNCTIONS'));
clc
katepath = '/Volumes/Kate-Research/Data/Eddy_Extraction/';
loadfile = strcat(basepath, 'EXTRACTION/ConfigFile_fulltime.mat');
load(loadfile);  %ConfigFile2_fulltime.mat');
input_dir = Directory_Of_Extracted_Files;
list = dir([ input_dir '/' '*.mat']);
filename = [input_dir list(1).name];
load(filename,'X','Y')
disp('')
min_amp = minamp/100;
disp('')

Rmin = ceil(2*max([diff(unique(X(:)));diff(unique(Y(:)))])*111.1/pi);
disp(['Considering the actual resolution and in order to include']);
disp(['at least 4 grid points inside the detected eddies,'])
disp(['the minimum eddy radius will be of ' num2str(Rmin) ' km']);
disp('')
min_R = max([minrad Rmin]);


disp('')
minDist_deg = 2;
disp('')
Minimum_Amplitude = min_amp;
Minimum_Radius = min_R;

savepath = strcat(basepath, 'EXTRACTION/Configfile_fulltime.mat');

save(savepath, 'Minimum_Amplitude', 'Minimum_Radius', '-append'); %Configfile2_fulltime.mat Minimum_Amplitude Minimum_Radius
if fullextract == 1
    Extraction_Type ='Basic Extraction';
elseif fullextract == 2
    Extraction_Type = 'Full Extraction';
end
h = waitbar(0,'Eddy Extraction');

switch Extraction_Type
    case 'Basic Extraction'
        
        
        Fields = [
            '1.Longitude of the Eddy Center          '; ...
            '2.Latitude of the Eddy Center           '; ...
            '3.Longitude of the Eddy Centroid        '; ...
            '4.Latitude of the Eddy Centroid         '; ...
            '5.Longitudes of the Eddy Edge           '; ...
            '6.Latitudes of the Eddy Edge            '; ...
            '7.Equivalent Radius [km]                '; ...
            '8.Eddy Area [km^2]                      '; ...
            '9.Amplitude [m]                         '; ...
            '10.Mean EKE [(m/s)^2]                   ';];
        
        
        % Loop on each AVISO file (thus on each dates)
        for i=[1:length(list)]
            %for i=[1:1096]
            waitbar(i/length(list),h)
            
            % Display file treated name and create associated filename
            disp(['Extracting eddies for:    ' list(i).name ''])
            filename = [input_dir  list(i).name];
            load(filename)
            EKE = double(EKE);
            if exist('ADT','var')==1 %'SLA','var')==1
                ADT = double(ADT);
                %SLA = double(SLA);
            elseif exist('MADT','var')==1
                MADT = double(MADT);
            elseif exist('SLA','var') ==1
                SLA = double(SLA);
            end
            U = double(U);
            V = double(V);
            X = double(X);
            Y = double(Y);
            
            type_SSH=Metadata.Type_of_SSH;
            
            Label_anti = single(zeros(length(X),length(Y)));
            Label_cyclo = single(zeros(length(X),length(Y)));
            
            if strcmp(type_SSH,'adt')%'sla')
                SSH=ADT;
                %SSH=SLA;
                clear ADT
                %clear SLA
            elseif strcmp(type_SSH,'madt')
                SSH=MADT;
                clear MADT
            elseif strcmp(type_SSH, 'sla')
                SSH=SLA;
                clear SLA
            elseif strcmp(type_SSH, 'surf_el')
                SSH=surf_el;
                clear surf_el
            end
            
            % Compute the Windows size in degree to number of point during the
            % first iteration and create domain size matrices of coordinate
            if i==1
                minDist=ceil(minDist_deg/(X(2)-X(1)));
                disp(['The windows size is ' , num2str(minDist*(X(2)-X(1))), ' degrees, corresponding to ',num2str(minDist), ' pixels '])
                clear minDist_deg;
                [y,x] = meshgrid(Y,X);
            end
            SSH=double(SSH);
            % Call Function that extract eddy extraction from closed countours
            % technic. Return position of centers, values of amplitude, contour
            % position and number of eddies
            [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo] = extraction_eddies_new(X,Y,SSH,minDist,min_amp,min_R,1);
            
            
            % Compute and store parameters of each Cyclones if there is at least
            % one detected
            if Ncyclo>0
                for j=1:Ncyclo %Loop on each Cyclone
                    % Edge coordinates
                    disp(['(' int2str(i) ')    ' int2str(j) '   of   ' int2str(Ncyclo)])
                    Xcon = double(Cyclonic_Cell{j,5});
                    Ycon = double(Cyclonic_Cell{j,6});
                    
                    % Change x & y to reduce search area
                    ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);
                    
                    % Coordinate index of grid points inside the eddy
                    in = find(inpoly([x(ind1),y(ind1)]',[Xcon(:),Ycon(:)]')==1);
                    INDinside = ind1(in);
                    
                    % Xcenter and Ycenter
                    [~,iiind] = min(SSH(INDinside));
                    INDcenter = INDinside(iiind(1));
                    
                    Cyclonic_Cell{j,1}= single(x(INDcenter));
                    Cyclonic_Cell{j,2}= single(y(INDcenter));
                    
                    % Mean EKE
                    Cyclonic_Cell{j,10}=single(nanmean(EKE(INDinside)));
                    
                    Label_cyclo(ind1(in)) = j;
                end
            end
            
            % Compute and store parameters of each Anticyclones if there is at least
            % one detected
            if Nanti>0
                for j=1:Nanti %Loop on each Cyclone
                    disp(['(' int2str(i) ')    '  'Cyclone loop:   ' int2str(j)  ' of '  int2str(Nanti)]);
                    % Edge coordinates
                    Xcon = double(Anticyclonic_Cell{j,5});
                    Ycon = double(Anticyclonic_Cell{j,6});
                    
                    % Change x & y to reduce search area
                    ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);
                    
                    % Coordinate index of grid points inside the eddy
                    a1=[x(ind1),y(ind1)]';b1=[Xcon(:),Ycon(:)]';%pause(20);
                    in = find(inpoly(a1,b1)==1);%pause(5)
                    INDinside = ind1(in);
                    
                    % Xcenter and Ycenter
                    [~,iiind] = max(SSH(INDinside));
                    INDcenter = INDinside(iiind(1));
                    
                    Anticyclonic_Cell{j,1}= single(x(INDcenter));
                    Anticyclonic_Cell{j,2}= single(y(INDcenter));
                    
                    % Mean EKE
                    Anticyclonic_Cell{j,10}=single(nanmean(EKE(INDinside)));
                    
                    
                    Label_anti(ind1(in)) = j;
                end
            end
            
            % Save Metada parameters
            %Metadata().Min_Amplitude_cm = min_amp;
            %Metadata().Windows_size_pixel = minDist;
            % Metada().Min_Radius_km = min_R;
            
            % Saves cells and variables in .mat file
            eval([' save -append ' filename ' Metadata Anticyclonic_Cell Cyclonic_Cell Nanti Ncyclo Fields Label_cyclo Label_anti'])
            
        end
        close(h)
        
    case 'Full Extraction'
        
        % Create Field of cell variable name and fille Metadata structure
        Fields = [
            '1.Longitude of the Eddy Center          '; ...
            '2.Latitude of the Eddy Center           '; ...
            '3.Longitude of the Eddy Centroid        '; ...
            '4.Latitude of the Eddy Centroi          '; ...
            '5.Longitudes of the Eddy Edge           '; ...
            '6.Latitudes of the Eddy Edge            '; ...
            '7.Equivalent Radius [km]                '; ...
            '8.Eddy Area [km^2]                      '; ...
            '9.Amplitude [m]                         '; ...
            '10.Mean EKE [(m/s)^2]                   '; ...
            '11.Mean Speed [m/s]                     '; ...
            '12.Mean Vorticity [1/s]                 '; ...
            '13.Vorticity at center  [1/s]           '; ...
            '14.Normalized Vorticity at center [1/s] '; ...
            '15.Mean Straining deformation rate [1/s]'; ...
            '16.Mean shearing deformation rate [1/s] '; ...
            '17.Mean Okubo_Weiss parameter  [1/s]    ';...
            '18.Longitudes inside the Eddy           ';...
            '19.Latitudes inside the Eddy            ';...
            '20.U inside the Eddy                    ';...
            '21.V inside the Eddy                    '];
        
        
        
        % Loop on each AVISO file (thus on each dates)
        %         for i=1:length(list)
        %             waitbar(i/length(list),h)
        for i=1:366
            waitbar(i/366,h)
            
            % Display file treated name and create associated filename
            disp(['Extracting eddies for:    ' list(i).name ''])
            filename = [input_dir '/EXTRACTION/INDIVIDUAL_FILES/' list(i).name];
            load(filename)
            EKE = double(EKE);
            if exist('ADT','var')==1 %'SLA','var')==1
                ADT = double(ADT);
                %SLA = double(SLA);
            elseif exist('MADT','var')==1
                MADT = double(MADT);
            end
            U = double(U);
            V = double(V);
            
            type_SSH=Metadata.Type_of_SSH;
            Uinside_cyclo = NaN*U;
            Uinside_anti = NaN*U;
            Vinside_cyclo = NaN*U;
            Vinside_anti = NaN*U;
            
            Label_anti = zeros(length(X),length(Y));
            Label_cyclo = zeros(length(X),length(Y));
            
            if strcmp(type_SSH,'adt') %'sla')
                SSH=ADT;
                %SSH=SLA;
                clear ADT
                %clear SLA
            elseif strcmp(type_SSH,'madt')
                SSH=MADT;
                clear MADT
            elseif strcmp(type_SSH,'surf_el')
                SSH=surf_el;
                clear surf_el
            end
            
            % Compute the Windows size in degree to number of point during the
            % first iteration and create domain size matrices of coordinate
            if i==1
                minDist=ceil(minDist_deg/(X(2)-X(1)));
                disp(['The windows size is ' , num2str(minDist*(X(2)-X(1))), ' degrees, corresponding to ',num2str(minDist), ' pixels '])
                clear minDist_deg;
                [y,x] = meshgrid(Y,X);
            end
            
            % Call Function that extract eddy extraction from closed countours
            % technic. Return position of centers, values of amplitude, contour
            % position and number of eddies
            [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo] = extraction_eddies_new(X,Y,SSH,minDist,min_amp,min_R,2);
            
            
            % Compute and store parameters of each Cyclones if there is at least
            % one detected
            if Ncyclo>0
                for j=1:Ncyclo %Loop on each Cyclone
                    % Edge coordinates
                    Xcon = Cyclonic_Cell{j,5};
                    Ycon = Cyclonic_Cell{j,6};
                    
                    % Change x & y to reduce search area
                    ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);
                    %hold on
                    % Coordinate index of grid points inside the eddy
                    a1=[x(ind1),y(ind1)]'; b1=[Xcon(:),Ycon(:)]';
                    %hold on
                    inp1=inpoly(a1,b1);%pause(2)
                    in = find(inp1==1);%pause(2);
                    %                    max(in)
                    %                     if isempty(in)==1
                    %                         inp1=inpoly(a1,b1);%pause(2)
                    %                         in = find(inp1==1);%pause(2);
                    %                     end
                    INDinside = ind1(in);
                    
                    % Xcenter and Ycenter
                    [~,iiind] = min(SSH(ind1(in)));
                    INDcenter = INDinside(iiind(1));
                    
                    Cyclonic_Cell{j,1}= x(INDcenter);
                    Cyclonic_Cell{j,2}= y(INDcenter);
                    % (Normalized) Vorticity at the center
                    Cyclonic_Cell{j,13}=Vorticity(INDcenter);
                    Cyclonic_Cell{j,14}=Vorticity(INDcenter)./coriolis(y(INDcenter));
                    
                    % Mean EKE
                    Cyclonic_Cell{j,10}=nanmean(EKE(INDinside));
                    
                    % Mean Speed
                    Cyclonic_Cell{j,11}=nanmean(Speed(INDinside));
                    
                    % Mean Vorticity
                    Cyclonic_Cell{j,12}=nanmean(Vorticity(INDinside));
                    
                    % Mean straining, shearing deformation rates and Okubo-Weiss parameter
                    Cyclonic_Cell{j,15}=nanmean(Sn(INDinside));
                    Cyclonic_Cell{j,16}=nanmean(Ss(INDinside));
                    Cyclonic_Cell{j,17}=nanmean(OW(INDinside));
                    
                    % Values of X,Y,U,V inside the eddies
                    Cyclonic_Cell{j,18}=x(INDinside);
                    Cyclonic_Cell{j,19}=y(INDinside);
                    Cyclonic_Cell{j,20}=U(INDinside);
                    Cyclonic_Cell{j,21}=V(INDinside);
                    
                    Uinside_cyclo(INDinside) = U(INDinside);
                    Vinside_cyclo(INDinside) = V(INDinside);
                    
                    Label_cyclo(ind1(in)) = j;
                end
            end
            
            % Compute and store parameters of each Anticyclones if there is at least
            % one detected
            if Nanti>0
                for j=1:Nanti %Loop on each Cyclone
                    % Edge coordinates
                    Xcon = Anticyclonic_Cell{j,5};
                    Ycon = Anticyclonic_Cell{j,6};
                    
                    % Change x & y to reduce search area
                    ind1 = find(x>=min(Xcon)-0.5 & x<=max(Xcon)+0.5 & y>=min(Ycon)-0.5 & y<=max(Ycon)+0.5);
                    
                    % Coordinate index of grid points inside the eddy
                    in = find(inpoly([x(ind1),y(ind1)]',[Xcon(:),Ycon(:)]')==1);
                    INDinside = ind1(in);
                    
                    % Xcenter and Ycenter
                    [~,iiind] = max(SSH(ind1(in)));
                    INDcenter = INDinside(iiind1(1));
                    
                    Anticyclonic_Cell{j,1}= x(INDcenter);
                    Anticyclonic_Cell{j,2}= y(INDcenter);
                    % (Normalized) Vorticity at the center
                    Anticyclonic_Cell{j,13}=Vorticity(INDcenter);
                    Anticyclonic_Cell{j,14}=Vorticity(INDcenter)./coriolis(y(INDcenter));
                    
                    % Mean EKE
                    Anticyclonic_Cell{j,10}=nanmean(EKE(INDinside));
                    
                    % Mean Speed
                    Anticyclonic_Cell{j,11}=nanmean(Speed(INDinside));
                    
                    % Mean Vorticity
                    Anticyclonic_Cell{j,12}=nanmean(Vorticity(INDinside));
                    
                    % Mean straining, shearing deformation rates and Okubo-Weiss parameter
                    Anticyclonic_Cell{j,15}=nanmean(Sn(INDinside));
                    Anticyclonic_Cell{j,16}=nanmean(Ss(INDinside));
                    Anticyclonic_Cell{j,17}=nanmean(OW(INDinside));
                    
                    % Values of X,Y,U,V inside the eddies
                    Anticyclonic_Cell{j,18}=x(INDinside);
                    Anticyclonic_Cell{j,19}=y(INDinside);
                    Anticyclonic_Cell{j,20}=U(INDinside);
                    Anticyclonic_Cell{j,21}=V(INDinside);
                    
                    Uinside_anti(INDinside) = U(INDinside);
                    Vinside_snti(INDinside) = V(INDinside);
                    
                    Label_anti(ind1(in)) = j;
                end
            end
            
            % Save Metada parameters
            %Metadata().Min_Amplitude_cm = min_amp;
            %Metadata().Windows_size_pixel = minDist;
            % Metada().Min_Radius_km = min_R;
            
            % Saves cells and variables in .mat file
            eval([' save -append ' filename ' Metadata Anticyclonic_Cell Cyclonic_Cell Nanti Ncyclo Fields Uinside_cyclo Uinside_anti Vinside_cyclo Vinside_anti Label_cyclo Label_anti'])
            
        end
        close(h)
end
end