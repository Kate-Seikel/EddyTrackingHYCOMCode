%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        CETTE FONCTION MERGE LES MATRICES JOURNALIERES EN MATRICES       %
%                CYCLONIQUES, ANTICYCLONIQUES ET CONTOURS                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Step3_Merge_eddymat_HYCOM(basepath, fullextract)
addpath(strcat(basepath, 'FUNCTIONS'));
% function Step3_Merge_eddymat(input_dir,output_dir)
%Merge daily matrices of dected eddies in 2 matrices : Cyclones and Antyclones
%Last modification 06/05/14
%Automatic demand to open repertories and specify main parameters
katepath = '/Volumes/Kate-Research/Data/Eddy_Extraction';
load(strcat(basepath, 'EXTRACTION/ConfigFile_fulltime.mat')); %ConfigFile2_fulltime.mat');
output_dir = [strcat(katepath, '/EDDY_PROPERTIES')]; %EDDY_PROPERTIES/2018'];
[status,message,messageid] = mkdir(output_dir);

input_dir= [katepath]; %INDIVIDUAL_FILES/'];
list = dir([ input_dir '/*.mat']);

max_eddies = 0;
ii = 0;
while max_eddies==0;
ii = ii+1;
filename = [input_dir '/' list(ii).name];
load(filename,'Ncyclo','Nanti');
max_eddies = 3*max(Ncyclo,Nanti);
end

if fullextract == 1
    Extraction_Type ='Basic Extraction';
elseif fullextract == 2
    Extraction_Type = 'Full Extraction';
end

% Allow memory to variables

n_day=length(list);

% Common parameters. Add indice _t for total as each matrice has
% variable with same name
date_num_t=NaN(1,n_day);

% Creates uniqu id for Anti and Cyclo

id_anti=NaN(max_eddies,n_day);
id_anti(:,1)=(1:max_eddies);

for i=1:1096 %732:1096 %367:731 %1:366
    id_anti(:,i+1)=(id_anti(max_eddies,i)+1:(id_anti(max_eddies,i)+max_eddies));
end

id_cyclo=id_anti;

switch Extraction_Type
    case 'Basic Extraction'
        %Cyclones parameters
        
        %Eddies
        Xcenter_c=NaN(max_eddies,n_day);
        Ycenter_c=NaN(max_eddies,n_day);
        Xcentroid_c=NaN(max_eddies,n_day);
        Ycentroid_c=NaN(max_eddies,n_day);
        Radius_c=NaN(max_eddies,n_day);
        Area_c=NaN(max_eddies,n_day);
        Amplitude_c=NaN(max_eddies,n_day);
        EKE_c=NaN(max_eddies,n_day);
        Ncyclo_c=NaN(1,n_day);
        CEs=cell(max_eddies,n_day,2);  %Create cells of cyclonic contours CEs{3,2,1}=tables of Xcon(1) of the 3rd eddies the second days
        
        
        %Antiyclones parameters
        
        %Eddies
        Xcenter_a=NaN(max_eddies,n_day);
        Ycenter_a=NaN(max_eddies,n_day);
        Xcentroid_a=NaN(max_eddies,n_day);
        Ycentroid_a=NaN(max_eddies,n_day);
        Radius_a=NaN(max_eddies,n_day);
        Area_a=NaN(max_eddies,n_day);
        Amplitude_a=NaN(max_eddies,n_day);
        EKE_a=NaN(max_eddies,n_day);
        Nanti_a=NaN(1,n_day);
        AEs=cell(max_eddies,n_day,2);   %Idem for Anticyclones
        
        
        Fields_cont = ('Eddies Contours. CEs{3,2,1}=tables of Xcon(1) of the 3rd Cyclones the 2nd days and AEs{4,5,2}=tables of Ycon(2) of the 4th Anticyclones the 5th days');
        

        
        
        %Cyclonic (CEs) and Anticyclonic (AEs) Contours
        
        % Opening and fill new files
        h = waitbar(0,'Extraction eddies');
        
        
        for day_id=1:length(list) %732:1096 %367:731 %1:366
            
%             waitbar(day_id/n_day,h)
            disp(['Extracting variables for:    ' list(day_id).name '     (' int2str(day_id) ')  '])
            filename = [input_dir '/' list(day_id).name]; %change '\' to '/' %!!!!!!!!!
            load(filename)
            
            
            type_SSH=Metadata.Type_of_SSH;
            
            if strcmp(type_SSH,'sla') %'sla')
                renvar SLA SSH
                %renvar SLA SSH
            elseif strcmp(type_SSH,'madt')
                renvar MADT SSH
            elseif strcmp(type_SSH, 'adt')
                renvar ADT SSH
            elseif strcmp(type_SSH, 'surf_el')
                renvar surf_el SSH
            end
                        
            
            if Ncyclo>max_eddies || Nanti>max_eddies
                disp('Number of eddies superior of memory allocated. The program stop prematurely to prevent error');
                break
            end
            
            if day_id==1
                ssh=NaN(size(SSH,1),size(SSH,2),n_day);
                
            end
            
            ssh(:,:,day_id)=SSH;
            
            
            %Map parameters
            date_num_t(1,day_id)=date_num;
            Ncyclo_c(1,day_id)=Ncyclo;
            Nanti_a(1,day_id)=Nanti;
            
            %Cyclones parameters and contour
            
            
            if Ncyclo>0
                for j=1:Ncyclo %pour chaque tourbillon cyclonique
                    
                    %Parameters
                    Xcenter_c(j,day_id) = Cyclonic_Cell{j,1};
                    Ycenter_c(j,day_id) = Cyclonic_Cell{j,2};
                    Xcentroid_c(j,day_id)=Cyclonic_Cell{j,3};
                    Ycentroid_c(j,day_id)=Cyclonic_Cell{j,4};
                    
                    
                    Radius_c(j,day_id)=Cyclonic_Cell{j,7};
                    Area_c(j,day_id)=Cyclonic_Cell{j,8};
                    Amplitude_c(j,day_id)=Cyclonic_Cell{j,9};
                    EKE_c(j,day_id)=Cyclonic_Cell{j,10};
                    
                    %Contours
                    CEs{j,day_id,1} = Cyclonic_Cell{j,5};
                    CEs{j,day_id,2} = Cyclonic_Cell{j,6};
                    
                    
                end
            end
            
            if Nanti>0
                for i=1:Nanti %pour chaque tourbillon anti-cyclonique
                    
                    %Parameters
                    Xcenter_a(i,day_id) = Anticyclonic_Cell{i,1};
                    Ycenter_a(i,day_id) = Anticyclonic_Cell{i,2};
                    Xcentroid_a(i,day_id)=Anticyclonic_Cell{i,3};
                    Ycentroid_a(i,day_id)=Anticyclonic_Cell{i,4};
                    
                    
                    Radius_a(i,day_id)=Anticyclonic_Cell{i,7};
                    Area_a(i,day_id)=Anticyclonic_Cell{i,8};
                    Amplitude_a(i,day_id)=Anticyclonic_Cell{i,9};
                    EKE_a(i,day_id)=Anticyclonic_Cell{i,10};
                    
                    %Contours
                    AEs{i,day_id,1} = Anticyclonic_Cell{i,5};
                    AEs{i,day_id,2} = Anticyclonic_Cell{i,6};
                    
                    
                end
            end
            
        end
        
        B = all(isnan(Area_a),2);
        Xcenter_a(B==1,:) = [];
        Ycenter_a(B==1,:) = [];
        Xcentroid_a(B==1,:) = [];
        Ycentroid_a(B==1,:) = [];
        Radius_a(B==1,:) = [];
        Area_a(B==1,:) = [];
        Amplitude_a(B==1,:) = [];
        EKE_a(B==1,:) = [];
        AEs(B==1,:,:) = [];
        
        
        B = all(isnan(Area_c),2);
        Xcenter_c(B==1,:) = [];
        Ycenter_c(B==1,:) = [];
        Xcentroid_c(B==1,:) = [];
        Ycentroid_c(B==1,:) = [];
        Radius_c(B==1,:) = [];
        Area_c(B==1,:) = [];
        Amplitude_c(B==1,:) = [];
        EKE_c(B==1,:) = [];
        CEs(B==1,:,:) = [];
        
        %Cyclones
        
        %Save matrices
        filename_out_c = [output_dir '/CE_Properties'];
        disp('Save Cyclonic Eddy Properties')
        
        
        %Change name removing type_SSH of eddy indices to fit with needed variables
        Xcenter=single(Xcenter_c);
        Ycenter=single(Ycenter_c);
        Xcentroid=single(Xcentroid_c);
        Ycentroid=single(Ycentroid_c);
        Radius=single(Radius_c);
        Area=single(Area_c);
        Amplitude=single(Amplitude_c);
        EKE=single(EKE_c);
        
        date_num=single(date_num_t);
        Ncyclo=single(Ncyclo_c);
        
        if strcmp(type_SSH,'sla') %'sla')
             SLA=single(ssh);
            %SLA=single(ssh);
        elseif strcmp(type_SSH,'madt')
            MADT=single(ssh);
        elseif strcmp(type_SSH, 'surf_el')
            surf_el=single(ssh);
        end
        
        save(filename_out_c,'Metadata','X','Y','Xcenter','Ycenter','Xcentroid','Ycentroid','Radius','Area','Amplitude','EKE','date_num','Ncyclo','id_cyclo');
        
        %Anticyclones
        
        %Save matrices
        filename_out_a = [output_dir '/AE_Properties'];
        disp('Save Anticyclonic Eddy Properties')
        
        %Change name removing type_SSH of eddy indices to fit with needed variables
        Xcenter=single(Xcenter_a);
        Ycenter=single(Ycenter_a);
        Xcentroid=single(Xcentroid_a);
        Ycentroid=single(Ycentroid_a);
        Radius=single(Radius_a);
        Area=single(Area_a);
        Amplitude=single(Amplitude_a);
        EKE=single(EKE_a);        
        Nanti=single(Nanti_a);
        
        save(filename_out_a,'Metadata','X','Y','Xcenter','Ycenter','Xcentroid','Ycentroid','Radius','Area','Amplitude','EKE','date_num','Nanti','id_anti');
        
        
        %Contours
        %Save matrices
        filename_cont = [output_dir '/Eddy_Contours'];
        disp('Save Eddy Contours')
        save(filename_cont,'-v7.3','X','Y','CEs','AEs','Fields_cont','id_cyclo','Ncyclo','id_anti','Nanti','date_num','Metadata');
        
        close(h);

        
        
        
        
    case 'Full Extraction'
        %Cyclones parameters
        
        %Eddies
        Xcenter_c=NaN(max_eddies,n_day);
        Ycenter_c=NaN(max_eddies,n_day);
        Xcentroid_c=NaN(max_eddies,n_day);
        Ycentroid_c=NaN(max_eddies,n_day);
        Radius_c=NaN(max_eddies,n_day);
        Area_c=NaN(max_eddies,n_day);
        Amplitude_c=NaN(max_eddies,n_day);
        EKE_c=NaN(max_eddies,n_day);
        Speed_c=NaN(max_eddies,n_day);
        Vorticity_c=NaN(max_eddies,n_day);
        Vorticity_center_c=NaN(max_eddies,n_day);
        Vorticity_center_normalized_f_c=NaN(max_eddies,n_day);
        Sn_c=NaN(max_eddies,n_day);
        Ss_c=NaN(max_eddies,n_day);
        OW_c=NaN(max_eddies,n_day);
        Ncyclo_c=NaN(1,n_day);
        CEs=cell(max_eddies,n_day,2);  %Create cells of cyclonic contours CEs{3,2,1}=tables of Xcon(1) of the 3rd eddies the second days
        
        
        %Antiyclones parameters
        
        %Eddies
        Xcenter_a=NaN(max_eddies,n_day);
        Ycenter_a=NaN(max_eddies,n_day);
        Xcentroid_a=NaN(max_eddies,n_day);
        Ycentroid_a=NaN(max_eddies,n_day);
        Radius_a=NaN(max_eddies,n_day);
        Area_a=NaN(max_eddies,n_day);
        Amplitude_a=NaN(max_eddies,n_day);
        EKE_a=NaN(max_eddies,n_day);
        Speed_a=NaN(max_eddies,n_day);
        Vorticity_a=NaN(max_eddies,n_day);
        Vorticity_center_a=NaN(max_eddies,n_day);
        Vorticity_center_normalized_f_a=NaN(max_eddies,n_day);
        Sn_a=NaN(max_eddies,n_day);
        Ss_a=NaN(max_eddies,n_day);
        OW_a=NaN(max_eddies,n_day);
        Nanti_a=NaN(1,n_day);
        AEs=cell(max_eddies,n_day,2);   %Idem for Anticyclones
        
        
        Fields_cont = ('Eddies Contours. CEs{3,2,1}=tables of Xcon(1) of the 3rd Cyclones the 2nd days and AEs{4,5,2}=tables of Ycon(2) of the 4th Anticyclones the 5th days');
        
        
        
        %Cyclonic (CEs) and Anticyclonic (AEs) Contours
        
        % Opening and fill new files
        h = waitbar(0,'Extraction eddies');
        
        for day_id=1:1096 %732:1096 %367:731 %1:366
            
            waitbar(day_id/n_day,h)
            disp(['Extracting variables for:    ' list(day_id).name ''])
            filename = [input_dir '/' list(day_id).name]; %change '\' to '/' %!!!!!!!!!
            load(filename)
            
            
            type_SSH=Metadata.Type_of_SSH;
            
            if strcmp(type_SSH,'sla') %'sla')
                 renvar SLA SSH
                %renvar SLA SSH
            elseif strcmp(type_SSH,'madt')
                renvar MADT SSH
            elseif strcmp(type_SSH,'surf_el')
                renvar surf_el SSH
            end
            
            
            
            if Ncyclo>max_eddies || Nanti>max_eddies
                disp('Number of eddies superior of memory allocated. The program stop prematurely to prevent error');
                break
            end
            
            if day_id==1
                ssh=NaN(size(SSH,1),size(SSH,2),n_day);
                
            end
            
            ssh(:,:,day_id)=SSH;
            
            
            %Map parameters
            date_num_t(1,day_id)=date_num;
            Ncyclo_c(1,day_id)=Ncyclo;
            Nanti_a(1,day_id)=Nanti;
            
            %Cyclones parameters and contour
            
            
            if Ncyclo>0
                for j=1:Ncyclo %pour chaque tourbillon cyclonique
                    
                    %Parameters
                    Xcenter_c(j,day_id) = Cyclonic_Cell{j,1};
                    Ycenter_c(j,day_id) = Cyclonic_Cell{j,2};
                    Xcentroid_c(j,day_id)=Cyclonic_Cell{j,3};
                    Ycentroid_c(j,day_id)=Cyclonic_Cell{j,4};
                    
                    
                    Radius_c(j,day_id)=Cyclonic_Cell{j,7};
                    Area_c(j,day_id)=Cyclonic_Cell{j,8};
                    Amplitude_c(j,day_id)=Cyclonic_Cell{j,9};
                    EKE_c(j,day_id)=Cyclonic_Cell{j,10};
                    Speed_c(j,day_id)=Cyclonic_Cell{j,11};
                    Vorticity_c(j,day_id)=Cyclonic_Cell{j,12};
                    Vorticity_center_c(j,day_id)=Cyclonic_Cell{j,13};
                    Vorticity_center_normalized_f_c(j,day_id)=Cyclonic_Cell{j,14};
                    Sn_c(j,day_id)=Cyclonic_Cell{j,15};
                    Ss_c(j,day_id)=Cyclonic_Cell{j,16};
                    OW_c(j,day_id)=Cyclonic_Cell{j,17};
                    
                    %Contours
                    CEs{j,day_id,1} = Cyclonic_Cell{j,5};
                    CEs{j,day_id,2} = Cyclonic_Cell{j,6};
                    
                    
                end
            end
            
            if Nanti>0
                for i=1:Nanti %pour chaque tourbillon anti-cyclonique
                    
                    %Parameters
                    Xcenter_a(i,day_id) = Anticyclonic_Cell{i,1};
                    Ycenter_a(i,day_id) = Anticyclonic_Cell{i,2};
                    Xcentroid_a(i,day_id)=Anticyclonic_Cell{i,3};
                    Ycentroid_a(i,day_id)=Anticyclonic_Cell{i,4};
                    
                    
                    Radius_a(i,day_id)=Anticyclonic_Cell{i,7};
                    Area_a(i,day_id)=Anticyclonic_Cell{i,8};
                    Amplitude_a(i,day_id)=Anticyclonic_Cell{i,9};
                    EKE_a(i,day_id)=Anticyclonic_Cell{i,10};
                    Speed_a(i,day_id)=Anticyclonic_Cell{i,11};
                    Vorticity_a(i,day_id)=Anticyclonic_Cell{i,12};
                    Vorticity_center_a(i,day_id)=Anticyclonic_Cell{i,13};
                    Vorticity_center_normalized_f_a(i,day_id)=Anticyclonic_Cell{i,14};
                    Sn_a(i,day_id)=Anticyclonic_Cell{i,15};
                    Ss_a(i,day_id)=Anticyclonic_Cell{i,16};
                    OW_a(i,day_id)=Anticyclonic_Cell{i,17};
                    
                    %Contours
                    AEs{i,day_id,1} = Anticyclonic_Cell{i,5};
                    AEs{i,day_id,2} = Anticyclonic_Cell{i,6};
                    
                    
                end
            end
            
        end
        
        B = all(isnan(Area_a),2);
        Xcenter_a(B==1,:) = [];
        Ycenter_a(B==1,:) = [];
        Xcentroid_a(B==1,:) = [];
        Ycentroid_a(B==1,:) = [];
        Radius_a(B==1,:) = [];
        Area_a(B==1,:) = [];
        Amplitude_a(B==1,:) = [];
        EKE_a(B==1,:) = [];
        Speed_a(B==1,:) = [];
        Vorticity_a(B==1,:) = [];
        Vorticity_center_a(B==1,:) = [];
        Vorticity_center_normalized_f_a(B==1,:) = [];
        Sn_a(B==1,:) = [];
        Ss_a(B==1,:) = [];
        OW_a(B==1,:) = [];
        AEs(B==1,:,:) = [];
        
        
        B = all(isnan(Area_c),2);
        Xcenter_c(B==1,:) = [];
        Ycenter_c(B==1,:) = [];
        Xcentroid_c(B==1,:) = [];
        Ycentroid_c(B==1,:) = [];
        Radius_c(B==1,:) = [];
        Area_c(B==1,:) = [];
        Amplitude_c(B==1,:) = [];
        EKE_c(B==1,:) = [];
        Speed_c(B==1,:) = [];
        Vorticity_c(B==1,:) = [];
        Vorticity_center_c(B==1,:) = [];
        Vorticity_center_normalized_f_c(B==1,:) = [];
        Sn_c(B==1,:) = [];
        Ss_c(B==1,:) = [];
        OW_c(B==1,:) = [];
        CEs(B==1,:,:) = [];
        
        %Cyclones
        
        %Save matrices
        filename_out_c = [output_dir '/CE_Properties'];
        disp('Save Cyclonic Eddy Properties')
        
        
        %Change name removing type_SSH of eddy indices to fit with needed variables
        Xcenter=Xcenter_c;
        Ycenter=Ycenter_c;
        Xcentroid=Xcentroid_c;
        Ycentroid=Ycentroid_c;
        Radius=Radius_c;
        Area=Area_c;
        Amplitude=Amplitude_c;
        EKE=EKE_c;
        Speed=Speed_c;
        Vorticity=Vorticity_c;
        Vorticity_center=Vorticity_center_c;
        Vorticity_center_normalized_f=Vorticity_center_normalized_f_c;
        Sn=Sn_c;
        Ss=Ss_c;
        OW=OW_c;
        
        date_num=date_num_t;
        Ncyclo=Ncyclo_c;
        
        if strcmp(type_SSH,'sla') %'sla')
             SLA=ssh;
            %SLA=ssh;
        elseif strcmp(type_SSH,'madt')
            MADT=ssh;
        elseif strcmp(type_SSH,'surf_el')
            surf_el=ssh;
        end
        
        save(filename_out_c,'-v7.3','Metadata','X','Y','Xcenter','Ycenter','Xcentroid','Ycentroid','Radius','Area','Amplitude','EKE','Speed','Vorticity','Vorticity_center','Vorticity_center_normalized_f','Ss','Sn','OW','date_num','Ncyclo','id_cyclo');
        
        %Anticyclones
        
        %Save matrices
        filename_out_a = [output_dir '/AE_Properties'];
        disp('Save Antiyclonic Eddy Properties')
        
        %Change name removing type_SSH of eddy indices to fit with needed variables
        Xcenter=Xcenter_a;
        Ycenter=Ycenter_a;
        Xcentroid=Xcentroid_a;
        Ycentroid=Ycentroid_a;
        Radius=Radius_a;
        Area=Area_a;
        Amplitude=Amplitude_a;
        EKE=EKE_a;
        Speed=Speed_a;
        Vorticity=Vorticity_a;
        Vorticity_center=Vorticity_center_a;
        Vorticity_center_normalized_f=Vorticity_center_normalized_f_a;
        Sn=Sn_a;
        Ss=Ss_a;
        OW=OW_a;
        
        Nanti=Nanti_a;
        
        save(filename_out_a,'-v7.3','Metadata','X','Y','Xcenter','Ycenter','Xcentroid','Ycentroid','Radius','Area','Amplitude','EKE','Speed','Vorticity','Vorticity_center','Vorticity_center_normalized_f','Ss','Sn','OW','date_num','Nanti','id_anti');
        
        
        %Contours
        %Save matrices
        filename_cont = [output_dir '/Eddy_Contours'];
        disp('Save Eddy Contours')
        save(filename_cont,'-v7.3','X','Y','CEs','AEs','Fields_cont','id_cyclo','Ncyclo','id_anti','Nanti','date_num','Metadata');
        
        close(h);
        
end
end


