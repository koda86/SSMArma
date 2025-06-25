clear all;
close all;

% chdir('C:\Users\algona_workstation\Documents\GitHub\CVPR24_PointSetReg\matlab_code');

addpath('./src');
addpath('./utils/');

%%
%##########################################################################
%Einlesen

%Ort wo Daten sind
% f = dir('D:\Tu Chemnitz\Arbeiten\Masterarbeit\Matlab\Daten DS even2');
f = dir('C:\Users\algona_workstation\Documents\SSM_Arma\Daten_ArmaSuisse\Export_Nicklas');
f = dir('C:\Users\algona_workstation\Documents\SSM_Arma\Daten_ArmaSuisse\unprocessed_test');
names={f(~[f.isdir]).name};
names=string(names);

Pointclouds = {};
[m n] = size(names);

%Ort wo Daten sind
% cd 'D:\Tu Chemnitz\Arbeiten\Masterarbeit\Matlab\Daten DS even2';
cd 'C:\Users\algona_workstation\Documents\SSM_Arma\Daten_ArmaSuisse\Export_Nicklas';
cd 'C:\Users\algona_workstation\Documents\SSM_Arma\Daten_ArmaSuisse\unprocessed_test';

for i=1:n
    Pointclouds(i,1) = {pcread(names(1,i))};
    Pointclouds{i,1}.Color=[0.706 0.706 1];
end

%%
%%#########################################################################
%Referenz auswaehlen

indx_ref_pointcloud = 2;
reference_pointcloud = Pointclouds{indx_ref_pointcloud,1};
Pointclouds_plot = Pointclouds;
Pointclouds_plot(indx_ref_pointcloud,:)= [];


%%
%##########################################################################
%Downsampling der besten Referenz

%reference_pointcloud = pcdownsample(reference_pointcloud,"gridAverage",13);


%%
%##########################################################################
%Kontrollplot vor Registration

figure;

s1 = "Referenz und Punktwolke ";

pos=1;
for i=1:n
    if i~=indx_ref_pointcloud
    subplot(2,4,pos)
    s2= num2str(i);

    pcshowpair(reference_pointcloud,Pointclouds{i,1}); hold on
    pos=pos+1;
    s3=strcat(s1,s2);

    title(s3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    end
end

%%
%##########################################################################
% Non Rigid Registration der Referenz auf Punktewolken

% tic;
% registered_pointclouds = {};
% pos=1;
% 
% for i = 1:n
% 
%         if i~=indx_ref_pointcloud
% 
% 
%         [tform,pointcloud_reg,rmse] = pcregistercpd(reference_pointcloud,Pointclouds{i,1},"Transform","Nonrigid","InteractionSigma",1.5,"SmoothingWeight",0.1);
% 
%         pointcloud_reg.Color = [1 0 0];
%         registered_pointclouds(pos,1)={pointcloud_reg};
% 
%         pos=pos+1;
%         end
% 
% end
% toc;
%%
%##########################################################################
% Alternative Registration der Referenz auf Punktewolken
% tic;
% registered_pointclouds = {};
% pos=1;
% 
% for i=1:n
%     if i~=indx_ref_pointcloud
% 
%     tform = registerPoints3d_affine(reference_pointcloud.Location,Pointclouds{i,1}.Location,20);
%     temp = transformVector3d(reference_pointcloud.Location,tform);
%     temp_pt = pointCloud(temp, Color=[1 0 0]);
%     registered_pointclouds(pos,1)={temp_pt};
%     pos=pos+1;
%     end
% end
% toc;
%%
%##########################################################################
% Alternative Registration der Referenz auf Punktewolken

registered_pointclouds = {};
pos=1;
tic;

for i = 1:n-1
    src=reference_pointcloud;
    tgt=Pointclouds_plot{i,1};
    
    
    % Normalize the point sets
    src_pt=src.Location;
    tgt_pt=tgt.Location;
    
    [src_pt_normal,src_pre_normal]=data_normalize_input(src_pt);
    [tgt_pt_normal,tgt_pre_normal]=data_normalize_input(tgt_pt);
    
    % Downsample point sets make their size less than 5,000
    src_pt_normal=pointCloud(src_pt_normal);
    tgt_pt_normal=pointCloud(tgt_pt_normal);
    
    % gridStep=0.03; 
    % 
    % src_pt_normal=pcdownsample(src_pt_normal,'gridAverage',gridStep); 
    % tgt_pt_normal=pcdownsample(tgt_pt_normal,'gridAverage',gridStep);
    
    src_pt_normal=double(src_pt_normal.Location);
    tgt_pt_normal=double(tgt_pt_normal.Location);
    
    % Show the normalized source and target point clouds
    % figure;
    % subplot(1,2,1)
    % scatter3(src_pt_normal(:,1),src_pt_normal(:,2),src_pt_normal(:,3),'filled');
    % title("source")
    % subplot(1,2,2)
    % scatter3(tgt_pt_normal(:,1),tgt_pt_normal(:,2),tgt_pt_normal(:,3),'filled');
    % title("target")
    % hold off;
    
    src_pt_normal_gpu=gpuArray(src_pt_normal);
    tgt_pt_normal_gpu=gpuArray(tgt_pt_normal);
    [alpha,T_deformed]=fuzzy_cluster_reg(src_pt_normal_gpu,tgt_pt_normal_gpu);
    
    % Denormalize the deformed point cloud
    T_deformed_denormal=denormalize(tgt_pre_normal,T_deformed);
    tgt_pt_denormal=denormalize(tgt_pre_normal,tgt_pt_normal);
    
    % Show the original source and target point clouds
    % figure;
    % subplot(1,2,1)
    % scatter3(src_pt(:,1),src_pt(:,2),src_pt(:,3),'filled');
    % title("source")
    % subplot(1,2,2)
    % scatter3(tgt_pt(:,1),tgt_pt(:,2),tgt_pt(:,3),'filled');
    % title("target")
    % hold off;
    
    % Show the target and deformed point clouds
    % figure;
    % hold on;
    % scatter3(tgt_pt_denormal(:,1),tgt_pt_denormal(:,2),tgt_pt_denormal(:,3),'filled');
    % scatter3(T_deformed_denormal(:,1),T_deformed_denormal(:,2),T_deformed_denormal(:,3),'filled');
    % title("Registration")
    % hold off;
    
    temp = gather(T_deformed_denormal);
    
    tgt_pt_denormal_pt=pointCloud(tgt_pt_denormal,Color=[0.706 0.706 1]);
    T_deformed_denormal_pt = pointCloud(temp,Color = [0.706 0.706 1]);
    
    registered_pointclouds(i,1)={T_deformed_denormal_pt};
    registered_pointclouds(i,2)={tgt_pt_denormal_pt};
end

toc;
%%
%##########################################################################
% Kontrollplot für Registration 


figure;

s1 = "Referenz und Punktwolke ";
s3 =  " Nach Registration";
pos=1;
for i=1:n-1

    
    subplot(2,4,pos)
    s2= num2str(i);

    pcshowpair(registered_pointclouds{i,1},Pointclouds_plot{i,1}); hold on
    %pcshow(registered_pointclouds{i,1},"MarkerSize",50); hold on
    %pcshow(Pointclouds_plot{i,1});
    pos=pos+1;
    s4=strcat(s1,s2,s3);

    title(s4);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    

end


s5 = 'Referenz Ursprl. mit Registrierter Referenz von Punktewolke ';
pos=1;
figure;

for i=1:n-1

    
    subplot(2,4,pos)
    s2= num2str(i);

    pcshowpair(reference_pointcloud,registered_pointclouds{i,1}); hold on
    %pcshow(registered_pointclouds{i,1},"MarkerSize",50); hold on
    %pcshow(Pointclouds_plot{i,1});
    pos=pos+1;
    s4=strcat(s5,s2);

    title(s4);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    

end

s6 = "Abgleich Ursprl. Punktewolke ";
s8 = " Denomalen Punktwolke" ;
pos=1;

figure;
for i = 1:n-1
subplot(2,4,pos)
s7= num2str(i);
s9= strcat(s6,s7,s8);
pcshowpair(registered_pointclouds{i,2},Pointclouds_plot{i,1}); hold on
pos=pos+1;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(s9);
end
%%
%##########################################################################
%Zufallsindex als Kontrolle

anz= 15;
[mr nr] = size(reference_pointcloud.Location);
indx = randi([1 mr],1,anz);
%indx2 = [38	224	1902	512	1907	1015	877	494	1736	779	182	663	1553	1656	112];



pos=2;
figure;

subplot(3,3,1);

pcshow(reference_pointcloud); hold on
pcshow(reference_pointcloud.Location(indx,:),'r','Markersize',100);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Indexe an der urlspr. Referenz');

s10 = "Indexe an der registrierten Referenz auf Punktewolke ";
for i = 1:n-1
    
subplot(3,3,pos);


pcshow(registered_pointclouds{i,1}); hold on
pcshow(registered_pointclouds{i,1}.Location(indx,:),'r','Markersize',100);
s11=num2str(i);
s12= strcat(s10,s11);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(s12);
pos=pos+1;


end

l=length(indx);
color_plot = ["r","g","b","y","m","c","k","w","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];


figure;

pcshow(reference_pointcloud); hold on
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Superimposed Registered Pointclouds')

for i=1:l
pcshow(reference_pointcloud.Location(indx(1,i),:),color_plot(i),'Markersize',100)

end

for i=1:n-1
    pcshow(registered_pointclouds{i,1});

    for j=1:l
    pcshow(registered_pointclouds{i,1}.Location(indx(1,j),:),color_plot(j),'Markersize',100)


    end


end

%%
%##########################################################################
%Export

% WICHTIG: Export als obj-Format, das von importWavefront() gelesen werden
% kann (d. h. inkl. faces)!!!

EM = {};
EM(1,1) = {reference_pointcloud};

for i=1:n-1
    EM(i+1,1)= registered_pointclouds(i,1); 
end

se1 = "Test";

for i = 1:n
    se2 = num2str(i);
    se3 = strcat(se1,se2);
    
    pcwrite(EM{i,1},se3)
end

