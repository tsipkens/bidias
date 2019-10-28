
clear;
close all;
fclose('all');
clc;

fid = fopen('..\data\Soot Data FlareNet 18\20180601_E_CPMA+SMPS.txt');

foo = '';
while ~strcmp(foo,'Diameter Midpoint')
    foo = fgetl(fid);
    t0 = textscan(foo,'%s','delimiter','\t');
    
    if strcmp(t0{1}{1},'DMA Inner Radius(cm)'); prop.R1=str2double(t0{1}{2}); end
    if strcmp(t0{1}{1},'DMA Outer Radius(cm)'); prop.R2=str2double(t0{1}{2}); end
    if strcmp(t0{1}{1},'DMA Characteristic Length(cm)'); prop.L=str2double(t0{1}{2}); end
    if strcmp(t0{1}{1},'Reference Gas Temperature (K)'); prop.T=str2double(t0{1}{2}); end
    if strcmp(t0{1}{1},'Reference Gas Pressure (kPa)'); prop.p=str2double(t0{1}{2})/101.325; end
    if strcmp(t0{1}{1},'Sample #'); num=str2double(t0{1}{end}); end
end

ii=0;
foo = fgetl(fid);
t2 = repmat('%f',1,num+1);
t0 = textscan(foo,t2,'delimiter','\t');
data = [];
while ~contains(foo,'Scan Up Time(s)')
    ii=ii+1;
    data(ii,:) = [t0{:}];
    foo = fgetl(fid);
    t0 = textscan(foo,t2,'delimiter','\t');
end

d = data(:,1);
data = data(:,2:end)';

foo = fgetl(fid);
while foo~=-1
    t0 = textscan(foo,'%s','delimiter','\t');
    
    if strcmp(t0{1}{1},'Sheath Flow(lpm)'); prop.Q_c=str2double(t0{1}{2})/60/1000; end
    if strcmp(t0{1}{1},'Aerosol Flow(lpm)'); prop.Q_a=str2double(t0{1}{2})/60/1000; end
    
    foo = fgetl(fid);
end

fclose(fid);

