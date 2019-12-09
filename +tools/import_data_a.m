
fname = '..\data\Soot Data FlareNet 18\20180601_A_CPMA+SMPS.txt';

fid = fopen(fname);

bool_data = 0;
data = [];
data_d = [];

tline = fgetl(fid);
while ischar(tline)
    
    if bool_data % if reading data
        if strcmp(tline(1:4),'Scan') % stop reading daata
            bool_data = 0;
        
        else % format data on the line
            t0 = split(tline); % split tab-delimited string

            t1 = []; % convert strings to numbers
            for ii=1:length(t0)
                t1(ii) = str2double(t0{ii});
            end

            if isnan(t1(1)); t1(1) = []; end % if initial characters are whitespace
            data = [data;t1(2:end)];
            data_d = [data_d;t1(1)];
        end
        
    elseif strcmp(tline(1:8),'Diameter') % find start of the data
        bool_data = 1;
    end
    
    tline = fgetl(fid);
end

fclose(fid);

data_m = [];

