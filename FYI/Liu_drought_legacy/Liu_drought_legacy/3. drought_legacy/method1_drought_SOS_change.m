path_root = 'K:\Liu_drought_legacy\drought_identification_NDVI_SM';
path_root1 = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI';
% 定义变量

num_file_drought = 33;
num_file_SOS = 34; 

% 初始化所需数组
drought_TOTAL = zeros(720, 4320, num_file_drought);
SOS_TOTAL = zeros(720, 4320, num_file_SOS);
SOS_TOTAL0 = zeros(720, 4320, num_file_SOS);
number = zeros(720, 4320);

% 定义文件路径
search_SM_file = fullfile(path_root, 'gs_NDVI_SM_drought_type_', '*.tif');
SM_filenames = dir(search_SM_file);

% 读取文件
for file_index = 1:numel(SM_filenames)
    [drought_data,R] = readgeoraster(strcat(path_root ,SM_filenames(file_index).name));  

    drought_TOTAL(:,:,file_index-1) = drought_data;
end

search_SOS_file = fullfile(path_root1, 'SOS_20_', '*.tif');
SOS_filenames = dir(search_SOS_file);

% 读取文件
for file_index = 1:numel(SOS_filenames)
    [SOS_data,R] = readgeoraster(strcat(path_root1 ,SOS_filenames(file_index).name));     
    SOS_TOTAL(:,:,file_index-1) = SOS_data;
end

% 计算SOS变化
diff = zeros(720, 4320, num_file_drought);
dif0 = zeros(720, 4320, num_file_drought);

for column = 1:720
    for line = 1:4320
        if SOS_TOTAL(column, line, 1) ~= SOS_TOTAL(1, 1, 1)
            for year = 1:num_file_drought-1
                if drought_TOTAL(column, line, year) ~= 1 && drought_TOTAL(column, line, year+1) == 1
                    diff(column, line, year) = SOS_TOTAL(column, line, year+1) - SOS_TOTAL(column, line, year);
                elseif drought_TOTAL(column, line, year) ~= 1 && drought_TOTAL(column, line, year+1) ~= 1
                    dif0(column, line, year) = SOS_TOTAL(column, line, year+1) - SOS_TOTAL(column, line, year);
                else
                    diff(column, line, year) = NaN;
                    dif0(column, line, year) = NaN;
                end
            end
        else
            diff(column, line, :) = NaN;
            dif0(column, line, :) = NaN;
        end
    end
end

% 计算平均值
dif_SOS_value = mean(diff, 3, 'omitnan');
dif0_SOS_value = mean(dif0, 3, 'omitnan');

% 写入TIFF文件
query_status = query_tiff('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif', 'infor', 'geotiff', geoinfor);
write_tiff(strcat('K:\Liu_drought_legacy\SOS_change\Method_1_Drought_SOS_change(delect multiple drought)_type all.tif'), dif_SOS_value, 'float', 'geotiff', geoinfor);
write_tiff(strcat('K:\Liu_drought_legacy\SOS_change\Method_1_Nodrought_SOS_change.tif'), dif0_SOS_value, 'float', 'geotiff', geoinfor);



% 写入每个年份的差异
for xx = 0:31
    gs_drought_num_outfile = fullfile('K:\Liu_drought_legacy\SOS_change\Method_1_dif_SOS', num2str(xx+1984, 1), '-', num2str(xx+1983, 1), '.tif');
    write_tiff(gs_drought_num_outfile, diff(:, :, xx), 'float', 'geotiff', geoinfor);
end


path_root = 'K:\Liu_drought_legacy\drought_identification_NDVI_SM\'; 


% 初始化所需数组
type1_TOTAL = zeros(720, 4320, 33);
type2_TOTAL = zeros(720, 4320, 33);
type3_TOTAL = zeros(720, 4320, 33);
type_TOTAL = zeros(720, 4320, 33);

% 定义文件路径
search_type1_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type1_', '*.tif');
search_type2_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type2_', '*.tif');
search_type3_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type3_', '*.tif');

% 读取文件
type1_filenames = dir(search_type1_file);
type2_filenames = dir(search_type2_file);
type3_filenames = dir(search_type3_file);

% 处理每个文件
for file_index = 1:numel(type1_filenames)
    % 读取TIFF文件
    type1_data = readgeoraster(strcat(path_root ,type1_filenames(file_index).name));
    type1_TOTAL(:,:,file_index) = type1_data;
end

for file_index = 1:numel(type2_filenames)
    % 读取TIFF文件
    
    type2_data = readgeoraster(strcat(path_root ,type2_filenames(file_index).name));
    type2_TOTAL(:,:,file_index) = type2_data;
end

for file_index = 1:numel(type3_filenames)
    % 读取TIFF文件
    type3_data = readgeoraster(strcat(path_root ,type3_filenames(file_index).name));
    type3_TOTAL(:,:,file_index) = type3_data;
end

% 计算总类型
type_TOTAL = type1_TOTAL + 2*type2_TOTAL + 3*type3_TOTAL;

% 写入TIFF文件
for file_index = 1:numel(type1_filenames)
    % 创建文件路径
    gs_drought_num_outfile = fullfile(path_root, 'gs_NDVI_SM_l_drought_types(123)_', num2str(file_index+1981, 1), '.tif');
    % 写入TIFF文件
    imwrite(type_TOTAL(:,:,file_index), gs_drought_num_outfile, 'TIFF', 'GeoKeyDirectory', geoinfor);
end

% 处理SOS变化

path_root = 'K:\Liu_drought_legacy\SOS_change\';
search_dif_SOS_file= fullfile('K:\Liu_drought_legacy\SOS_change\Method_1_dif_SOS', '*.tif');
dif_SOS_filenames = dir(search_dif_SOS_file);
dif_SOS_TOTAL = zeros(720, 4320, numel(dif_SOS_filenames));

% 读取文件
for file_index = 1:numel(dif_SOS_filenames)
    % 读取TIFF文件
    dif_SOS_data = readgeoraster(strcat(path_root ,dif_SOS_filenames(file_index).name));
    dif_SOS_TOTAL(:,:,file_index) = dif_SOS_data;
end

% 计算SOS变化
dif_SOS_type1 = zeros(720, 4320);
dif_SOS_type2 = zeros(720, 4320);
dif_SOS_type3 = zeros(720, 4320);

for i = 1:720
    for j = 1:4320
        aa = 0; bb = 0; cc = 0;
        a = 0; b = 0; c = 0;
        for k = 2:33
            if type_TOTAL(i,j,k) == 1 && dif_SOS_TOTAL(i,j,k-1) > -1000
                a = a + dif_SOS_TOTAL(i,j,k-1);
                aa = aa + 1;
            elseif type_TOTAL(i,j,k) == 2 && dif_SOS_TOTAL(i,j,k-1) > -1000
                b = b + dif_SOS_TOTAL(i,j,k-1);
                bb = bb + 1;
            elseif type_TOTAL(i,j,k) == 3 && dif_SOS_TOTAL(i,j,k-1) > -1000
                c = c + dif_SOS_TOTAL(i,j,k-1);
                cc = cc + 1;
            end
        end

        if aa > 0
            dif_SOS_type1(i,j) = a / aa;
        else
            dif_SOS_type1(i,j) = NaN;
        end

        if bb > 0
            dif_SOS_type2(i,j) = b / bb;
        else
            dif_SOS_type2(i,j) = NaN;
        end

        if cc > 0
            dif_SOS_type3(i,j) = c / cc;
        else
            dif_SOS_type3(i,j) = NaN;
        end
        
    end
end

% 写入TIFF文件

info=geotiffinfo('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif');
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_1_Drought_SOS_change(delect multiple drought)_type1l.tif',dif_SOS_type1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_1_Drought_SOS_change(delect multiple drought)_type12.tif',dif_SOS_type2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_1_Drought_SOS_change(delect multiple drought)_type13.tif',dif_SOS_type3,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)


