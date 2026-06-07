sos_obs = readgeoraster('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif');
info = geotiffinfo('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif');
[sos_rows, sos_cols] = size(sos_obs);

% 读取1983-2015年间的SOS观测值
SOS_filenames = dir('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_*.tif');
num_file_SOS = length(SOS_filenames);
SOS_obs_series = zeros(sos_rows, sos_cols, num_file_SOS - 1);

for year = 2:num_file_SOS
    SOS_obs_series(:, :, year - 1) = readgeoraster(fullfile(SOS_filenames(year).folder, SOS_filenames(year).name));
end

% 读取1983-2015年间的预测SOS值
predict_SOS_filenames = dir('K:\Liu_drought_legacy\predicted SOS\sos_predicted_*.tif');
num_predict_SOS_file = length(predict_SOS_filenames);
predict_SOS_series =zeros(sos_rows, sos_cols, num_predict_SOS_file);

for year = 1:num_predict_SOS_file
    predict_SOS_series(:, :, year) = readgeoraster(fullfile(predict_SOS_filenames(year).folder, predict_SOS_filenames(year).name));
end


% 读取1982-2014年间的干旱数据
drought_filenames = dir('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_drought_type_*.tif');
num_file_drought = length(drought_filenames);
drought_series = NaN(sos_rows, sos_cols, num_file_drought);

for year = 1:num_file_drought
    drought_series(:, :, year) = readgeoraster(fullfile(drought_filenames(year).folder, drought_filenames(year).name));
end

% 初始化差异数组
dif_drought_mean = NaN(sos_rows, sos_cols);
dif_no_drought_mean = NaN(sos_rows, sos_cols);

% 计算差异
for i = 1:sos_rows
    for j = 1:sos_cols
        k = 0;
        kk = 0;
        l = 0;
        ll = 0;
        for year = 1:31
            if (SOS_obs_series(i, j, year + 1) > 0 && predict_SOS_series(i, j, year + 1) > 0 && drought_series(i, j, year) ~= 1 && drought_series(i, j, year + 1) >= 1)
                l = l + SOS_obs_series(i, j, year + 1) - predict_SOS_series(i, j, year + 1);
                ll = ll + 1;
            end
            if ll >= 1
                dif_drought_mean(i, j) = l / ll;
            else
                dif_drought_mean(i, j) = NaN;
            end
            
            if (SOS_obs_series(i, j, year) > 0 && predict_SOS_series(i, j, year) > 0 && drought_series(i, j, year) ~= 1 && drought_series(i, j, year + 1) ~= 1)
                k = k + SOS_obs_series(i, j, year + 1) - predict_SOS_series(i, j, year + 1);
                kk = kk + 1;
            end
            if kk >= 1
                dif_no_drought_mean(i, j) = k / kk;
            else
                dif_no_drought_mean(i, j) = NaN;
            end
        end
    end
end

% 写入TIFF文件
query_status = query_tiff('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif', 'infor', 'geotiff', geoinfor);
geotiffwrite('K:\Liu_drought_legacy\SOS_change\method2_Nodrought_SOS_change(obser-predict).tif', dif_no_drought_mean, info.RefMatrix, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('K:\Liu_drought_legacy\SOS_change\method2_Drought_SOS_change(obser-predict).tif', dif_drought_mean, info.RefMatrix, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);





% 写入每个年份的差异
for xx = 0:31
    gs_drought_num_outfile = fullfile('K:\Liu_drought_legacy\SOS_change\Method_2_dif_SOS', num2str(xx+1984, 1), '-', num2str(xx+1983, 1), '.tif');
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
search_dif_SOS_file= fullfile('K:\Liu_drought_legacy\SOS_change\Method_2_dif_SOS', '*.tif');
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
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_2_Drought_SOS_change(delect multiple drought)_type1l.tif',dif_SOS_type1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_2_Drought_SOS_change(delect multiple drought)_type12.tif',dif_SOS_type2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('K:\Liu_drought_legacy\SOS_change\Method_2_Drought_SOS_change(delect multiple drought)_type13.tif',dif_SOS_type3,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag)


