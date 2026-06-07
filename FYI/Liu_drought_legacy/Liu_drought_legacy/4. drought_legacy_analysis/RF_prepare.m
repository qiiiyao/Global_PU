clc
close all
clear all


% 定义文件路径
search_DIF_file = fullfile('K:\Liu_drought_legacy\SOS_change\Method_1_dif_SOS', '*.tif');
path_root = 'K:\Liu_drought_legacy\RF\static variable\';
search_RF_file = fullfile(path_root, '*.tif');

path_root = 'K:\Liu_drought_legacy\RF\dynamic variable\';
search_EOS_file = fullfile(path_root, 'detrend_EOS_', '*.tif');
search_FLUXCOM_GPP_file = fullfile(path_root, 'detrend_GPP_', '*.tif');
search_spring_tmp_file = fullfile(path_root, 'detrend_spring_tmp_', '*.tif');
search_winter_tmp_file = fullfile(path_root, 'detrend_winter_tmp_', '*.tif');
search_spring_sm_file = fullfile(path_root, 'detrend_spring_sm_', '*.tif');

path_root = 'K:\Liu_drought_legacy\drought_identification_NDVI_SM\';
search_start_file = fullfile(path_root, 'difference_start_month_POS_', '*.tif');
search_end_file = fullfile(path_root, 'difference_end_month_l_', '*.tif');
search_sm_l_sum_file = fullfile(path_root, 'gs_drought_sm_l_sum_', '*.tif');
search_type1_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type1_', '*.tif');
search_type2_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type2_', '*.tif');
search_type3_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type3_', '*.tif');

% 读取文件
DIF_filenames = dir(search_DIF_file);
RF_filenames = dir(search_RF_file);
EOS_filenames = dir(search_EOS_file);
FLUXCOM_GPP_filenames = dir(search_FLUXCOM_GPP_file);
spring_tmp_filenames = dir(search_spring_tmp_file);
winter_tmp_filenames = dir(search_winter_tmp_file);
spring_sm_filenames = dir(search_spring_sm_file);
start_filenames = dir(search_start_file);
end_filenames = dir(search_end_file);
sm_l_sum_filenames = dir(search_sm_l_sum_file);
type1_filenames_l = dir(search_type1_file);
type2_filenames_l = dir(search_type2_file);
type3_filenames_l = dir(search_type3_file);

% 处理每个文件

dif_TOTAL = zeros(720, 4320, numel(DIF_filenames));
for file_index = 1:numel(DIF_filenames)
    % 读取TIFF文件
    dif_data = readgeoraster(fullfile(DIF_filenames(i).folder, DIF_filenames(i).name));
    dif_TOTAL(:,:,file_index-1) = dif_data;
end

RF_TOTAL = zeros(720, 4320, numel(RF_filenames));
for file_index = 1:numel(RF_filenames)
    % 读取TIFF文件
    RF_data = readgeoraster(fullfile(RF_filenames(i).folder, RF_filenames(i).name));
    RF_TOTAL(:,:,file_index-1) = RF_data(1:4320, 1:720);
end

EOS_TOTAL = zeros(720, 4320, numel(EOS_filenames)-2);
for file_index = 2:numel(EOS_filenames)-1
    % 读取TIFF文件
    EOS_data = readgeoraster(fullfile(EOS_filenames(i).folder, EOS_filenames(i).name));
    EOS_TOTAL(:,:,file_index-1) = EOS_data(1:4320, 1:720);
end

FLUXCOM_GPP_TOTAL = zeros(720, 4320, numel(FLUXCOM_GPP_filenames)-2);
for file_index = 2:numel(FLUXCOM_GPP_filenames)-1
    % 读取TIFF文件
    FLUXCOM_GPP_data = readgeoraster(fullfile(FLUXCOM_GPP_filenames(i).folder, FLUXCOM_GPP_filenames(i).name));
    FLUXCOM_GPP_TOTAL(:,:,file_index-1) = FLUXCOM_GPP_data(1:4320, 1:720);
end

spring_tmp_TOTAL = zeros(720, 4320, numel(spring_filenames)-2);
for file_index = 3:numel(spring_tmp_filenames)
    % 读取TIFF文件
    spring_tmp_data =  readgeoraster(fullfile(spring_tmp_filenames(i).folder,spring_tmp_filenames(i).name));
    spring_tmp_TOTAL(:,:,file_index-2) = spring_tmp_data(1:4320, 1:720);
end

winter_tmp_TOTAL = zeros(720, 4320, numel(winter_filenames)-2);
for file_index = 3:numel(winter_tmp_filenames)
    % 读取TIFF文件 
    winter_tmp_data = readgeoraster(fullfile(winter_tmp_filenames(i).folder,winter_tmp_filenames(i).name));
    winter_tmp_TOTAL(:,:,file_index-2) = winter_tmp_data(1:4320, 1:720);
end

spring_sm_TOTAL = zeros(720, 4320, numel(spring_sm_filenames)-2);
for file_index = 3:numel(spring_sm_filenames)
    % 读取TIFF文件
    spring_sm_data = readgeoraster(fullfile(spring_sm_filenames(i).folder,spring_sm_filenames(i).name));
    spring_sm_TOTAL(:,:,file_index-2) = spring_sm_data(1:4320, 1:720);
end


start_month_TOTAL = zeros(720, 4320, numel(start_filenames)-1);
for file_index = 2:numel(start_filenames)
    % 读取TIFF文件
    start_data = readgeoraster(fullfile(start_filenames(i).folder,start_filenames(i).name));
    start_month_TOTAL(:,:,file_index-1) = start_data(1:4320, 1:720);
end

end_month_TOTAL = zeros(720, 4320, numel(end_filenames)-1);
for file_index = 2:numel(end_filenames)
    % 读取TIFF文件
    end_data = readgeoraster(fullfile(end_filenames(i).folder,end_filenames(i).name));
    end_month_TOTAL(:,:,file_index-1) = end_data(1:4320, 1:720);
    length_month_TOTAL(:,:,file_index-1) = end_month_TOTAL(:,:,file_index-1) - start_month_TOTAL(:,:,file_index-1);
end

sm_l_sum_TOTAL = zeros(720, 4320, numel(sm_l_sum_filenames)-1);
for file_index = 2:numel(sm_l_sum_filenames)
    % 读取TIFF文件 
    sm_l_sum_data = readgeoraster(fullfile(sm_l_sum_filenames(i).folder,sm_l_sum_filenames(i).name));
    sm_l_sum_TOTAL(:,:,file_index-1) = sm_l_sum_data(1:4320, 1:720);
end

drought_type_TOTAL1_l = zeros(720, 4320, numel(type1_filenames_l)-1);
drought_type_TOTAL2_l = zeros(720, 4320, numel(type1_filenames_l)-1);
drought_type_TOTAL3_l = zeros(720, 4320, numel(type1_filenames_l)-1);
drought_type_TOTAL_l = zeros(720, 4320, numel(type1_filenames_l)-1);

for file_index = 2:numel(type1_filenames_l)
    % 读取TIFF文件
    type1_data = readgeoraster(fullfile(type1_filenames_l(i).folder,type1_filenames(i).name));
    type2_data = readgeoraster(fullfile(type2_filenames_l(i).folder,type2_filenames(i).name));
    type3_data = readgeoraster(fullfile(type3_filenames_l(i).folder,type3_filenames(i).name));
    drought_type_TOTAL1_l(:,:,file_index-1) = type1_data(1:4320, 1:720);
    drought_type_TOTAL2_l(:,:,file_index-1) = type2_data(1:4320, 1:720);
    drought_type_TOTAL3_l(:,:,file_index-1) = type3_data(1:4320, 1:720);
    drought_type_TOTAL_l(:,:,file_index-1) = drought_type_TOTAL1_l(:,:,file_index-1) + 2*drought_type_TOTAL2_l(:,:,file_index-1) + 3*drought_type_TOTAL3_l(:,:,file_index-1);
end


% 处理RF_element
pixel = 0;
pixel0 = 0;
pixel_type1_l = 0;
pixel_type2_l = 0;
pixel_type3_l = 0;


for column = 1:720
    for line = 1:4320
        for year = 1:32
            if dif_TOTAL(column, line, year) > -100 && dif_TOTAL(column, line, year) < 100
                RF_element(1:9, pixel) = RF_TOTAL(column, line,:);
                RF_element(10, pixel) = spring_tmp_TOTAL(column, line, year);
                RF_element(11, pixel) = winter_tmp_TOTAL(column, line, year);
                RF_element(12, pixel) = spring_sm_TOTAL(column, line, year);
                RF_element(13, pixel) = FLUXCOM_GPP_TOTAL(column, line, year);
                RF_element(14, pixel) = EOS_TOTAL(column, line, year);
                RF_element(15, pixel) = start_month_TOTAL(column, line, year);
                RF_element(16, pixel) = end_month_TOTAL(column, line, year);
                RF_element(17, pixel) = length_month_TOTAL(column, line, year);
                RF_element(18, pixel) = sm_l_sum_TOTAL(column, line, year);
                RF_element(19, pixel) = dif_TOTAL(column, line, year);

                if max(RF_element(:, pixel)) < 1e8 && min(RF_element(:, pixel)) > -1e9 && RF_element(17, pixel) > 0.00001
                    RF_element1(:, pixel0) = RF_element(:, pixel);
                    if drought_type_TOTAL_l(column, line, year) == 1
                        RF_element_type1_l(:, pixel_type1_l) = RF_element1(:, pixel0);
                        pixel_type1_l = pixel_type1_l + 1;
                    elseif drought_type_TOTAL_l(column, line, year) == 2
                        RF_element_type2_l(:, pixel_type2_l) = RF_element1(:, pixel0);
                        pixel_type2_l = pixel_type2_l + 1;
                    elseif drought_type_TOTAL_l(column, line, year) == 3
                        RF_element_type3_l(:, pixel_type3_l) = RF_element1(:, pixel0);
                        pixel_type3_l = pixel_type3_l + 1;
                    end
                    pixel0 = pixel0 + 1;
                end
                pixel = pixel + 1;
            end
        end
    end
end

% 写入CSV文件
header = ['MAT', 'MAP', 'AGB', 'Biodiversity', 'Biomes', 'LGS', 'Root depth', 'Soil clay', 'Soil sand', 'isohydricity', 'ST anomaly', 'WT anomaly', 'SSM anomaly', 'GPP anomaly', 'EOS anomaly', 'Timing of drought', 'end of drought', 'Drought duration', 'SM loss', 'SOS_change'];
fid = fopen('K:\Liu_drought_legacy\RF\RF_variables.csv', 'w');
fprintf(fid, '%s\n', header);
for i = 1:numel(RF_element1)
    fprintf(fid, '%s\n', sprintf('%g', RF_element1(i, :)));
end
fclose(fid);

fid = fopen('K:\Liu_drought_legacy\RF\RF_variables_drought_type1.csv', 'w');
fprintf(fid, '%s\n', header);
for i = 1:numel(RF_element_type1_l)
    fprintf(fid, '%s\n', sprintf('%g', RF_element_type1_l(i, :)));
end
fclose(fid);


fid = fopen('K:\Liu_drought_legacy\RF\RF_variables_drought_type2.csv', 'w');
fprintf(fid, '%s\n', header);
for i = 1:numel(RF_element_type2_l)
    fprintf(fid, '%s\n', sprintf('%g', RF_element_type2_l(i, :)));
end
fclose(fid);

fid = fopen('K:\Liu_drought_legacy\RF\RF_variables_drought_type3.csv', 'w');
fprintf(fid, '%s\n', header);
for i = 1:numel(RF_element_type3_l)
    fprintf(fid, '%s\n', sprintf('%g', RF_element_type3_l(i, :)));
end
fclose(fid);



