
path_root = 'K:\Liu_drought_legacy\sem\variable\';


% 定义文件路径
search_DIF_file = fullfile(path_root, 'Method_1_dif_SOS', '*.tif');
search_sm_l_sum_file = fullfile(path_root, 'gs_drought_sm_l_sum_', '*.tif');

search_EOS_file = fullfile(path_root, 'detrend_EOS_', '*.tif');
search_FLUXCOM_GPP_file = fullfile(path_root, 'detrend_GPP_', '*.tif');
search_spring_tmp_file = fullfile(path_root, 'detrend_spring_tmp_', '*.tif');
search_winter_tmp_file = fullfile(path_root, 'detrend_winter_tmp_', '*.tif');
search_spring_sm_file = fullfile(path_root, 'detrend_spring_sm_', '*.tif');

search_type1_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type1_', '*.tif');
search_type2_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type2_', '*.tif');
search_type3_file = fullfile(path_root, 'gs_NDVI_SM_l_drought_type3_', '*.tif');

% 读取文件
DIF_filenames = dir(search_DIF_file);
sm_l_sum_filenames = dir(search_sm_l_sum_file);


EOS_filenames = dir(search_EOS_file);
FLUXCOM_GPP_filenames = dir(search_FLUXCOM_GPP_file);
spring_tmp_filenames = dir(search_spring_tmp_file);
winter_tmp_filenames = dir(search_winter_tmp_file);
spring_sm_filenames = dir(search_spring_sm_file);

type1_filenames_l = dir(search_type1_file);
type2_filenames_l = dir(search_type2_file);
type3_filenames_l = dir(search_type3_file);

dif_TOTAL = zeros(720, 4320, numel(DIF_filenames));
% 处理每个文件
for file_index = 1:numel(DIF_filenames)
    % 读取TIFF文件
    dif_data = readgeoraster(fullfile(DIF_filenames(i).folder, DIF_filenames(i).name));
    dif_TOTAL(:,:,file_index-1) = dif_data;
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


% 处理sem
pixel = 0;
pixel0 = 0;
pixel_type1_l = 0;
pixel_type2_l = 0;
pixel_type3_l = 0;

for column = 1:720
    for line = 1:4320
        for year = 1:32
            if dif_TOTAL(column, line, year) > -100 && dif_TOTAL(column, line, year) < 100
                sem(1, pixel) = sm_l_sum_TOTAL(column, line, year);
                sem(2, pixel) = FLUXCOM_GPP_TOTAL(column, line, year);
                sem(3, pixel) = EOS_TOTAL(column, line, year);
                sem(4, pixel) = spring_tmp_TOTAL(column, line, year);
                sem(5, pixel) = winter_tmp_TOTAL(column, line, year);
                sem(6, pixel) = spring_sm_TOTAL(column, line, year);
                sem(7, pixel) = dif_TOTAL(column, line, year);

                if max(sem(:, pixel)) < 1e8 && min(sem(:, pixel)) > -1e9 && sem(17, pixel) > 0.00001
                    sem1(:, pixel0) = sem(:, pixel);
                    if drought_type_TOTAL_l(column, line, year) == 1
                        sem_type1_l(:, pixel_type1_l) = sem1(:, pixel0);
                        pixel_type1_l = pixel_type1_l + 1;
                    elseif drought_type_TOTAL_l(column, line, year) == 2
                        sem_type2_l(:, pixel_type2_l) = sem1(:, pixel0);
                        pixel_type2_l = pixel_type2_l + 1;
                    elseif drought_type_TOTAL_l(column, line, year) == 3
                        sem_type3_l(:, pixel_type3_l) = sem1(:, pixel0);
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

header = ['SM loss' , 'GPP anomaly', 'EOS anomaly','ST anomaly', 'WT anomaly','SSM anomaly', 'SOS_change'];
    
fid = fopen('K:\Liu_drought_legacy\sem\sem_variables.csv', 'w');
fprintf(fid, '%s\n', header);
for i = 1:numel(sem1)
    fprintf(fid, '%s\n', sprintf('%g', sem1(i, :)));
end
fclose(fid);




