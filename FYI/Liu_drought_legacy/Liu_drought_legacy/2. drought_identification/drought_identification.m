
sos_search_pattern = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI\month_SOS_20_.tif';
eos_search_pattern = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI\month_EOS_50_.tif';

sos_filenames = dir(sos_search_pattern);
eos_filenames = dir(eos_search_pattern);

num_file_sos = length(sos_filenames);
num_file_eos = length(eos_filenames);

sos_total = zeros(720, 4320, num_file_sos);
eos_total = zeros(720, 4320, num_file_eos);

SM_total = zeros(720, 4320, 24); % 记录两年的数据
ndvi_total = zeros(720, 4320, 24);

gs_drought = zeros(720, 4320, 33); % 记录33年数据
gs_drought1 = zeros(720, 4320, 33);
gs_drought2 = zeros(720, 4320, 33);
gs_drought3 = zeros(720, 4320, 33);

gs_drought1_l = zeros(720, 4320, 33);
gs_drought2_l = zeros(720, 4320, 33);
gs_drought3_l = zeros(720, 4320, 33);

difference_start_month_eos = zeros(720, 4320, 33);
difference_end_month_s_eos = zeros(720, 4320, 33);
difference_end_month_l_eos = zeros(720, 4320, 33);

gs_drought_sm_length_s = zeros(720, 4320, 33); % 记录33年数据
gs_drought_sm_length_l = zeros(720, 4320, 33); % 记录33年数据

gs_drought_sm_s_min = zeros(720, 4320, 33);
gs_drought_sm_s_sum = zeros(720, 4320, 33);
gs_drought_ndvi_s_min = zeros(720, 4320, 33);

gs_drought_sm_l_min = zeros(720, 4320, 33);
gs_drought_sm_l_sum = zeros(720, 4320, 33);
gs_drought_ndvi_l_min = zeros(720, 4320, 33);


SM_name = SM_filenames(file_index).name;% 图像名
drought_data = readgeoraster(strcat(path_root ,SM_name));

% 读取所有SOS和EOS文件
for file_index = 1:num_file_sos
    sos_name = sos_filenames(file_index).name;% 图像名
    sos_total(:,:,file_index) = readgeoraster(strcat(sos_search_pattern ,sos_name));
    eos_name = eos_filenames(file_index).name;% 图像名
    eos_total(:,:,file_index) = readgeoraster(strcat(eos_search_pattern ,eos_name));
end

stdev_sm = readgeoraster('F:\20221206_gleam_sm\SMroot_gimms\dede_stdev\0detrend_deseasonal_SMroot_stdev.tif');
stdev_ndvi = readgeoraster('F:\NDVI_LACC_month\dede_stdev\0detrend_deseasonal_ndvi_stdev.tif');

% 逐年处理
for year = 1982:2014
    SM_root='F:\gleam_sm\SMroot_gimms\dede_stdev\';
    ndvi_root='F:\NDVI_LACC_month\dede_stdev\';

    search_SM_file = sprintf('F:\gleam_sm\SMroot_gimms\dede_stdev\deseasonalized_detrended_SMroot_%d*.tif', year);
    search_ndvi_file = sprintf('F:\NDVI_LACC_month\dede_stdev\deseasonalized_detrended_NDVI_%d*.tif', year);

    SM_filenames = dir(search_SM_file);
    ndvi_filenames = dir(search_ndvi_file);

    % 读取连续2年的SM和NDVI文件
    for file_index = 1:12
        SM_total(:,:,file_index) = readgeoraster(strcat(SM_root ,SM_filenames(file_index).name));
        ndvi_total(:,:,file_index) = readgeoraster(strcat(ndvi_root ,ndvi_filenames(file_index).name));
    end
    search_SM_file1 = sprintf('F:\gleam_sm\SMroot_gimms\dede_stdev\deseasonalized_detrended_SMroot_%d*.tif', year+1);
    search_ndvi_file1 = sprintf('F:\NDVI_LACC_month\dede_stdev\deseasonalized_detrended_NDVI_%d*.tif', year+1);

    SM_filenames1 = dir(search_SM_file1);
    ndvi_filenames1 = dir(search_ndvi_file1);

    for file_index = 13:24
        SM_total(:,:,file_index) = readgeoraster(strcat(SM_root ,SM_filenames1(file_index-12).name));
        ndvi_total(:,:,file_index) = readgeoraster(strcat(ndvi_root ,ndvi_filenames1(file_index-12).name));
    end
    % % % % % % % % % % % % % % % % % % % % % % % %
    sos = sos_total(:,:,year-1981);
    eos = eos_total(:,:,year-1981);
    sos1 = sos_total(:,:,year-1980); % 第2年的SOS

    % 遍历所有像素点
    for i = 1:720
        for j = 1:4320

            if sos(i,j) > 0 && eos(i,j) > 0 && sos1(i,j) > 0 && sos(i,j) <= eos(i,j)
                asos = sos(i,j);
                asos1 = sos1(i,j);
                aeos = eos(i,j);
                lsd = aeos - asos + 1;

                SM_gs = SM_total(i,j,asos:asos1+12-1);
                ndvi_gs = ndvi_total(i,j,asos:asos1+12-1);

                SM_gs_location = find(SM_gs < -0.5 * stdev_sm(i,j));
                SM_gs_num = length(SM_gs_location);
                start_month=nan;
                length=nan;
                start_month_location=nan;



                if SM_gs_num >= 2
                    sum_gs = zeros(SM_gs_num-1, 1);
                    for ii = 1:SM_gs_num-1
                        k = SM_gs_location(ii+1) - SM_gs_location(ii);
                        if length_s > 1 && k > 1
                            break;
                        end

                        if ii == 1
                            if k == 1
                                start_month = SM_gs_location(1);
                                start_month_location = ii;
                                length_s = 2;
                                sum_gs(ii) = 1;
                            end
                        else
                            if sum_gs(ii-1) == 0
                                if k == 1
                                    start_month = SM_gs_location(ii);
                                    start_month_location = ii;
                                    length_s = 2;
                                    sum_gs(ii) = 1;
                                end
                            elseif sum_gs(ii-1) > 0
                                length_s = length_s + 1;
                                sum_gs(ii) = 1;
                            end
                        end
                    end
                end



                % 判断连续干旱的月份是否在生长季内
                if ~isnan(start_month) && ~isnan(length_s) && (start_month + sos(i,j) < eos(i,j))
                    if (length_s + start_month + sos(i,j) - 1) <= eos(i,j)
                        ndvi_drought = ndvi_gs(start_month:start_month+length_s-1);
                    else
                        ndvi_drought = ndvi_gs(start_month:eos(i,j)-sos(i,j));
                    end
                else
                    ndvi_drought = nan;
                end

                event = find(ndvi_drought < -0.5 * stdev_ndvi(i,j));

                if ~isempty(event)
                    difference_start_month_eos(i,j,year-1982) = (sos(i,j) + start_month) - eos(i,j);
                    difference_end_month_s_eos(i,j,year-1982) = (sos(i,j) + start_month + length_s) - eos(i,j);

                    gs_drought_sm_length_s(i,j,year-1982) = length_s;
                    gs_drought_sm_s_min(i,j,year-1982) = min(SM_gs(start_month:start_month+length_s-1));
                    gs_drought_sm_s_sum(i,j,year-1982) = sum(SM_gs(start_month:start_month+length_s-1));
                    gs_drought_ndvi_s_min(i,j,year-1982) = min(ndvi_gs(start_month:start_month+length_s-1));
                    gs_drought_ndvi_s_min(i,j,year-1982) = min(ndvi_gs(start_month:start_month+length_s-1));
                end



                
                % 判断生长季末尾的干旱情况
                end_month = NaN;
                length_l = NaN;
                end_month_location = NaN;

                if sm_gs_num >= 2
                    sum_gs_l = zeros(sm_gs_num-1, 1);
                    for ii = sm_gs_num:-1:2
                        k = sm_gs_location(ii) - sm_gs_location(ii-1);
                        if length_l > 1 && k > 1
                            break;
                        end
                        if ii == sm_gs_num
                            if k == 1
                                end_month = sm_gs_location(sm_gs_num);
                                end_month_location = ii;
                                length_l = 2;
                                sum_gs_l(ii) = 1;
                            end
                        else
                            if sum_gs_l(ii+1) == 0
                                if k == 1
                                    end_month = sm_gs_location(ii);
                                    end_month_location = ii;
                                    length_l = 2;
                                    sum_gs_l(ii) = 1;
                                end
                            elseif sum_gs_l(ii+1) > 0
                                length_l = length_l + 1;
                                sum_gs_l(ii) = 1;
                            end
                        end
                    end
                end

                if ~isnan(end_month) && ~isnan(length_l)
                    if (end_month + sos(i,j) - length_l + 1) >= sos(i,j)
                        if (end_month + sos(i,j)) <= eos(i,j)
                            ndvi_drought_l = ndvi_gs(end_month-length_l+1:end_month);
                        else
                            ndvi_drought_l = ndvi_gs(end_month-length_l+1:eos(i,j)-sos(i,j)+1);
                        end
                    else
                        ndvi_drought_l = NaN;
                    end
                else
                    ndvi_drought_l = NaN;
                end


                counts = find(ndvi_drought_l < -0.5 * stdev_ndvi(i,j));

                % 更新结果矩阵
                if counts >= 1
                    difference_start_month_eos(i, j, year-1982) = (SOS(i, j) + start_month) - EOS(i, j);
                    difference_end_month_s_eos(i, j, year-1982) = (SOS(i, j) + start_month + length) - EOS(i, j);

                    gs_drought_sm_length_s(i, j, year-1982) = length;
                    gs_drought_sm_s_min(i, j, year-1982) = min(SM_gs(start_month:start_month + length - 1));
                    gs_drought_sm_s_sum(i, j, year-1982) = sum(SM_gs(start_month:start_month + length - 1));
                    gs_drought_ndvi_s_min(i, j, year-1982) = min(NDVI_gs(start_month:start_month + length - 1));

                    gs_drought(i, j, year-1982) = 1;

                    if (start_month + length - 1 + SOS(i, j)) <= EOS(i, j)
                        gs_drought3(i, j, year-1982) = 1;
                    elseif (start_month + length - 1 + SOS(i, j)) <= (SOS1(i, j) + 12 - 2)
                        gs_drought2(i, j, year-1982) = 1;
                    else
                        gs_drought1(i, j, year-1982) = 1;
                    end

                    length_l = 0;
                    for xx = start_month:ASOS1 + 12 - 2 - ASOS + 1
                        if SM_gs(xx) > 0
                            break;
                        else
                            length_l = length_l + 1;
                        end
                    end

                    difference_end_month_l_eos(i, j, year-1982) = (SOS(i, j) + start_month + length_l) - EOS(i, j);
                    gs_drought_sm_length_l(i, j, year-1982) = length_l;
                    gs_drought_sm_l_min(i, j, year-1982) = min(SM_gs(start_month:start_month + length_l - 1));
                    gs_drought_sm_l_sum(i, j, year-1982) = sum(SM_gs(start_month:start_month + length_l - 1));
                    gs_drought_ndvi_l_min(i, j, year-1982) = min(NDVI_gs(start_month:start_month + length_l - 1));

                    if (start_month + length_l - 1 + SOS(i, j)) <= EOS(i, j)
                        gs_drought3_l(i, j, year-1982) = 1;
                    elseif (start_month + length_l - 1 + SOS(i, j)) <= (SOS1(i, j) + 12 - 2)
                        gs_drought2_l(i, j, year-1982) = 1;
                    else
                        gs_drought1_l(i, j, year-1982) = 1;
                    end

                else
                    % 设置缺失值为NaN
                    difference_start_month_eos(i, j, year-1982) = NaN;
                    difference_end_month_s_eos(i, j, year-1982) = NaN;
                    gs_drought_sm_length_s(i, j, year-1982) = NaN;
                    gs_drought_sm_s_min(i, j, year-1982) = NaN;
                    gs_drought_sm_s_sum(i, j, year-1982) = NaN;
                    gs_drought_ndvi_s_min(i, j, year-1982) = NaN;

                    difference_end_month_l_eos(i, j, year-1982) = NaN;
                    gs_drought_sm_length_l(i, j, year-1982) = NaN;
                    gs_drought_sm_l_min(i, j, year-1982) = NaN;
                    gs_drought_sm_l_sum(i, j, year-1982) = NaN;
                    gs_drought_ndvi_l_min(i, j, year-1982) = NaN;

                    gs_drought(i, j, year-1982) = NaN;
                    gs_drought3(i, j, year-1982) = NaN;
                    gs_drought2(i, j, year-1982) = NaN;
                    gs_drought1(i, j, year-1982) = NaN;

                    gs_drought3_l(i, j, year-1982) = NaN;
                    gs_drought2_l(i, j, year-1982) = NaN;
                    gs_drought1_l(i, j, year-1982) = NaN;
                end
            end
        end
    end

    info=geotiffinfo('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_1982.tif');
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\difference_start_month_EOS_%d.tif', year), difference_start_month_eos(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\difference_end_month_s_EOS_%d.tif', year), difference_end_month_s_eos(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\difference_end_month_l_EOS_%d.tif', year), difference_end_month_l_eos(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_drought_type_%d.tif', year), gs_drought(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_s_drought_type1_%d.tif', year), gs_drought1(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_s_drought_type2_%d.tif', year), gs_drought2(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_s_drought_type3_%d.tif', year), gs_drought3(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_l_drought_type1_%d.tif', year), gs_drought1_l(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_l_drought_type2_%d.tif', year), gs_drought2_l(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_NDVI_SM_l_drought_type3_%d.tif', year), gs_drought3_l(:,:,year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);


    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_length_s_%d.tif', year), gs_drought_sm_length_s(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_s_min_%d.tif', year),gs_drought_sm_s_min(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_s_sum_%d.tif', year), gs_drought_sm_s_sum(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_ndvi_s_min_%d.tif', year), gs_drought_ndvi_s_min(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_length_l_%d.tif', year),  gs_drought_sm_length_l(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_l_min_%d.tif', year), gs_drought_sm_l_min(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_sm_l_sum_%d.tif', year), gs_drought_sm_l_sum(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    geotiffwrite(sprintf('K:\Liu_drought_legacy\drought_identification_NDVI_SM\gs_drought_ndvi_l_min_%d.tif', year), gs_drought_ndvi_l_min(:, :, year-1982),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);



end
