
NDVI = importdata('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\GIMMS_NDVI.mat');
NDVI_NH=NDVI(1:720,:,:);
clear NDVI

folderPath = 'K:\soil_temperature\';
folderPath1 = 'K:\soil_temperature\soil_temperature_avg_';
filelist=dir([folderPath1 ,'*.tif']);
ST=zeros(720,4320,816);
for i=1:816
    filename = [folderPath, filelist(i).name];
    ST(:,:,i)=imread(filename);
end


NDVI_nosnow=zeros(720,4320,816);
NDVI_NH_nosnow_sgf=zeros(720,4320,816);
NDVI_background=zeros(720,4320);

for xx=1:720
    for yy=1:4320
        if NDVI_NH(xx,yy,1) ~= 0 && ST(xx,yy,1) > 0
           
            st_1=reshape(ST(xx,yy,:)\100-273.15,816,1);
            ndvi_1=reshape(double(NDVI_NH(xx,yy,:))\10000,816,1);
            % figure;
            % plot(st_1,'.-');
            % figure;
            % plot(ndvi_1,'.-');
            %温度大于0保持原始值
            ndvi_3=ndvi_1;
            ndvi_2=ndvi_1;
            %location=st_1>1;
             location=st_1>0;
            ndvi_2(location)= NaN;
            figure;
            plot(ndvi_2,'ro');
            hold on
            plot(ndvi_1,'b.-');
% 22,46,70,94
            figure;
            plot(ndvi_2((2009-1982)*24+1:(2012-1982+1)*24),'ro');
            hold on
            plot(ndvi_1((2009-1982)*24+1:(2012-1982+1)*24),'b.-');
           
           % figure;
            %plot(st_1((2009-1982)*24+1:(2012-1982+1)*24),'.-');

            %location=find(st_1>0.05);
            location=find(st_1>0 & ndvi_1>0.1);
           % location=find(st_1>1 & ndvi_1>0.1);
            %figure;
            %plot(ndvi_1(location),'.-');

            if length(location) > 82 %
            % 将数组按升序排列a
            sorted_data = sort(ndvi_1(location));
            figure;
            plot(sorted_data,'.-');

            % 计算数组的长度
            array_length = length(sorted_data);

            % 找到5%和10%位置的索引 [取一段数据的平均值]
            index_5_percent = ceil(0.05 * array_length);
            index_10_percent = ceil(0.05 * array_length);

            % 截取在5%到10%之间的子数组
            sub_array = sorted_data(index_5_percent:index_10_percent);

            % 计算子数组的平均值
            mean_value = mean(sub_array);
            NDVI_background(xx,yy)=mean_value;

            % 显示结果
            %disp(['5%到10%数值的平均值为：', num2str(mean_value)]);

            ndvi_1(ndvi_1<mean_value)=mean_value;
            %hold on
            %plot(ndvi_1(1:49),'.-');
            NDVI_nosnow(xx,yy,:)=ndvi_1*10000; %%%%%%%%%%%%%%%%没有雪的ndvi           
            order = 4;
            framelen =7;%*****************
            sgf_NDVI_series = sgolayfilt(double(ndvi_1),order,framelen);
            %hold on
            % plot(sgf_NDVI_series,'.-');
            %plot(sgf_NDVI_series(24*28+1:24*32),'.-');
            NDVI_NH_nosnow_sgf(xx,yy, :) = sgf_NDVI_series*10000;%%%%%%%%%%%%%%%%没有雪的sg  
            end
        end
    end
end


filename = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI\NDVI_NH_soiltmp_sgf.mat';
save(filename, 'NDVI_NH_nosnow_sgf');


filename = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI\NDVI_NH_soiltmp.mat';
save(filename, 'NDVI_nosnow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot(reshape(NDVI_NH_nosnow_sgf(300,3500,:),816,1),'.-')

ndvi_mean=mean(NDVI_NH_nosnow_sgf,3);

cc=zeros(720,4320,24);
cc_min=zeros(720,4320,34);
cc_max=zeros(720,4320,34);
for year=1:34 
    cc=NDVI_NH_nosnow_sgf(:,:,24*(year-1)+1:24*year);
    cc_max(:,:,year)=max(cc,[],3);
    cc_min(:,:,year)=min(cc,[],3);
end
ndvi_max=mean(cc_max,3);
ndvi_min=mean(cc_min,3);


var3=NDVI_NH_nosnow_sgf;
for i=1:720
    for j=1:4320
        if (ndvi_mean(i,j)<1000)|| (ndvi_max(i,j)<2000)||(ndvi_max(i,j)<1.2*ndvi_min(i,j))
            var3(i,j,:)=0;
           % disp([num2str(i), ' ', num2str(j)]);
        end
    end
end

% 定义要保存的文件名
filename = 'K:\Liu_drought_legacy\phenology_GIMMS_NDVI\NDVI_NH_soiltmp_sgf_effective.mat';
% 使用 save 函数保存矩阵到.mat文件
save(filename, 'var3');


NDVI=importdata('0427NDVI_NH_soiltmp_sgf_effective.mat');

NDVImax=zeros(720,4320,34);
POS=zeros(720,4320,34);
SOS_20=zeros(720,4320,34);
SOS_30=zeros(720,4320,34);
SOS_50=zeros(720,4320,34);
SOS_90=zeros(720,4320,34);

EOS_20=zeros(720,4320,34);
EOS_30=zeros(720,4320,34);
EOS_50=zeros(720,4320,34);
EOS_90=zeros(720,4320,34);

LGS_20_20=zeros(720,4320,34);
LGS_20_50=zeros(720,4320,34);
LPS_90_90=zeros(720,4320,34);

green_rate_20_90=zeros(720,4320,34);
green_rate_50_90=zeros(720,4320,34);
brown_rate_90_50=zeros(720,4320,34);
brown_rate_90_20=zeros(720,4320,34);
GS_90_90=zeros(720,4320,34);
SS_90_90=zeros(720,4320,34);
GS_20_90=zeros(720,4320,34);
SS_90_50=zeros(720,4320,34);
%第1天到第15天
 NDVI_year=zeros(720,4320,24);


for year=1:34
    NDVI_year=NDVI(:,:,1+24*(year-1):24*year);
    for i=1:720
        for j=1:4320
            if all((NDVI_year(i,j,1) ~= 0)) && all(find(NDVI_year(i,j,:)==max(NDVI_year(i,j,:)))>8)&& all(find(NDVI_year(i,j,:)==max(NDVI_year(i,j,:)))<20)
                DOY = 8:15:24*15+8;
                day=  1:1:365;
                series=double(reshape(NDVI_year(i,j,:),24,1))\10000.0;
                series(25)=series(24);
                daily_NDVI = spline(transpose(DOY),series,day);%时间，数值，差值到每一天
        
                %plot(daily_NDVI ,'.-')
                NDVImax(i,j,year)=max(daily_NDVI(60:300));
                POS_date=find(daily_NDVI==max(daily_NDVI(60:300)));%最大值限定范围
                POS(i,j,year)=POS_date(1);

                NDVI_min_position1=find(daily_NDVI(1:POS_date(1))==min(daily_NDVI(1:POS_date(1))));
                NDVI_min_1=NDVI_min_position1(1);
                NDVI_min_position2=find(daily_NDVI(POS_date(1)+1:365)==min(daily_NDVI(POS_date(1)+1:365)));
                NDVI_min_2=NDVI_min_position2(1)+POS_date(1);


                threshold_20_SOS=(max(daily_NDVI(60:300))-min(daily_NDVI(1:POS_date(1))))*0.2+min(daily_NDVI(1:POS_date(1)));
                threshold_30_SOS=(max(daily_NDVI(60:300))-min(daily_NDVI(1:POS_date(1))))*0.3+min(daily_NDVI(1:POS_date(1)));
                threshold_50_SOS=(max(daily_NDVI(60:300))-min(daily_NDVI(1:POS_date(1))))*0.5+min(daily_NDVI(1:POS_date(1)));
                threshold_90_SOS=(max(daily_NDVI(60:300))-min(daily_NDVI(1:POS_date(1))))*0.9+min(daily_NDVI(1:POS_date(1)));

                threshold_20_EOS=(max(daily_NDVI(60:300))-min(daily_NDVI(POS_date(1):365)))*0.2+min(daily_NDVI(POS_date(1):365));
                threshold_30_EOS=(max(daily_NDVI(60:300))-min(daily_NDVI(POS_date(1):365)))*0.3+min(daily_NDVI(POS_date(1):365));
                threshold_50_EOS=(max(daily_NDVI(60:300))-min(daily_NDVI(POS_date(1):365)))*0.5+min(daily_NDVI(POS_date(1):365));
                threshold_90_EOS=(max(daily_NDVI(60:300))-min(daily_NDVI(POS_date(1):365)))*0.9+min(daily_NDVI(POS_date(1):365));


                if all(POS_date(1) >60)  &&  all(POS_date(1) <300)
                    idx_ratio1=find(daily_NDVI(NDVI_min_1:POS_date(1))>=threshold_20_SOS);
                    SOS_20(i,j,year)=idx_ratio1(1)+NDVI_min_1-1;
                    %if (idx_ratio1(1)+48) <= 0 %  确定是属于上升趋势的曲线上取得数值
                    %SOS_20(i,j,year)=idx_ratio1(1)+NDVI_min_1;
                    %end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%防止出现小的波动在夏季******算出来的是最后的EOS如果是双峰或三峰只看最后一个下降
                    idx_ratio2=find(daily_NDVI(POS_date(1)+1:NDVI_min_2)<=threshold_20_EOS);
                    EOS_20(i,j,year)=idx_ratio2(1)+POS_date(1);

                    idx_ratio1=find(daily_NDVI(NDVI_min_1:POS_date(1))>=threshold_30_SOS);
                    SOS_30(i,j,year)=idx_ratio1(1)+NDVI_min_1-1;

                    idx_ratio2=find(daily_NDVI(POS_date(1)+1:NDVI_min_2)<=threshold_30_EOS);
                    EOS_30(i,j,year)=idx_ratio2(1)+POS_date(1);


                    idx_ratio1=find(daily_NDVI(NDVI_min_1:POS_date(1))>=threshold_50_SOS);
                    SOS_50(i,j,year)=idx_ratio1(1)+NDVI_min_1-1;

                    idx_ratio2=find(daily_NDVI(POS_date(1)+1:NDVI_min_2)>=threshold_50_EOS);
                    EOS_50(i,j,year)=idx_ratio2(size(idx_ratio2,2))+POS_date(1);


                    idx_ratio1=find(daily_NDVI(NDVI_min_1:POS_date(1))>=threshold_90_SOS);
                    SOS_90(i,j,year)=idx_ratio1(1)+NDVI_min_1-1;

                    idx_ratio2=find(daily_NDVI(POS_date(1)+1:NDVI_min_2)>=threshold_90_EOS);
                    EOS_90(i,j,year)=idx_ratio2(size(idx_ratio2,2))+POS_date(1);

                    %disp([num2str(SOS_20(i,j,year)),' ',num2str(SOS_30(i,j,year)), ' ',num2str(SOS_50(i,j,year)),' ',num2str(SOS_90(i,j,year)),' ',num2str(EOS_90(i,j,year)),' ',num2str(EOS_50(i,j,year)),' ',num2str(EOS_30(i,j,year)),' ',num2str(EOS_20(i,j,year))]);

                    LGS_20_20(i,j,year)=EOS_20(i,j,year)-SOS_20(i,j,year)+1;
                    LGS_20_50(i,j,year)=EOS_50(i,j,year)-SOS_20(i,j,year)+1;
                    LPS_90_90(i,j,year)=EOS_90(i,j,year)-SOS_90(i,j,year)+1;
                    GS_20_90(i,j,year)=SOS_90(i,j,year)-SOS_20(i,j,year)+1;
                    SS_90_50(i,j,year)=EOS_50(i,j,year)-EOS_90(i,j,year)+1;

                    green_rate_20_90(i,j,year)=(daily_NDVI(SOS_90(i,j,year))-daily_NDVI(SOS_20(i,j,year)))\(SOS_90(i,j,year)-SOS_20(i,j,year)+1);
                    green_rate_50_90(i,j,year)=(daily_NDVI(SOS_90(i,j,year))-daily_NDVI(SOS_20(i,j,year)))\(SOS_90(i,j,year)-SOS_50(i,j,year)+1);
                    brown_rate_90_50(i,j,year)=(daily_NDVI(EOS_90(i,j,year))-daily_NDVI(EOS_50(i,j,year)))\(EOS_50(i,j,year)-EOS_90(i,j,year)+1);
                    brown_rate_90_20(i,j,year)=(daily_NDVI(EOS_90(i,j,year))-daily_NDVI(EOS_20(i,j,year)))\(EOS_20(i,j,year)-EOS_90(i,j,year)+1);


                    if LGS_20_20(i,j,year) < 0
                        LGS_20_20(i,j,year)=0;
                    end

                    if LGS_20_50(i,j,year) < 0
                        LGS_20_50(i,j,year)=0;
                    end
                    if LPS_90_90 < 0
                        LPS_90_90(i,j,year)=0;
                    end

                    if GS_20_90 < 0
                        GS_20_90(i,j,year)=0;
                    end

                    if  SS_90_50 < 0
                        SS_90_50(i,j,year)=0;
                    end

                    if green_rate_20_90(i,j,year) < 0 || green_rate_20_90(i,j,year) > 10000
                        green_rate_20_90(i,j,year)=0;
                    end
                    if green_rate_50_90(i,j,year) < 0 || green_rate_50_90(i,j,year) > 10000
                        green_rate_50_90(i,j,year)=0;
                    end
                    if brown_rate_90_50(i,j,year) < 0 || brown_rate_90_50(i,j,year) > 10000
                        brown_rate_90_50(i,j,year)=0;
                    end
                    if brown_rate_90_20(i,j,year) < 0 || brown_rate_90_20(i,j,year) > 10000
                        brown_rate_90_20(i,j,year)=0;
                    end
                end
            end
        end
        disp([num2str(year), ' ', num2str(i)]);
    end
    [A,R]= readgeoraster('K:\Liu_drought_legacy\phenology_GIMMS_NDVI\soiltmp_backgroud_ndvi_0426.tif');
    suffix='.tif';
    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_20_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, SOS_20(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_30_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, SOS_30(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_50_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, SOS_50(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SOS_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, SOS_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\EOS_20_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, EOS_20(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\EOS_30_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, EOS_30(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\EOS_50_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, EOS_50(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\EOS_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, EOS_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\LGS_20_20_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, LGS_20_20(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\LGS_20_50_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, LGS_20_50(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\LPS_90_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, LPS_90_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\GS_20_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, GS_20_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\SS_90_50_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, SS_90_50(:,:,year), R);


    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\green_rate_20_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, green_rate_20_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\green_rate_50_90_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, green_rate_50_90(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\brown_rate_90_50_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, brown_rate_90_50(:,:,year), R);

    prefix='K:\Liu_drought_legacy\phenology_GIMMS_NDVI\brown_rate_90_20_';%change
    newtif=[prefix,num2str(year+1981),suffix];
    geotiffwrite(newtif, brown_rate_90_20(:,:,year), R);


end


