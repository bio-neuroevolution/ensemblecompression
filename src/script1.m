%测试百分位数和外力因子对压缩率的影响
%% 初始化路径
datasetname = ['车辆','船舶'];
paths = ["/Volumes/KODAK/data/Taxi_070220/","/Volumes/KODAK/data/tianchi/VIS/hy_round2_train_20200225/"];
varNames = [["ID","time","x","y","speed","direction","type"];["ID", "x", "y", "VarName4", "VarName5", "time", "type"]];
varTypes = [["double","string","double","double","double","double","double"];["double", "double", "double", "double", "double", "string", "categorical"]];
ranges  = [[2, Inf];[2, Inf]];
files = [["Taxi_105"];["20003.csv"]];  %文件集，[]则为所有文件

% 算法参数
E = 2;
Po = 0.2;
density = 1.15;
maxit = 100;
tol = 1e-2;
thick = 1;

percentage = [25,35,45,55, 65, 75, 85,95];
force = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];

[p0,p1]=meshgrid(percentage,force);
warning off;

% 遍历每一个数据集
for i = 1:length(paths)
    path = paths(i);
    varname = varNames(i,:);
    vartype = varTypes(i,:);
    range = ranges(i,:);
    file = files(i,:);
    file = strcat(path,file);
    % 得到数据集里所有文件名
    if isempty(file)
        dirOutput=dir(fullfile(path,'*'));
        file = {dirOutput.name};
    end
    
    % 初始化指标
    indictors = zeros([size(p0),5]);
    count = 1;
    % 遍历每一个数据文件
    for j = 1:size(file)
        % 读取数据
        [x0,y0] = importfile(file(j),range,varname,vartype);
        s = size(x0);
        orgin_len = sum(sqrt(diff(x0).^2 + diff(y0).^2));
        orgin_angle = 0;
        for t = 2:size(x0)-1
            t0x = x0(t-1);t1x = x0(t);t2x=x0(t+1);
            t0y = y0(t-1);t1y = y0(t);t2y=y0(t+1);
            if ((t0x == t1x && t0y==t1y) || (t1x==t2x && t1y==t2y))
                continue;
            end
            tt = (norm([t0x-t1x,t0y-t1y])*norm([t0x-t1x,t2y-t1y]));
            if tt ~= 0
                orgin_angle = orgin_angle + dot([t0x-t1x,t0y-t1y],[t2x-t1x,t2y-t1y])/tt;
            end
        end                             
        
        % 对每一个参数,调用算法获得结果
        [r,c] = size(p0);
        for k1 = 1:r
            for k2 = 1:c
                begintime = clock;
                [x1,y1,ss,Range,time] = compress_finite_element(x0,y0,s,E,Po,density,p0(k1,k2),maxit,tol,thick,p1(k1,k2));
                spantime = etime(clock,begintime);
                indictors(k1,k2,1) = indictors(k1,k2,1) + length(x0)*56/spantime; %花费时间
                indictors(k1,k2,2) = indictors(k1,k2,2) + 1 - length(x1) / length(x0);  %压缩率
                indictors(k1,k2,3) = indictors(k1,k2,3) + sum(sqrt(diff(x1).^2 + diff(y1).^2))/orgin_len;   %长度比       
                
                if indictors(k1,k2,2)>=0.2 && indictors(k1,k2,2)<=0.6
                    save(strcat(path,"compress",num2str(indictors(k1,k2,2)),".dat"),'x0','y0','x1','y1');
                end
                
                compress_angle = 0;
                for t = 2:size(x1')-1
                    t0x = x1(t-1);t1x = x1(t);t2x=x1(t+1);
                    t0y = y1(t-1);t1y = y1(t);t2y=y1(t+1);
                    tt = (norm([t0x-t1x,t0y-t1y])*norm([t0x-t1x,t2y-t1y]));
                    if tt ~= 0
                        compress_angle = compress_angle + dot([t0x-t1x,t0y-t1y],[t2x-t1x,t2y-t1y])/tt;
                    end
                end  
                indictors(k1,k2,4) = indictors(k1,k2,4) + compress_angle/orgin_angle; %曲度比
                [d1,~,~] = dtw([x0';y0'],[x1;y1]);
                indictors(k1,k2,5) = indictors(k1,k2,5) + d1;   %误差率   
            end
        end  
    end
    
    
    % 计算指标平均值
    indictors = indictors ./ length(file);
    % 误差率,长度比,曲度比归一化
    indictors(:,:,3) = mapminmax(indictors(:,:,3),0,1);
    indictors(:,:,4) = mapminmax(indictors(:,:,4),0,1);
    indictors(:,:,5) = mapminmax(indictors(:,:,5),0,1);
    % 打印输出
    fprintf('数据集%s处理完毕.\n',path);
    fprintf('\t压缩速率=%fbps\n',mean2(indictors(:,:,1)));
    fprintf('\t压缩率=%f \n',mean2(indictors(:,:,2)));
    fprintf('\t长度比=%f \n',mean2(indictors(:,:,3)));
    fprintf('\t曲度比=%f\n',mean2(indictors(:,:,4)));
    fprintf('\t相对误差率=%f\n',mean2(indictors(:,:,5)));
    
    save(strcat(path,"results.dat"),'p0','p1','indictors');

    subplot(1,2,i)
    mesh(p0,p1,indictors(:,:,5));
    xlabel("percentile"); 
    ylabel("stress factor");
    colormap Parula;
    hold off;
    
    

end