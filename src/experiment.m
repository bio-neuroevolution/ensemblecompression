function experiment(fig)

datasetname = ['车辆','船舶'];
paths = ["/Volumes/KODAK/data/Taxi_070220/","/Volumes/KODAK/data/tianchi/VIS/hy_round2_train_20200225/"];

if fig == 5  % 不同参数的压缩率
    for i = [1:length(paths)]
        load(strcat(paths(i),"results.dat"),'-mat','p0','p1','indictors');
        
        subplot(1,2,i)
        mesh(p0,p1,indictors(:,:,2));
        xlabel("percentile"); 
        ylabel("stress factor");
        colormap Parula;
        hold off;    
    end
    fprintf('\t压缩速率=%fbps\n',mean2(indictors(:,:,1)));
    fprintf('\t压缩率=%f \n',mean2(indictors(:,:,2)));
    fprintf('\t长度比=%f \n',mean2(indictors(:,:,3)));
    fprintf('\t曲度比=%f\n',mean2(indictors(:,:,4)));
    fprintf('\t相对误差率=%f\n',mean2(indictors(:,:,5)));   
       
elseif fig == 6  % 不同参数的压缩速率
    for i = [1:length(paths)]
        load(strcat(paths(i),"results.dat"),'-mat','p0','p1','indictors');
        
        subplot(1,2,i)
        mesh(p0,p1,indictors(:,:,1));
        xlabel("percentile"); 
        ylabel("stress factor");
        colormap Parula;
        hold off;    
    end

elseif fig == 7  % 不同参数的相对错误率
    for i = [1:length(paths)]
        load(strcat(paths(i),"results.dat"),'-mat','p0','p1','indictors');
        
        subplot(1,2,i)
        mesh(p0,p1,indictors(:,:,5));
        xlabel("percentile"); 
        ylabel("stress factor");
        colormap Parula;
        hold off;    
    end

elseif fig == 8 %压缩率与相对错误率关系曲线图
    for i = [1:length(paths)]
        load(strcat(paths(i),"results.dat"),'-mat','p0','p1','indictors');
        subplot(1,2,i);
        indictors = indictors(6,:,:);
        size(indictors)
        indictors(1,:,2)
        indictors(1,:,5)
        plot(indictors(1,:,2),indictors(1,:,5));
        xlabel('compression ratio');
        ylabel('error ratio');
    end
elseif fig==9  %不同比例尺下的轨迹压缩显示
    load("/Volumes/KODAK/data/Taxi_070220/compress0.2607.dat",'-mat','x0','y0','x1','y1');
    xl = x1;
    yl = y1;
    load("/Volumes/KODAK/data/Taxi_070220/compress0.53677.dat",'-mat','x0','y0','x1','y1');
    xh = x1;
    yh = y1;
    subplot(1,3,1);plot(xh,yh,'-*b');legend('0.53677');
    subplot(1,3,2);plot(xl,yl,'-*b');legend('0.2607');
    subplot(1,3,3);plot(x0,y0,'-*b');legend('raw');
elseif fig==10 %长度比和曲度比
    load(strcat(paths(1),"results.dat"),'-mat','p0','p1','indictors');
    indictors = indictors(6,:,:);
    subplot(1,2,1);
    plot(indictors(1,:,2),indictors(1,:,3),indictors(1,:,2),indictors(1,:,4));
    legend('length ratio','tortuosity ratio');
    xlabel('compression ratio');
    ylabel('length & curvature ratio');
    
    
end