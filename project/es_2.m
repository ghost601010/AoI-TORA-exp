% % 初始化
% clear;
% clc;
function[alltask1,md_tau,tau,dealnum,tbl_num,d,T,t_remain]=es_2()
MD_num=25;%移动端用户数目 N设为
MD_tau_num=5;%每个移动端用户总共产生的任务数 设为
for a=1:1
    %1区域任务读取
    alltask1=xlsread('param','Sheet5','A2:I126');
    for b=1:MD_num
        md_tau(b,:,:)=alltask1(MD_tau_num*(b-1)+1:MD_tau_num*b,:);
    end
    %{
    md1_tau=alltask1(MD_tau_num*0+1:MD_tau_num*1,:);
    md2_tau=alltask1(MD_tau_num*1+1:MD_tau_num*2,:);
    md3_tau=alltask1(MD_tau_num*2+1:MD_tau_num*3,:);
    md4_tau=alltask1(MD_tau_num*3+1:MD_tau_num*4,:);
    md5_tau=alltask1(MD_tau_num*4+1:MD_tau_num*5,:);
    %}
    dealnum=ones(1,MD_num);
    for b=1:MD_num
        tau(b,:)=md_tau(b,dealnum(b),:);
    end
    %s=xlsread('param','Sheet1','N2:O7');
    d=tau(:,2);%任务数据量(MB)
    %s_low=tau(:,4);%任务需要达到的最低安全等级
    T=tau(:,7);%截止时间
    tbl_num=tau(:,9);%任务所属的移动端序号
    %s_time=s(:,1);%每个安全等级处理单位数据需要的时间
    %s_energy=s(:,2);%每个安全等级处理单位数据需要的时间
    %B=50;%带宽
end
end

