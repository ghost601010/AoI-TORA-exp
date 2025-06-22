global alltask1;
%global tau;
global p;
global alltask2;
%global tauu;
%global s;
global f_all;
global B_all;
f_all=randi([35,45]);
B_all=9.97*randi([5,10]);
global f_remain;
global B_remain;
global gamma;
gamma=0.5;
global decision;

ES_num=2;
MD_num=25;
MD_tau_num=5;
global t_last_task;
global g_obj_aoi;
global g_obj_lastt;
global g_obj_energy;
global loc_aoi;
global loc_lastt;

global final_t1;
global final_e1;
global final_t2;
global final_e2;
tnum=1;
result1=zeros(tnum,1);
result2=zeros(tnum,1);
result3=zeros(tnum,1);
op_obj=zeros(1,120);
len_cloud=zeros(1,120);
wtmax=0.8;
wtmin=0;
global a;
vwtmax=wtmax;
vwtmin=-wtmax;
MG = 1000; 
particlesize = 50; 
c1 = 1.49618;
c2 = 1.49618; 
g_obj_cloud=0;

for a=1:tnum
    clc;
    loc_aoi=zeros(1,MD_num);
    loc_lastt=zeros(1,MD_num);
    t_last_task=zeros(1,MD_num);
    g_obj_aoi=zeros(1,MD_num);
    g_obj_lastt=zeros(1,MD_num);
    g_obj_energy=zeros(1,MD_num);
    final_t1=zeros(MD_num*MD_tau_num,1);
    final_e1=zeros(MD_num*MD_tau_num,1);
    final_t2=zeros(MD_num*MD_tau_num,1);
    final_e2=zeros(MD_num*MD_tau_num,1);
    decision=0;
    [alltask1,md1_tau,tau,dealnum1]=es_1();
    [alltask2,md2_tau,tauu,dealnum2]=es_2();
    d1=tau(:,2);
    tbl_num1=tau(:,9);
    t_remain1=tau(:,8);
    d2=tauu(:,2);
    T2=tauu(:,7);
    t_remain2=tauu(:,8);
    d1_bits=d1*1024*1024*8;
    d2_bits=d2*1024*1024*8;
    owega1=tau(:,3);
    p_trans1=tau(:,6);
    tbl_num2=tauu(:,9);
    owega2=tauu(:,3);
    p_trans2=tauu(:,6);
    r_cloud_1=99.7;
    r_cloud_2=99.7;
    f_cloud=18*10^9;
    cloud_sqe=[];
    cloud_num=0;
    
    alld1=alltask1(:,2);
    all_d_bit1=alld1*1024*1024*8;
    all_f_loc1=alltask1(:,5);
    all_owega1=alltask1(:,3);
    t_local_total1=all_d_bit1.*all_owega1./(all_f_loc1*10^9);
    e_loc_total1=5*10^(-27).*(all_f_loc1*10^9).^3.*t_local_total1;
    
    alld2=alltask2(:,2);
    all_d_bit2=alld2*1024*1024*8;
    all_f_loc2=alltask2(:,5);
    all_owega2=alltask2(:,3);
    t_local_total2=all_d_bit2.*all_owega2./(all_f_loc2*10^9);
    e_loc_total2=5*10^(-27).*(all_f_loc2*10^9).^3.*t_local_total2;
    
    cloud_t=0;
    cloud_e=0;
    cloud_gbest_aoi=g_obj_aoi;
    cloud_obj=zeros(MD_num,MD_tau_num);
    firstflag1=1;
    firstflag2=1;
    t0=cputime;
    tic;
    while sum(dealnum1(:)) ~= (MD_tau_num+1)*MD_num
        decision=decision+1;
        if firstflag1==1 && sum(dealnum1(:))==2*MD_num
            firstflag1=0;
        end
        if firstflag1~=1
            for h=1:MD_num
                tau(h,:)=md1_tau(h,dealnum1(h),:);
            end
        end
        tau_num=size(tau,1);
        [gest,each_t,each_e,each_gbest_aoi,each_obj,g_obj,e_local_dj,aoi_loc_best,T,dealnum1]=pso(dealnum1,alltask1,tau,f_all,B_all,0,MD_num);
        [tau,dealnum1,g_obj,final_t1,final_e1,redoflag]=cloud_pso(d1,dealnum1,tau,e_local_dj,g_obj,each_t,each_e,each_gbest_aoi,aoi_loc_best,gest,T,final_t1,final_e1,MD_num);
        if redoflag>0
            tau_num=size(tau,1);
            [gest,each_t,each_e,each_gbest_aoi,each_obj,g_obj,e_local_dj,aoi_loc_best,T,dealnum1]=pso(dealnum1,alltask1,tau,f_all,B_all,1,tau_num);
        end
    end 
    t1=cputime-t0;
    disp(t1);
    sum_aoi=sum(g_obj_aoi);
    sum_energy=sum(g_obj_energy);
    sum_loc_aoi=sum(loc_aoi);
    sum_loc_energy=sum(e_local_dj(:,MD_tau_num));
    obj_final1=gamma*((-sum_aoi+sum_loc_aoi)/sum_loc_aoi)+(1-gamma)*((-sum_energy+sum_loc_energy)/sum_loc_energy);
    result1(a)=obj_final1;
    
    disp(loc_aoi);
    disp(g_obj_aoi);
    decision=0;
    loc_aoi=zeros(1,MD_num);
    loc_lastt=zeros(1,MD_num);
    t_last_task=zeros(1,MD_num);
    g_obj_aoi=zeros(1,MD_num);
    g_obj_lastt=zeros(1,MD_num);
    g_obj_energy=zeros(1,MD_num);
    while sum(dealnum2(:)) ~= (MD_tau_num+1)*MD_num
        decision=decision+1;
        if firstflag2==1 && sum(dealnum2(:))==2*MD_num
            firstflag2=0;
        end
        if firstflag2~=1
            for h=1:MD_num
                tauu(h,:)=md2_tau(h,dealnum2(h),:);
            end
        end
            tauu_num=size(tauu,1);
            [gest2,each_t2,each_e2,each_gbest_aoi2,each_obj2,g_obj2,e_local_dj2,aoi_loc_best2,T,dealnum2]=pso(dealnum2,alltask2,tauu,f_all,B_all,0,MD_num);
            [tauu,dealnum2,g_obj2,final_t2,final_e2,redoflag]=cloud_pso(d2,dealnum2,tauu,e_local_dj2,g_obj2,each_t2,each_e2,each_gbest_aoi2,aoi_loc_best2,gest2,T,final_t2,final_e2,MD_num);
            if redoflag>0
                tauu_num=size(tauu,1);
                [gest2,each_t2,each_e2,each_gbest_aoi2,each_obj2,g_obj2,e_local_dj2,aoi_loc_best2,T,dealnum2]=pso(dealnum2,alltask2,tauu,f_all,B_all,1,tauu_num,e_local_dj2);
            end
        disp(['a=',num2str(a)]);
        disp(['ES2中第',num2str(decision),'次任务卸载决策后任务处理的平均AOI和能耗加权和优化情况如下：']);
        disp(g_obj2);
        disp(g_obj_aoi);
    end
    disp(loc_aoi);
    disp(g_obj_aoi);
    sum_aoi=sum(g_obj_aoi);
    sum_energy=sum(g_obj_energy);
    sum_loc_aoi=sum(loc_aoi);
    sum_loc_energy=sum(e_local_dj2(:,MD_tau_num));
    obj_final2=gamma*((-sum_aoi+sum_loc_aoi)/sum_loc_aoi)+(1-gamma)*((-sum_energy+sum_loc_energy)/sum_loc_energy);
    result2(a)=obj_final2;
end
toc;
