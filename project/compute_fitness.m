function [obj_val, each_t, each_e, each_gbest_aoi, aoi_loc_best,t_last_task,g_obj_aoi,g_obj_lastt,g_obj_energy,loc_aoi,loc_lastt] = compute_fitness(particle, tau, alltask, dealnum,f_loc,p_trans,T,tbl_num,channelgain,t_local_total,t_last_task,g_obj_aoi,g_obj_lastt,g_obj_energy,loc_aoi,loc_lastt,owega,e_local_dj)
d=tau(:,2);
d_bit=d*1024*1024*8;
owega=tau(:,3);
f_loc=tau(:,5);
p_trans=tau(:,6);
gamma=0.5;
MD_tau_num=1;
tau_num=size(tau,1);
beta = particle(1:tau_num);
r = particle(tau_num+1:2*tau_num);
f = particle(2*tau_num+1:3*tau_num);
wt = particle(3*tau_num+1:4*tau_num);
aoi_loc_now=zeros(1,50);
each_t=zeros(1,tau_num);
each_e=zeros(1,tau_num);
each_gbest_aoi=zeros(1,50);
aoi_loc_best=zeros(1,50);
for i=1:tau_num
    spx=0;
    ssj=0;         
    mdnum=tbl_num(i);
    j=dealnum(mdnum);
    if j>1
        spx=(t_local_total(MD_tau_num*(mdnum-1)+j-1)+wt(i))*t_local_total(MD_tau_num*(mdnum-1)+j);
    end
    try
        if j==MD_tau_num
            ssj=t_local_total(MD_tau_num*(mdnum-1)+j)*t_local_total(MD_tau_num*(mdnum-1)+j)/2;
        else
            ssj=(t_local_total(MD_tau_num*(mdnum-1)+j)+wt(i))*(t_local_total(MD_tau_num*(mdnum-1)+j)+wt(i))/2;
        end
    catch 
        disp(size(t_local_total));
    end
    aoi_loc_now(mdnum)=(loc_aoi(mdnum)*(loc_lastt(mdnum))+ssj+spx)/(loc_lastt(mdnum)+t_local_total(MD_tau_num*(mdnum-1)+j)+wt(i));
end
particlesize=40;
rate=zeros(tau_num,particlesize);
pc=p_trans.*channelgain;
backin=10^(-13); 
for i=1:tau_num
        rate(i)= r(i)*1024*1024*8*log2(1+pc(i)/(backin+sum(pc.*beta)-pc(i)));
end
feaflag=1;
t_exe_loc=(1-beta).*owega.*d_bit./(f_loc*10^9);
e_exe=5*10^(-27).*(f_loc).^3.*t_exe_loc;
t_trans=beta.*d_bit./rate;
e_trans=p_trans.*t_trans;
t_exe_ES=beta.*owega.*d_bit./(f*10^9);
t_total=t_exe_loc+(t_trans+t_exe_ES);
e_total=beta.*e_trans+(1-beta).*e_exe;


for i=1:tau_num
    if t_total(i)>T(i)
        feaflag=0;
        break;
    end
end
if feaflag==0
    obj_val=-1;
else
    aoi_MD_now=zeros(1,50);
    e_MEC=zeros(1,50);
    for i=1:size(tau,1)
        mdnum=tbl_num(i);
        if dealnum(mdnum)==MD_tau_num
            slsj=t_total(i)*t_total(i)/2;
            spx=t_last_task(mdnum)*t_total(i);
            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+slsj)/(g_obj_lastt(mdnum)+wt(i)+t_total(i));
            aoi_MD_now(mdnum)=aoi_new;
            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
        elseif dealnum(mdnum)==1
            ssj=(t_total(i)+wt(i))*(t_total(i)+wt(i))/2;
            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+ssj)/(g_obj_lastt(mdnum)+wt(i)+t_total(i));
            aoi_MD_now(mdnum)=aoi_new;
            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
        else
            ssj=(t_total(i)+wt(i))*(t_total(i)+wt(i))/2;    
            spx=t_last_task(mdnum)*t_total(i);
            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+ssj)/(g_obj_lastt(mdnum)+wt(i)+t_total(i));
            aoi_MD_now(mdnum)=aoi_new;
            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
        end
    end      
    sum_aoi_loc=0;
    sum_e_loc=0;
    sum_aoi_MEC=0;
    sum_e_MEC=0;
    for i=1:tau_num
        sum_aoi_loc=sum_aoi_loc+aoi_loc_now(tbl_num(i));
        sum_e_loc=sum_e_loc+e_local_dj(tbl_num(i),dealnum(tbl_num(i)));
        sum_aoi_MEC=sum_aoi_MEC+aoi_MD_now(tbl_num(i));
        sum_e_MEC=sum_e_MEC+e_MEC(tbl_num(i));
    end
    obj_val=gamma*((sum_aoi_loc-sum_aoi_MEC)/sum_aoi_loc)+(1-gamma)*((sum_e_loc-sum_e_MEC)/sum_e_loc);
    each_t=t_total+wt;
    each_e=e_total;
    each_gbest_aoi=aoi_MD_now;
    aoi_loc_best=aoi_loc_now;
end
end