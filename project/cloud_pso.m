function[tau,dealnum,g_obj,final_t,final_e,redoflag]=cloud_pso(d,dealnum,tau,e_local_dj,g_obj,each_t,each_e,each_gbest_aoi,aoi_loc_best,gest,T,final_t,final_e,aa)
global t_last_task;
global g_obj_aoi;
global g_obj_lastt;
global g_obj_energy;
global gamma;
global loc_aoi;
global loc_lastt;
global t_local_total;
global decision;
MD_num=aa;
MD_tau_num=5;
MG = 800;  
particlesize = 40;
c1 = 2.05; 
c2 = 2.05;
wmax = 0.9; 
wmin=0.4;
timeout=0;
sourceout=0;
wtmax=0.5;
wtmin=0;
vwtmax=wtmax*randi([10,20])/100;
vwtmin=-vwtmax;
g_obj_cloud=0;

tau_num=size(tau,1);
d_bits=d*1024*1024*8;
tbl_num=tau(:,9);
owega=tau(:,3);
p_trans=tau(:,6);
r_cloud_1=1.52;
r_cloud_2=1.52;
f_cloud=10*10^9;
cloud_sqe=[];
cloud_num=0;
cloud_t=0;
cloud_e=0;
cloud_gbest_aoi=g_obj_aoi;
cloud_obj=zeros(MD_num,MD_tau_num);
loc_aoi_new=0;
redoflag=0;
certainflag=zeros(tau_num,1);
    for i = 1:tau_num
        g_obj_cloud=0;
        mdnum=tbl_num(i);
        p=zeros(1,particlesize);
        p(1,:)=rand(1,particlesize);
        p(1,:)=wtmax.*p(1,:);
        obj_cloud=zeros(1,particlesize)
        p_obj_cloud=obj_cloud;
        pest_cloud=p;
        gest_cloud=p(1,randi([1 particlesize]));
        for t=1:particlesize
            spx=0;
            ssj=0;
            j=dealnum(mdnum);
            if j>1
                spx=(t_local_total(MD_tau_num*(mdnum-1)+j-1)+p(1,t))*t_local_total(MD_tau_num*(mdnum-1)+j);
            end
            if j==MD_tau_num
                ssj=t_local_total(MD_tau_num*(mdnum-1)+j)*t_local_total(MD_tau_num*(mdnum-1)+j)/2;
            else
                ssj=(t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t))*(t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t))/2;
            end
            aoi_loc_now=(loc_aoi(mdnum)*(loc_lastt(mdnum))+ssj+spx)/(loc_lastt(mdnum)+t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t));
               
            feaflag=1;
            t_trans_cloud=d(i)/(r_cloud_1*10^6/(8*1024*1024))+d(i)/(r_cloud_2*10^6/(8*1024*1024));
            e_trans_cloud=p_trans(i)*t_trans_cloud;
            t_exe_cloud=d_bits(i)*owega(i)/f_cloud;
            t_total=t_trans_cloud+t_exe_cloud;
            
            if t_total>T(i)
                feaflag=0;
            end
            if feaflag==0
                timeout=timeout+1;
                obj_cloud(t)=-1;
            else
                if dealnum(mdnum)==MD_tau_num
                    slsj=t_total*t_total/2;
                    spx=t_last_task(mdnum)*t_total;
                    aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+slsj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                    e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                elseif dealnum(mdnum)==1
                    ssj=(t_total+p(1,t))*(t_total+p(1,t))/2;
                    aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+ssj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                    e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                else
                    ssj=(t_total+p(1,t))*(t_total+p(1,t))/2;
                    spx=t_last_task(mdnum)*t_total;
                    aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+ssj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                    e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                end
                sum_aoi_loc=0;
                sum_e_loc=0;
                sum_aoi_cloud=0;
                sum_e_cloud=0;
                for h=1:tau_num
                    if h==i
                        sum_aoi_loc=sum_aoi_loc+aoi_loc_now;
                        sum_aoi_cloud=sum_aoi_cloud+aoi_new;
                        sum_e_cloud=sum_e_cloud+e_new;
                    else
                        sum_aoi_loc=sum_aoi_loc+aoi_loc_best(tbl_num(h));
                        sum_aoi_cloud=sum_aoi_cloud+each_gbest_aoi(tbl_num(h));
                        sum_e_cloud=sum_e_cloud+each_e(h);
                    end
                    sum_e_loc=sum_e_loc+e_local_dj(tbl_num(h),decision);
                end
                obj_cloud(t)=gamma*((sum_aoi_loc-sum_aoi_cloud)/sum_aoi_loc)+(1-gamma)*((sum_e_loc-sum_e_cloud)/sum_e_loc);
            end
            p_obj_cloud(t)=obj_cloud(t);
            if obj_cloud(t)>g_obj_cloud
                g_obj_cloud=obj_cloud(t);
                gest_cloud=p(:,t);
                cloud_t=t_total+p(1,t);
                cloud_e=e_trans_cloud;
                cloud_gbest_aoi(mdnum)=aoi_new;
                loc_aoi_new=aoi_loc_now;
            end
        end
        
        v = zeros(1,particlesize);
        v(1,:) = rand(1,particlesize)*(wtmax-wtmin)-p(1,:);
        g=1;
        while g<=MG
            for t =1:particlesize
                v(1,t)=0.729*v(1,t)+c1*rand*sign(pest_cloud(1,t)-p(1,t))*0.2+c2*rand*sign(gest_cloud-p(1,t))*0.2;
                for j = 1:1
                    if v(j,t)>vwtmax
                        v(j,t)=vwtmax;
                    elseif v(j,t)<vwtmin
                        v(j,t)=vwtmin;
                    end
                end
                p(:,t)=p(:,t)+v(:,t);
                for j = 1:1
                    if p(j,t)>wtmax
                        p(j,t)=wtmax;
                    elseif p(j,t)<wtmin
                        p(j,t)=wtmin;
                    end
                end
                
                spx=0;
                ssj=0;
                j=dealnum(mdnum);
                if j>1
                    spx=(t_local_total(MD_tau_num*(mdnum-1)+j-1)+p(1,t))*t_local_total(MD_tau_num*(mdnum-1)+j);
                end
                if j==MD_tau_num
                    ssj=t_local_total(MD_tau_num*(mdnum-1)+j)*t_local_total(MD_tau_num*(mdnum-1)+j)/2;
                else
                    ssj=(t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t))*(t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t))/2;
                end
                aoi_loc_now=(loc_aoi(mdnum)*(loc_lastt(mdnum))+ssj+spx)/(loc_lastt(mdnum)+t_local_total(MD_tau_num*(mdnum-1)+j)+p(1,t));
                
                feaflag=0;
                t_trans_cloud=d(i)/(r_cloud_1*10^6/(8*1024*1024))+d(i)/(r_cloud_2*10^6/(8*1024*1024));
                e_trans_cloud=p_trans(i)*t_trans_cloud;
                t_exe_cloud=d_bits(i)*owega(i)/f_cloud;
                t_total=t_trans_cloud+t_exe_cloud;
                
                if t_total>T(i)
                    feaflag=0;
                end
                if feaflag==0
                    timeout=timeout+1;
                    obj_cloud(t)=-1;
                else
                    if dealnum(mdnum)==MD_tau_num
                        slsj=t_total*t_total/2;
                        spx=t_last_task(mdnum)*t_total;
                        aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+slsj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                        e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                    elseif dealnum(mdnum)==1
                        ssj=(t_total+p(1,t))*(t_total+p(1,t))/2;
                        aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+ssj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                        e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                    else
                        ssj=(t_total+p(1,t))*(t_total+p(1,t))/2;
                        spx=t_last_task(mdnum)*t_total;
                        aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+ssj)/(g_obj_lastt(mdnum)+t_total+p(1,t));
                        e_new=g_obj_energy(mdnum)+ e_trans_cloud;
                    end
                    sum_aoi_loc=0;
                    sum_e_loc=0;
                    sum_aoi_cloud=0;
                    sum_e_cloud=0;
                    for h=1:tau_num
                        if h==i
                            sum_aoi_loc=sum_aoi_loc+aoi_loc_now;
                            sum_aoi_cloud=sum_aoi_cloud+aoi_new;
                            sum_e_cloud=sum_e_cloud+e_new;
                        else
                            sum_aoi_loc=sum_aoi_loc+aoi_loc_best(tbl_num(h));
                            sum_aoi_cloud=sum_aoi_cloud+each_gbest_aoi(tbl_num(h));
                            sum_e_cloud=sum_e_cloud+each_e(h);
                        end
                        sum_e_loc=sum_e_loc+e_local_dj(tbl_num(h),decision);
                    end
                    obj_cloud(t)=gamma*((sum_aoi_loc-sum_aoi_cloud)/sum_aoi_loc)+(1-gamma)*((sum_e_loc-sum_e_cloud)/sum_e_loc);
                end
                if obj_cloud(t)>p_obj_cloud(t)
                    pest_cloud(:,t)=p(:,t);
                    p_obj_cloud(t)=obj_cloud(t);
                end
                if obj_cloud(t)>g_obj_cloud
                    g_obj_cloud=obj_cloud(t);
                    gest_cloud=p(1,t);
                    cloud_t=t_total+p(1,t);
                    cloud_e=e_trans_cloud;
                    cloud_gbest_aoi(mdnum)=aoi_new;
                    loc_aoi_new=aoi_loc_now;
                end
            end
            g=g+1;
        end
        
        if g_obj_cloud-g_obj>0.01
            %disp(['本次决策中，移动端设备',num2str(mdnum),'的任务',num2str(dealnum(mdnum)),'选择卸载到云上。下一个任务产生的等待时间为',num2str(gest_cloud),'s']);
            t_last_task(mdnum)=cloud_t;
            g_obj_aoi(mdnum)=cloud_gbest_aoi(mdnum);
            g_obj_lastt(mdnum)= g_obj_lastt(mdnum)+cloud_t;
            g_obj_energy(mdnum)=g_obj_energy(mdnum)+cloud_e;
            loc_aoi(mdnum)=loc_aoi_new;
            loc_lastt(mdnum)=loc_lastt(mdnum)+t_local_total(MD_tau_num*(mdnum-1)+dealnum(mdnum))+gest_cloud;
            final_t(MD_tau_num*(mdnum-1)+dealnum(mdnum))=cloud_t-gest_cloud;
            final_e(MD_tau_num*(mdnum-1)+dealnum(mdnum))=cloud_e;
            dealnum(tbl_num(i))=dealnum(tbl_num(i))+1;
            certainflag(i,1)=1;
            redoflag=redoflag+1;
        end
    end
    
    if redoflag==0
        for i=1:tau_num
            mdnum=tbl_num(i);
            if gest(i)==1
                %disp(['本次决策中，移动端设备',num2str(mdnum),'的任务',num2str(dealnum(mdnum)),'选择卸载到边缘服务器上。下一个任务产生的等待时间为',num2str(gest(3*tau_num+i)),'s']);
                final_t(MD_tau_num*(mdnum-1)+dealnum(mdnum))=each_t(i)-gest(3*tau_num+i);
                final_e(MD_tau_num*(mdnum-1)+dealnum(mdnum))=each_e(i);
            else
                %disp(['本次决策中，移动端设备',num2str(mdnum),'的任务',num2str(dealnum(mdnum)),'选择在本地处理。下一个任务产生的等待时间为',num2str(gest(3*tau_num+i)),'s']);
                final_t(MD_tau_num*(mdnum-1)+dealnum(mdnum))=t_local_total(MD_tau_num*(mdnum-1)+dealnum(mdnum));
                final_e(MD_tau_num*(mdnum-1)+dealnum(mdnum))=each_e(i);
            end
            t_last_task(tbl_num(i))=each_t(i);
            g_obj_aoi(tbl_num(i))=each_gbest_aoi(tbl_num(i));
            g_obj_lastt(tbl_num(i))=g_obj_lastt(tbl_num(i))+t_last_task(tbl_num(i));
            g_obj_energy(tbl_num(i))=g_obj_energy(tbl_num(i))+each_e(i);
            loc_aoi(mdnum)=aoi_loc_best(mdnum);
            loc_lastt(mdnum)=loc_lastt(mdnum)+each_t(i);
            dealnum(tbl_num(i))=dealnum(tbl_num(i))+1;
        end
    else 
        newtau=zeros(tau_num-redoflag,9);
        j=1;
        for i=1:tau_num
            if certainflag(i,1)==0
                newtau(j,:)=tau(i,:);
                j=j+1;
            end
        end
        tau=newtau;
    end         
end
