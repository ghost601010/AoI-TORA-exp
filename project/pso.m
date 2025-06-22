function [gest, each_t, each_e, each_gbest_aoi, each_obj, g_obj, e_local_dj, aoi_loc_best, T, dealnum] = pso(dealnum, alltask, tau, f_avail, B_avail, againflag, aa)
global t_last_task g_obj_aoi g_obj_lastt g_obj_energy loc_aoi loc_lastt p ;                           
global gamma;
MG = 1800;
particlesize = 40;
betamin=0;
betamax=1;
fmin=1.5;
fmax=3*f_avail/(size(tau,1)*2);
rme=3*B_avail/(size(tau,1)*2);
rmin=2;
wtmax=0.5;
wtmin=0;
vbetamax=1*randi([20,40])/100;
vbetamin=-vbetamax;
vrmax=rme*randi([20,40])/100;
vrmin=-vrmax;
vfmax=fmax*randi([20,40])/100;
vfmin=-vfmax;
vwtmax=wtmax*randi([20,40])/100;
vwtmin=-vwtmax;


MD_num=aa;
MD_tau_num=5;
w = 0.8; 

tau_num=size(tau,1);
d=tau(:,2);
alld=alltask(:,2);
all_d_bit=alld*1024*1024*8;
tbl_num=tau(:,9);
d_bit=d*1024*1024*8;
owega=tau(:,3);
all_owega=alltask(:,3);
f_loc=tau(:,5);
all_f_loc=alltask(:,5);
p_trans=tau(:,6);
g_obj=0;
distance=(5+5*rand(tau_num,1));
channelgain=distance.^(-4);
backin=10^(-13); 
global t_local_total;
if againflag==0
    t_local_total=all_d_bit.*all_owega./(all_f_loc*10^9);
    e_loc_total=5*10^(-27).*(all_f_loc*10^9).^3.*t_local_total;
end
t_loc_ddl=owega.*d_bit./(f_loc*10^9);
T=t_loc_ddl.*(randi([15,20])/10);
global e_local_dj;
if againflag==0
    e_local_dj=zeros(MD_num,MD_tau_num);
    for i=1:MD_num
        for j=1:MD_tau_num
            if j==1
                e_local_dj(i,j)=e_loc_total(MD_tau_num*(i-1)+j);
            else
                e_local_dj(i,j)=e_local_dj(i,j-1)+e_loc_total(MD_tau_num*(i-1)+j);
            end
        end
    end
end

each_t=zeros(1,MD_num);
each_e=zeros(1,MD_num);
each_gbest_aoi=zeros(1,50);
each_obj=zeros(MD_num,MD_tau_num);
aoi_loc_best=zeros(1,50);
global gamma;
for k=1:1
    t0=cputime;
    p=zeros(tau_num*4,particlesize);
    pp=rand(tau_num*4,particlesize);
    deta=zeros(particlesize,1);
    deta=pp(:,1)-0.5;
    p(:,1)=(betamax+betamin)/2+(betamax-betamin).*sign(deta).*(abs(deta)).^(1+1/particlesize);
    for q=2:particlesize
        pp(:,q)=4.*pp(:,q-1).*(1-pp(:,q-1));
    end
    p(1:tau_num,:)=pp(1:tau_num,:)*(betamax-betamin)+betamin;
    p(1:tau_num,:)=round( p(1:tau_num,:));
    p(tau_num+1:tau_num*2,:)=pp(tau_num+1:tau_num*2,:)*(rme-rmin)+rmin;
    rate=zeros(tau_num,particlesize);
    pc=p_trans.*channelgain;
    for i=1:tau_num
        for u=1:particlesize
            rate(i,u)= p(tau_num+i,u)*1024*1024*8*log2(1+pc(i)/(backin+sum(pc.*p(1:tau_num,u))-pc(i)));
        end
    end
    p(2*tau_num+1:tau_num*3,:)=pp(2*tau_num+1:tau_num*3,:)*(fmax-fmin)+fmin;
    p(3*tau_num+1:tau_num*4,:)=pp(3*tau_num+1:tau_num*4,:)*(wtmax-wtmin)+wtmin;
  
    obj=zeros(1,particlesize)-1000;
    p_obj=obj;
    
    aoi_loc_now=zeros(1,50);
    for t=1:particlesize
        for i=1:tau_num
            spx=0;
            ssj=0;
            mdnum=tbl_num(i);
            j=dealnum(mdnum);
            if j>1
                spx=(t_local_total(MD_tau_num*(mdnum-1)+j-1)+p(3*tau_num+i,t))*t_local_total(MD_tau_num*(mdnum-1)+j);
            end
            if j==MD_tau_num
                ssj=t_local_total(MD_tau_num*(mdnum-1)+j)*t_local_total(MD_tau_num*(mdnum-1)+j)/2;
            else
                ssj=(t_local_total(MD_tau_num*(mdnum-1)+j)+p(3*tau_num+i,t))*(t_local_total(MD_tau_num*(mdnum-1)+j)+p(3*tau_num+i,t))/2;
            end
            aoi_loc_now(mdnum)=(loc_aoi(mdnum)*(loc_lastt(mdnum))+ssj+spx)/(loc_lastt(mdnum)+t_local_total(MD_tau_num*(mdnum-1)+j)+p(3*tau_num+i,t));
        end
    
        feaflag=1;
        t_exe_loc=(1-p(1:tau_num,t)).*owega.*d_bit./(f_loc*10^9);
        %disp(t_exe_loc);
        e_exe=5*10^(-27).*(f_loc).^3.*t_exe_loc;
        t_trans=p(1:tau_num,t).*d_bit./rate(:,t);
        %disp(t_trans);
        e_trans=p_trans.*t_trans;
        t_exe_ES=p(1:tau_num,t).*owega.*d_bit./(p(2*tau_num+1:tau_num*3,t)*10^9);
        %disp(t_exe_ES);
        t_total=t_exe_loc+(t_trans+t_exe_ES);
        e_total=p(1:tau_num,t).*e_trans+(1-p(1:tau_num,t)).*e_exe;
        %disp(T-t_total);
      
        for i=1:tau_num
            if t_total(i)>T(i)
                feaflag=0;
                break;
            end
        end
        if feaflag==0
            obj(t)=-1;
        else
            aoi_MD_now=zeros(1,50);
            e_MEC=zeros(1,50);
            for i=1:size(tau,1)
                mdnum=tbl_num(i);
                if dealnum(mdnum)==MD_tau_num
                    slsj=t_total(i)*t_total(i)/2;
                    spx=t_last_task(mdnum)*t_total(i);
                    aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+slsj)/(g_obj_lastt(mdnum)+t_total(i)+p(3*tau_num+i,t));
                    aoi_MD_now(mdnum)=aoi_new;
                    e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
                elseif dealnum(mdnum)==1
                    ssj=(t_total(i)+p(3*tau_num+i,t))*(t_total(i)+p(3*tau_num+i,t))/2;
                    aoi_new=ssj/(t_total(i)+p(3*tau_num+i,t));
                    aoi_MD_now(mdnum)=aoi_new;
                    e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
                else
                    ssj=(t_total(i)+p(3*tau_num+i,t))*(t_total(i)+p(3*tau_num+i,t))/2;    
                    spx=t_last_task(mdnum)*t_total(i);
                    aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+ssj)/(g_obj_lastt(mdnum)++p(3*tau_num+i,t)+t_total(i));
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
            obj(t)=gamma*((sum_aoi_MEC-sum_aoi_loc)/sum_aoi_loc)+(1-gamma)*((sum_e_MEC-sum_e_loc)/sum_e_loc);
        end 
        p_obj(t)=obj(t);
        if obj(t)>g_obj
            g_obj=obj(t);
            gest=p(:,t);
            each_t=t_total+p(3*tau_num+1:tau_num*4,t);
            each_e=e_total;
            each_gbest_aoi=aoi_MD_now;
            aoi_loc_best=aoi_loc_now;
        end
    end
    data(1,1)=g_obj;
    
    v = zeros(tau_num*4,particlesize);%粒子的移动速度
    v(1:tau_num,t)=w*v(1:tau_num,t)+c1*rand*sign(pest(1:tau_num,t)-p(1:tau_num,t))*0.25+c2*rand*sign(gest(1:tau_num)-p(1:tau_num,t))*0.25;
    v(tau_num+1:tau_num*2,t)=w*v(tau_num+1:tau_num*2,t)+c1*rand*sign(pest(tau_num+1:tau_num*2,t)-p(tau_num+1:tau_num*2,t))*0.6+c2*rand*sign(gest(tau_num+1:tau_num*2)-p(tau_num+1:tau_num*2,t))*0.6;
    v(2*tau_num+1:tau_num*3,t)=w*v(2*tau_num+1:tau_num*3,t)+c1*rand*sign(pest(2*tau_num+1:tau_num*3,t)-p(2*tau_num+1:tau_num*3,t))*0.5+c2*rand*sign(gest(2*tau_num+1:tau_num*3)-p(2*tau_num+1:tau_num*3,t))*0.5;
    v(3*tau_num+1:tau_num*4,t)=w*v(3*tau_num+1:tau_num*4,t)+c1*rand*sign(pest(3*tau_num+1:tau_num*4,t)-p(3*tau_num+1:tau_num*4,t))*0.2+c2*rand*sign(gest(3*tau_num+1:tau_num*4)-p(3*tau_num+1:tau_num*4,t))*0.2;
        
    disp('initial');
    disp(cputime-t0);
end
g = 1;
pest=p;
gest=p(:,randi([1 particlesize]));
t_last_task1 = t_last_task;
g_obj_aoi1=g_obj_aoi;
g_obj_lastt1=g_obj_lastt;
g_obj_energy1=g_obj_energy;
loc_aoi1=loc_aoi;
loc_lastt1=loc_lastt;
edj=e_local_dj;
tlt=t_local_total;
while g <= MG
    w = (0.9 - 0.4) * ((MG - g) / MG) + 0.4;
    c1 = 2 * (sin(pi/2 * (1 - g/MG))^2);
    c2 = 2 * (sin(pi * g / (2 * MG))^2);
    parfor t = 1:particlesize
        pp=pest(:,t);
        local_p = p(:, t);
        local_v = v(:, t);
        
        local_v(1:tau_num) = w * local_v(1:tau_num) + c1 * rand() * sign(pp(1:tau_num) - local_p(1:tau_num))*0.25 + c2 * rand() * sign(gest(1:tau_num) - local_p(1:tau_num))*0.25;
        local_v(tau_num+1:2*tau_num) = w * local_v(tau_num+1:2*tau_num) + c1 * rand() * sign(pp(tau_num+1:2*tau_num) - local_p(tau_num+1:2*tau_num))*0.6 + c2 * rand() * sign(gest(tau_num+1:2*tau_num) - local_p(tau_num+1:2*tau_num))*0.6;
        local_v(2*tau_num+1:3*tau_num) = w * local_v(2*tau_num+1:3*tau_num) + c1 * rand() * sign(pp(2*tau_num+1:3*tau_num) - local_p(2*tau_num+1:3*tau_num))*0.5 + c2 * rand() * sign(gest(2*tau_num+1:3*tau_num) - local_p(2*tau_num+1:3*tau_num))*0.5;
        local_v(3*tau_num+1:4*tau_num) = w * local_v(3*tau_num+1:4*tau_num) + c1 * rand() * sign(pp(3*tau_num+1:4*tau_num) - local_p(3*tau_num+1:4*tau_num))*0.2 + c2 * rand() * sign(gest(3*tau_num+1:4*tau_num) - local_p(3*tau_num+1:4*tau_num))*0.2;
        
        
        local_v(1:tau_num) = max(min(local_v(1:tau_num), vbetamax), vbetamin);
        local_v(tau_num+1:2*tau_num) = max(min(local_v(tau_num+1:2*tau_num), vrmax), vrmin);
        local_v(2*tau_num+1:3*tau_num) = max(min(local_v(2*tau_num+1:3*tau_num), vfmax), vfmin);
        local_v(3*tau_num+1:4*tau_num) = max(min(local_v(3*tau_num+1:4*tau_num), vwtmax), vwtmin);
        
        
        local_p = local_p + local_v;
        
       
        local_p(1:tau_num) = round(max(min(local_p(1:tau_num), betamax), betamin));
        local_p(tau_num+1:2*tau_num) = max(min(local_p(tau_num+1:2*tau_num), rme), rmin);
        local_p(2*tau_num+1:3*tau_num) = max(min(local_p(2*tau_num+1:3*tau_num), fmax), fmin);
        local_p(3*tau_num+1:4*tau_num) = max(min(local_p(3*tau_num+1:4*tau_num), wtmax), wtmin);
        
        
        [obj(t), ~, ~, ~, ~,~, ~, ~, ~,~,~] = compute_fitness(local_p,tau, alltask, dealnum,f_loc,p_trans,T,tbl_num,channelgain,tlt,t_last_task1,g_obj_aoi1,g_obj_lastt1,g_obj_energy1,loc_aoi1,loc_lastt1,owega,edj); % 传入必要参数
        
        
        if obj(t) > p_obj(t)
            pest(:, t) = local_p;
            p_obj(t) = obj(t);
        end
       
        p(:, t) = local_p;
        v(:, t) = local_v;
    end
    
    [current_best, idx] = max(obj);
    if current_best > g_obj
        g_obj = current_best;
        gest = p(:, idx);
        [~, each_t, each_e, each_gbest_aoi, aoi_loc_best,t_last_task,g_obj_aoi,g_obj_lastt,g_obj_energy,loc_aoi,loc_lastt] = compute_fitness(gest, tau, alltask, dealnum,f_loc,p_trans,T,tbl_num,channelgain,t_local_total,t_last_task1,g_obj_aoi1,g_obj_lastt1,g_obj_energy1,loc_aoi1,loc_lastt1,owega,edj);
    end
    data(g+1) = g_obj;
    g = g + 1;
end

if againflag==1
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
end
end

