% % 初始化
% clear;
% clc;
function[gest,each_t,each_e,each_gbest_aoi,each_obj,g_obj,e_local_dj,aoi_loc_best,T,dealnum]=pso(dealnum,alltask,tau,f_avail,B_avail,againflag,aa,e_local_dj)
global t_last_task;
global g_obj_aoi;
global g_obj_lastt;
global g_obj_energy;
global loc_aoi;
global loc_lastt;
global p;
global t_local_total;
timeout=0;
sourceout=0;               
MG =1800;  
data=zeros(MG+1,1);
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

%读取任务相关参数
%global s;
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

%计算信道增益
distance=(5+5*rand(tau_num,1));
channelgain=distance.^(-4);
backin=10^(-13);

if againflag==0
    t_local_total=all_d_bit.*all_owega./(all_f_loc*10^9);
    e_loc_total=5*10^(-27).*(all_f_loc*10^9).^3.*t_local_total;
end
t_loc_ddl=owega.*d_bit./(f_loc*10^9);
T=t_loc_ddl.*(randi([15,20])/10);

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
    
    obj=zeros(1,particlesize);
    p_obj=obj;
    pest=p;
    gest=p(:,randi([1 particlesize]));
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
            timeout=timeout+1;
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
    
    v = zeros(tau_num*4,particlesize);
   % v(1:tau_num,:) = rand(tau_num,particlesize)*(betamax-betamin)-p(1:tau_num,:);
   % v(tau_num+1:tau_num*2,:) = rand(tau_num,particlesize)*(rme-rmin)-p(tau_num+1:tau_num*2,:);
    %v(2*tau_num+1:tau_num*3,:)=rand(tau_num,particlesize)*(fmax-fmin)-p(2*tau_num+1:tau_num*3,:);
   % v(3*tau_num+1:tau_num*4,:)=rand(tau_num,particlesize)*(wtmax-wtmin)-p(3*tau_num+1:tau_num*4,:);
    v(1:tau_num,:) = rand(tau_num,particlesize)*(betamax-betamin);
    v(tau_num+1:tau_num*2,:) = rand(tau_num,particlesize)*(rme-rmin);
    v(2*tau_num+1:tau_num*3,:)=rand(tau_num,particlesize)*(fmax-fmin);
    v(3*tau_num+1:tau_num*4,:)=rand(tau_num,particlesize)*(wtmax-wtmin);
    %disp(v);
    disp('initial');
    disp(cputime-t0);
    g=1;
    wmax=0.9;
    wmin=0.4;
    t0=cputime;
    while g<=MG
        w=(wmax-wmin)*((MG-g)/MG)+wmin;
        c1=2*(sin(pi/2*(1-g/MG)))^2;
        c2=2*(sin(pi*g/(2*MG)))^2;
        for t =1:particlesize
            v(1:tau_num,t)=w*v(1:tau_num,t)+c1*rand*sign(pest(1:tau_num,t)-p(1:tau_num,t))*0.25+c2*rand*sign(gest(1:tau_num)-p(1:tau_num,t))*0.25;
            v(tau_num+1:tau_num*2,t)=w*v(tau_num+1:tau_num*2,t)+c1*rand*sign(pest(tau_num+1:tau_num*2,t)-p(tau_num+1:tau_num*2,t))*0.6+c2*rand*sign(gest(tau_num+1:tau_num*2)-p(tau_num+1:tau_num*2,t))*0.6;
            v(2*tau_num+1:tau_num*3,t)=w*v(2*tau_num+1:tau_num*3,t)+c1*rand*sign(pest(2*tau_num+1:tau_num*3,t)-p(2*tau_num+1:tau_num*3,t))*0.5+c2*rand*sign(gest(2*tau_num+1:tau_num*3)-p(2*tau_num+1:tau_num*3,t))*0.5;
            v(3*tau_num+1:tau_num*4,t)=w*v(3*tau_num+1:tau_num*4,t)+c1*rand*sign(pest(3*tau_num+1:tau_num*4,t)-p(3*tau_num+1:tau_num*4,t))*0.2+c2*rand*sign(gest(3*tau_num+1:tau_num*4)-p(3*tau_num+1:tau_num*4,t))*0.2;    
            for j = 1:tau_num
                s_jt=1/(1+exp(-v(j,t)));
                if s_jt>=rand
                    p(j,t)=1;
                else
                    p(j,t)=0;
                end
            end
            
            for j=(tau_num+1):tau_num*2
                if v(j,t)>vrmax
                    v(j,t)=vrmax;
                elseif v(j,t)<vrmin
                    v(j,t)=vrmin;
                end
            end
            for j=(2*tau_num+1):tau_num*3
                if v(j,t)>vfmax
                    v(j,t)=vfmax;
                elseif v(j,t)<vfmin
                    v(j,t)=vfmin;
                end
            end
            for j=(3*tau_num+1):tau_num*4
                if v(j,t)>vwtmax
                    v(j,t)=vwtmax;
                elseif v(j,t)<vwtmin
                    v(j,t)=vwtmin;
                end
            end
            p(tau_num+1:tau_num*4,t)=p(tau_num+1:tau_num*4,t)+v(tau_num+1:tau_num*4,t);
            
            r_total=0;
            for j=(tau_num+1):tau_num*2
                if p(j,t)>rme
                    p(j,t)=rme;
                elseif p(j,t)<rmin
                    p(j,t)=rmin;
                end
                r_total=r_total+p(j,t)*p(j-tau_num,t);
            end
            rate=zeros(tau_num,particlesize);
            pc=p_trans.*channelgain;
            for i=1:tau_num
                for u=1:particlesize
                    rate(i,u)= p(tau_num+i,u)*1024*1024*8*log2(1+pc(i)/(backin+sum(pc.*p(1:tau_num,u))-pc(i)));
                end
            end
            
            f_total=0;
            for j=(2*tau_num+1):tau_num*3
                if p(j,t)>fmax
                    p(j,t)=fmax;
                elseif p(j,t)<fmin
                   p(j,t)=fmin;
                end
                f_total=f_total+p(j,t)*p(j-2*tau_num,t);
            end
            
            for j=(3*tau_num+1):tau_num*4
                if p(j,t)>wtmax
                    p(j,t)=wtmax;
                elseif p(j,t)<wtmin
                    p(j,t)=wtmin;
                end
            end
            %         disp(v(:,t));
            %         disp(p(:,t));
            %         disp(f_total);
            if f_total>f_avail+10^(-10)|| r_total>B_avail
                sourceout=sourceout+1;
                obj(t)=-1;
            else
                aoi_loc_now=zeros(1,50);
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
                e_exe=5*10^(-27).*(f_loc).^3.*t_exe_loc;
                t_trans=p(1:tau_num,t).*d_bit./rate(:,t);
                e_trans=p_trans.*t_trans;
                t_exe_ES=p(1:tau_num,t).*owega.*d_bit./(p(2*tau_num+1:tau_num*3,t)*10^9);
                t_total=t_exe_loc+(t_trans+t_exe_ES);
                e_total=p(1:tau_num,t).*e_trans+(1-p(1:tau_num,t)).*e_exe;
                
                
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
                            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+slsj)/(g_obj_lastt(mdnum)+p(3*tau_num+i,t)+t_total(i));
                            aoi_MD_now(mdnum)=aoi_new;
                            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
                        elseif dealnum(mdnum)==1
                            ssj=(t_total(i)+p(3*tau_num+i,t))*(t_total(i)+p(3*tau_num+i,t))/2;
                            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+ssj)/(g_obj_lastt(mdnum)+p(3*tau_num+i,t)+t_total(i));
                            aoi_MD_now(mdnum)=aoi_new;
                            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
                        else
                            ssj=(t_total(i)+p(3*tau_num+i,t))*(t_total(i)+p(3*tau_num+i,t))/2;    
                            spx=t_last_task(mdnum)*t_total(i);
                            aoi_new=(g_obj_aoi(mdnum)*g_obj_lastt(mdnum)+spx+ssj)/(g_obj_lastt(mdnum)+p(3*tau_num+i,t)+t_total(i));
                            aoi_MD_now(mdnum)=aoi_new;
                            e_MEC(mdnum)=g_obj_energy(mdnum)+e_total(i);
                        end
                    end      
                    sum_aoi_loc=0;
                    sum_e_loc=0;sum_aoi_MEC=0;
                    sum_e_MEC=0;
                    for i=1:tau_num
                        sum_aoi_loc=sum_aoi_loc+aoi_loc_now(tbl_num(i));
                        sum_e_loc=sum_e_loc+e_local_dj(tbl_num(i),dealnum(tbl_num(i)));
                        sum_aoi_MEC=sum_aoi_MEC+aoi_MD_now(tbl_num(i));
                        sum_e_MEC=sum_e_MEC+e_MEC(tbl_num(i));
                    end
                    obj(t)=gamma*((sum_aoi_loc-sum_aoi_MEC)/sum_aoi_loc)+(1-gamma)*((sum_e_loc-sum_e_MEC)/sum_e_loc);
                end
            end
            if obj(t)>p_obj(t)
                pest(:,t)=p(:,t);
                p_obj(t)=obj(t);
            end
            if obj(t)>g_obj
                g_obj=obj(t);
                gest=p(:,t);
                each_t=t_total+p(3*tau_num+1:tau_num*4,t);
                each_e=e_total;
                each_gbest_aoi=aoi_MD_now;
                aoi_loc_best=aoi_loc_now;
            end
        end
        g=g+1;
        data(g,1)=g_obj;
    end
    disp('iter');
    disp(cputime-t0);
    k=k+1;
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
