tic;
clear all; 

% Mass Rb 87: 
m = 1.443161.*10^(-25);
% scattering length Rb 87; 
scat = 98*5.29*10^(-11); % from 98, 95 and 100 (aa, bb, and ab)
% hbar; 
hbar = 1.0546.*10.^(-34);

% From Campbell Expt Slides #11: 
%http://online.kitp.ucsb.edu/online/boptilatt10/campbell/pdf/Campbell_OpticalLattices_KITP.pdf
% Number of atoms: 
%Natom = 1.5*10^5; % 1.5*10^5; % using 1/2 that in gretchens expt so density the same
% Trap Frequencies: 
%omegaz = 500; % Hz
%omegar = 100; % Hz
%r0 = 20*10^(-6);

% Initial Grids setting up grids
xmax=0.9*10^(-5); Ngx=2^7; [x,dx,px,dpx]=fftdef(xmax,Ngx); 
zmax=0.55*10^(-4); Ngz=2^9; [z,dz,pz,dpz]=fftdef(zmax,Ngz); 
[pzm,pxm]=meshgrid(pz,px);
[zm,xm]=meshgrid(z,x);

%-------------------------------------------------------------------------  
% Constructing the Potential
% Harmonic Trap
%% 
omegar = 1500; 
a = 0.4e-5;%2*xmax/11;%0.5e-5;


n=21;
r= zeros(Ngz, Ngx, n);
Vtrap = zeros(Ngz, Ngx, n);

offset = floor(n/2)+1;% if odd addd one, else?
for ii=1:n
    r(:, :, ii) = sqrt((zm+(ii-offset)*a).^2 + xm.^2);
    Vtrap(:, :, ii) = 0.5*m*omegar*omegar.*r(ii).^2;
end

%V = sum(Vtrap);
%Vmin = min(Vtrap);

%%


% r_10=sqrt((zm-10*a).^2 + xm.^2); 
% Vtrap_10 = 0.5*m*omegar*omegar.*r_10.^2;
% 
% r_9=sqrt((zm-9*a).^2 + xm.^2); 
% Vtrap_9 = 0.5*m*omegar*omegar.*r_9.^2;
% 
% r_8=sqrt((zm-8*a).^2 + xm.^2); 
% Vtrap_8 = 0.5*m*omegar*omegar.*r_8.^2;
% 
% r_7=sqrt((zm-7*a).^2 + xm.^2); 
% Vtrap_7 = 0.5*m*omegar*omegar.*r_7.^2;
% 
% r_6=sqrt((zm-6*a).^2 + xm.^2); 
% Vtrap_6 = 0.5*m*omegar*omegar.*r_6.^2;
% 
% r_5=sqrt((zm-5*a).^2 + xm.^2); 
% Vtrap_5 = 0.5*m*omegar*omegar.*r_5.^2;
% 
% r_4=sqrt((zm-4*a).^2 + xm.^2); 
% Vtrap_4 = 0.5*m*omegar*omegar.*r_4.^2;
% 
% r_3=sqrt((zm-3*a).^2 + xm.^2); 
% Vtrap_3 = 0.5*m*omegar*omegar.*r_3.^2;

% r_2=sqrt((zm-2*a).^2 + xm.^2); 
% Vtrap_2 = 0.5*m*omegar*omegar.*r_2.^2;
% 
% r_1=sqrt((zm-a).^2 + xm.^2); 
% Vtrap_1 = 0.5*m*omegar*omegar.*r_1.^2;
% 
% r_0=sqrt((zm).^2 + xm.^2); 
% Vtrap_0 = 0.5*m*omegar*omegar.*r_0.^2;
% 
% r1=sqrt((zm+a).^2 + xm.^2); 
% Vtrap1 = 0.5*m*omegar*omegar.*r1.^2;
% 
% r2=sqrt((zm+2*a).^2 + xm.^2); 
% Vtrap2 = 0.5*m*omegar*omegar.*r2.^2;

% r3=sqrt((zm+3*a).^2 + xm.^2); 
% Vtrap3 = 0.5*m*omegar*omegar.*r3.^2;
% 
% r4=sqrt((zm+4*a).^2 + xm.^2); 
% Vtrap4 = 0.5*m*omegar*omegar.*r4.^2;
% 
% r5=sqrt((zm+5*a).^2 + xm.^2); 
% Vtrap5 = 0.5*m*omegar*omegar.*r5.^2;
% 
% r6=sqrt((zm+6*a).^2 + xm.^2); 
% Vtrap6 = 0.5*m*omegar*omegar.*r6.^2;
% 
% % r7=sqrt((zm+7*a).^2 + xm.^2); 
% Vtrap7 = 0.5*m*omegar*omegar.*r7.^2;
% 
% r8=sqrt((zm+8*a).^2 + xm.^2); 
% Vtrap8 = 0.5*m*omegar*omegar.*r8.^2;
% 
% r9=sqrt((zm+9*a).^2 + xm.^2); 
% Vtrap9 = 0.5*m*omegar*omegar.*r9.^2;
% 
% r10=sqrt((zm+10*a).^2 + xm.^2); 
% Vtrap10 = 0.5*m*omegar*omegar.*r10.^2;

ii = 20; 
Vfindmin = Vtrap(:,ii,21)+Vtrap(:,ii,20);
[Vminvalue,mid22]=min(Vfindmin);
mid2=floor(mid22);

Vfindmin = Vtrap(:,:,9)+Vtrap(:,:,8);
[Vminvalue,mid33]=min(Vfindmin);
mid3=floor(mid33);

mid4=mid3-mid2;


stop;
% mid2=floor(Ngx/(2^21));

V11 = min(Vtrap1,Vtrap_1);
V22 = min(Vtrap2,Vtrap_2);
V33 = min(Vtrap3,Vtrap_3);
V44 = min(Vtrap4,Vtrap_4);
V55 = min(Vtrap5,Vtrap_5);
V66 = min(Vtrap6,Vtrap_6);
V77 = min(Vtrap7,Vtrap_7);
V88 = min(Vtrap8,Vtrap_8);
V99 = min(Vtrap9,Vtrap_9);
V1010 = min(Vtrap10,Vtrap_10);

V01 = min(Vtrap_0,V11);
V23 = min(V22,V33);
V45 = min(V44,V55);
V67 = min(V66,V77);
V89 = min(V88,V99);

V03 = min(V01,V23);
V47 = min(V45,V67);
V810 = min(V89,V1010);

V07 = min(V03,V47);


V1 = min(V07,V810);
 
% So you can see what it looks like: 

    figure(1);
    subplot(3,2,1);
    plot(V1(mid2,:),'-*');  
 
    subplot(3,2,2); 
    pcolor(V1); shading interp;  
    
    subplot(3,2,[3,4]);
    mesh(V1);  
 
    subplot(3,2,[5,6]);
    hold; 
    plot(Vtrap_2(ii,:),'b-*')
    plot(Vtrap_1(ii,:),'r-*')
    plot(Vtrap_0(ii,:),'y-*')
    plot(Vtrap1(ii,:),'g-*')
    plot(Vtrap2(ii,:),'k-*')
    hold;
    drawnow;
    
%         figure(20);
%     hold; 
%     plot(Vtrap_2(ii,:),'b-*')
%     plot(Vtrap_1(ii,:),'r-*')
%     plot(Vtrap_0(ii,:),'y-*')
%     plot(Vtrap1(ii,:),'g-*')
%     plot(Vtrap2(ii,:),'k-*')
%     hold;  
% 
%     figure(18); 
%     pcolor(V1); shading interp;  
%     figure(10); mesh(V1);  
%     figure(13);   
%     plot(V1(mid2,:),'-*');
    
%-------------------------------------------------------------------------
% Initial wavefunction
alpha = sqrt((m*omegar)/hbar); 
theta = atan2(xm,zm);
% initial state
wf1 = (alpha/sqrt(pi)).*exp(-0.5*(alpha.^2).*r_0.^2);
% First excited state 
% wf1 = alpha.*r_0.*exp(1i*theta).*wf1;
rho= wf1.*conj(wf1);
  
   
   figure(2);clf; 
   pcolor(rho(:,:)); shading interp;
   
%-------------------------------------------------------------------------
% Number of steps and stepsize

steps=1000000; 
Dt=0.1*10.^(-4);

% Number of samples of wavefunction over the run

Samples = 3000; 
SampleInc = floor(steps/Samples);

% -------------------------------------------------------------------------
% Initialise stuff so we are not changing size every timestep
% Things we are outputting for a number of steps

Eprop1=zeros(1,Samples);
Natom1=zeros(1,Samples);
t_1 = zeros(1,Samples);

 
% Na_Stuff = zeros(1,Samples);

Na0 = zeros(1,Samples);
Na1 = zeros(1,Samples);
Na_1 = zeros(1,Samples);
Na2 = zeros(1,Samples);
Na_2 = zeros(1,Samples);
Na3 = zeros(1,Samples);
Na_3 = zeros(1,Samples);
Na4 = zeros(1,Samples);
Na_4 = zeros(1,Samples);
Na5 = zeros(1,Samples);
Na_5 = zeros(1,Samples);
Na6 = zeros(1,Samples);
Na_6 = zeros(1,Samples);
Na7 = zeros(1,Samples);
Na_7 = zeros(1,Samples);
Na8 = zeros(1,Samples);
Na_8 = zeros(1,Samples);
Na9 = zeros(1,Samples);
Na_9 = zeros(1,Samples);
Na10 = zeros(1,Samples);
Na_10 = zeros(1,Samples);


% -------------------------------------------------------------------------
% Kinetic Energy 
% +[(hbar^2)/(2m)]*del^2

K =((hbar).^2/(2.*m))*(pxm.^2 + pzm.^2);
EKp=exp(-1i.*K.*Dt./hbar);
EVr1= exp(-1i*V1.*Dt./hbar);
bb =0; 

% shoddy incrementation for stuff to output
for i2=1:steps;
    
    %--------------------------------------------------
    % step 1
    
    
    wf1=EVr1.*wf1;
    
    % step 2 
    
    wf1=ifft2(EKp.*fft2(wf1));  
    
    % step 3
    
    %EVr1= exp(-0.5*1i*V1.*Dt./hbar);
    %wf1=EVr1.*wf1;
    
    %--------------------------------------------------
    % Sampling things
    % Do this as often as you like :)
    if mod(i2,SampleInc)==0 
       bb = bb+1; 
        
       % Energy
       Epot = V1.*wf1; 
       Ekin=ifft2(K.*fft2(wf1));
       Eprop1(bb)=real(sum(sum(conj(wf1).*(Ekin+Epot))))*dx*dz;

       % Number of atoms - should be 1! 
       Natom1(bb)=(sum(sum(wf1.*conj(wf1)))*dx*dz);
%      N1=sum(sum(wf1.*conj(wf1()))*dx*dz;
       
       
       
       Na10(bb)=sum(sum(wf1(:,1:mid2).*conj(wf1(:,1:mid2)))*dx*dz);
       
       Na9(bb)=sum(sum(wf1(:,mid2:mid3).*conj(wf1(:,mid2:mid3)))*dx*dz);
       Na8(bb)=sum(sum(wf1(:,mid3:mid4+mid3).*conj(wf1(:,mid3:mid4+mid3)))*dx*dz);
       Na7(bb)=sum(sum(wf1(:,mid4+mid3:2*mid4+mid3).*conj(wf1(:,mid4+mid3:2*mid4+mid3)))*dx*dz);
       Na6(bb)=sum(sum(wf1(:,2*mid4+mid3:3*mid4+mid3).*conj(wf1(:,2*mid4+mid3:3*mid4+mid3)))*dx*dz);
       Na5(bb)=sum(sum(wf1(:,3*mid4+mid3:4*mid4+mid3).*conj(wf1(:,3*mid4+mid3:4*mid4+mid3)))*dx*dz);
       Na4(bb)=sum(sum(wf1(:,4*mid4+mid3:5*mid4+mid3).*conj(wf1(:,4*mid4+mid3:5*mid4+mid3)))*dx*dz); 
       Na3(bb)=sum(sum(wf1(:,5*mid4+mid3:6*mid4+mid3).*conj(wf1(:,5*mid4+mid3:6*mid4+mid3)))*dx*dz);
       Na2(bb)=sum(sum(wf1(:,6*mid4+mid3:7*mid4+mid3).*conj(wf1(:,6*mid4+mid3:7*mid4+mid3)))*dx*dz);
       Na1(bb)=sum(sum(wf1(:,7*mid4+mid3:8*mid4+mid3).*conj(wf1(:,7*mid4+mid3:8*mid4+mid3)))*dx*dz);
        
       Na0(bb)=sum(sum(wf1(:,8*mid4+mid3:9*mid4+mid3).*conj(wf1(:,8*mid4+mid3:9*mid4+mid3)))*dx*dz);
    
       Na_1(bb)=sum(sum(wf1(:,9*mid4+mid3:10*mid4+mid3).*conj(wf1(:,9*mid4+mid3:10*mid4+mid3)))*dx*dz);
       Na_2(bb)=sum(sum(wf1(:,10*mid4+mid3:11*mid4+mid3).*conj(wf1(:,10*mid4+mid3:11*mid4+mid3)))*dx*dz);
       Na_3(bb)=sum(sum(wf1(:,11*mid4+mid3:12*mid4+mid3).*conj(wf1(:,11*mid4+mid3:12*mid4+mid3)))*dx*dz);
       Na_4(bb)=sum(sum(wf1(:,12*mid4+mid3:13*mid4+mid3).*conj(wf1(:,12*mid4+mid3:13*mid4+mid3)))*dx*dz);
       Na_5(bb)=sum(sum(wf1(:,13*mid4+mid3:14*mid4+mid3).*conj(wf1(:,13*mid4+mid3:14*mid4+mid3)))*dx*dz);
       Na_6(bb)=sum(sum(wf1(:,14*mid4+mid3:15*mid4+mid3).*conj(wf1(:,14*mid4+mid3:15*mid4+mid3)))*dx*dz);
       Na_7(bb)=sum(sum(wf1(:,15*mid4+mid3:16*mid4+mid3).*conj(wf1(:,15*mid4+mid3:16*mid4+mid3)))*dx*dz);
       Na_8(bb)=sum(sum(wf1(:,16*mid4+mid3:17*mid4+mid3).*conj(wf1(:,16*mid4+mid3:17*mid4+mid3)))*dx*dz);
       Na_9(bb)=sum(sum(wf1(:,17*mid4+mid3:18*mid4+mid3).*conj(wf1(:,17*mid4+mid3:18*mid4+mid3)))*dx*dz);
       Na_10(bb)=sum(sum(wf1(:,18*mid4+mid3:Ngz).*conj(wf1(:,18*mid4+mid3:Ngz)))*dx*dz);

       
       

       % current time
       t_1(bb) = Dt*(i2); 
 
    end
    if rem(i2,100)==0
        figure(3)
        subplot(2,1,1);
        pcolor(wf1(:,:).*conj(wf1(:,:)));shading interp; axis equal; colormap bone; axis equal;  %shading interp;colormap bone; axis equal;
        axis([0 Ngz 0 Ngx]);
        title('Density')
        figure(3);
        subplot(2,1,2);
        pcolor(angle(wf1(:,:))); shading interp; axis equal; colormap bone; axis equal; 
        axis([0 Ngz 0 Ngx]);
        title('Phase')
        drawnow; 
        
        
%         % Plotting energy and atom number
%         figure(4);clf; 
%         subplot(1,2,1)
%         plot(t_1(1:bb),Eprop1(1:bb),'b -*')
%         ylabel('Energy')
%         xlabel('t')
%         subplot(1,2,2)
%         plot(t_1(1:bb),Natom1(1:bb),'b -*')
%         ylabel('N')
%         xlabel('t')
%         drawnow; 
        
        

        
%                 hold;
%                 figure(6);clf;
%                 hold;
%                 plot(Na10(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_9(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_8(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_7(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_6(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_5(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_4(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_3(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_2(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na_1(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na0(1:bb),t_1(1:bb),'r -*');shading interp;
%                 plot(Na1(1:bb),t_1(1:bb),'g -*');shading interp;  
%                 plot(Na2(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na3(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na4(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na5(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na6(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na7(1:bb),t_1(1:bb),'g -*');shading interp;
%                 plot(Na8(1:bb),t_1(1:bb),'g -*');shading interp;    
%                 plot(Na9(1:bb),t_1(1:bb),'g -*');shading interp;    
%                 plot(Na10(1:bb),t_1(1:bb),'g -*');shading interp;    
% 
%                 hold;
%                 drawnow;
                %;shading interp;%axis equal; colormap bone; axis equal;  %shading interp;colormap bone; axis equal;
%                 plot(t_1(1:bb),Na9(1:bb),'b -*')%
%                 plot(t_1(1:bb),Na8(1:bb),'b -*')%;shading interp;
%                hold;
%             
         
        
%---------------------------------------------------------------------------        
%         figure(5);clf;
        %title('Density @ well')
% 
%         for jj=1:21
%             if jj==1
%                 subplot(3,7,1)
%                 pcolor(wf1(:,1:(jj)*mid2).*conj(wf1(:,1:(jj)*mid2)));shading interp; %axis equal; colormap bone; axis equal;  %shading interp;colormap bone; axis equal;
%                 ylabel('Density')
%                 xlabel('t')
%                 
%             else
%                 subplot(3,7,jj)
%                 pcolor(wf1(:,(jj)*mid2:(jj+1)*mid2).*conj(wf1(:,(jj)*mid2:(jj+1)*mid2)));shading interp; %axis equal; colormap bone; axis equal;  %shading interp;colormap bone; axis equal;
%                 ylabel('Density')
%                 xlabel('t')
%             end
%             
%         end
%               drawnow; 
%  
         figure(6);clf;
%         title('Density')

%         for kk=1:20
%             if kk==1
%                 subplot(3,7,1)
% %               plot(wf1(:,1:(jj)*mid2).*conj(wf1(:,1:(jj)*mid2)));
% %                  plot(wf1(1:mid2).*conj(wf1(1:mid2)),'b -*');
%               plot(t_1(1:(bb)),wf1(:,1:(jj)*bb*mid2).*conj(wf1(1:(bb)*jj*mid2)),'b -*')
%                 ylabel('Density')
%                 xlabel('t')
%                 
%             else
%                 subplot(3,7,kk)
%                 plot(wf1((kk-1)*mid2:(kk)*mid2).*conj(wf1((kk-1)*mid2,(kk)*mid2)),'b -*');
% %               plot(wf1(:,(jj)*mid2:(jj+1)*mid2).*conj(wf1(:,(jj)*mid2:(jj+1)*mid2)));
%                 ylabel('Density')
%                 xlabel('t')
%             end
%             
%             
%             
%         end
%               drawnow;    
        
              
                subplot(7,3,1)
                plot(t_1(1:bb),Na10(1:bb),'b -*');shading interp;%axis equal; colormap bone; axis equal;  %shading interp;colormap bone; axis equal;
                ylabel('Density @ well 10')
                xlabel('t')
               
                subplot(7,3,2)
                plot(t_1(1:bb),Na9(1:bb),'b -*');shading interp;
                ylabel('Well 9')
                xlabel('t')
           
                subplot(7,3,3)
                plot(t_1(1:bb),Na8(1:bb),'b -*');shading interp;
                ylabel(' Well 8')
                xlabel('t')
                
                subplot(7,3,4)
                plot(t_1(1:bb),Na7(1:bb),'b -*');shading interp;
                ylabel('Well 7')
                xlabel('t')
                 
                subplot(7,3,5)
                plot(t_1(1:bb),Na6(1:bb),'b -*');shading interp;
                ylabel('Well 6')
                xlabel('t')
                
                subplot(7,3,6)
                plot(t_1(1:bb),Na5(1:bb),'b -*');shading interp;
                ylabel('Well 5')
                xlabel('t')
                
                subplot(7,3,7)
                plot(t_1(1:bb),Na4(1:bb),'b -*');shading interp;
                ylabel('Well 4')
                xlabel('t')
                
                subplot(7,3,8)
                plot(t_1(1:bb),Na3(1:bb),'b -*');shading interp;
                ylabel('Well 3')
                xlabel('t')
                
                subplot(7,3,9)
                plot(t_1(1:bb),Na2(1:bb),'b -*');shading interp;
                ylabel('Well 2')
                xlabel('t')
                
                subplot(7,3,10)
                plot(t_1(1:bb),Na1(1:bb),'b -*');shading interp;
                ylabel('Well 1')
                xlabel('t')
                
                subplot(7,3,11)
                plot(t_1(1:bb),Na0(1:bb),'r -*');shading interp;
                ylabel(' Middle trap')
                xlabel('t')
                
                subplot(7,3,12)
                plot(t_1(1:bb),Na_1(1:bb),'b -*');shading interp;
                ylabel('Well 1R')
                xlabel('t')
                
                subplot(7,3,13)
                plot(t_1(1:bb),Na_2(1:bb),'b -*');shading interp;
                ylabel('Well 2R')
                xlabel('t')
                
                subplot(7,3,14)
                plot(t_1(1:bb),Na_3(1:bb),'b -*');shading interp;
                ylabel('Well 3R')
                xlabel('t')
                
                subplot(7,3,15)
                plot(t_1(1:bb),Na_4(1:bb),'b -*');shading interp;
                ylabel('Well 4R')
                xlabel('t')
                
                subplot(7,3,16)
                plot(t_1(1:bb),Na_5(1:bb),'b -*');shading interp;
                ylabel('Well 5R')
                xlabel('t')
                
                subplot(7,3,17)
                plot(t_1(1:bb),Na_6(1:bb),'b -*');shading interp;
                ylabel('Well 6R')
                xlabel('t')
                
                subplot(7,3,18)
                plot(t_1(1:bb),Na_7(1:bb),'b -*');shading interp;
                ylabel('Well 7R')
                xlabel('t')
                
                subplot(7,3,19)
                plot(t_1(1:bb),Na_8(1:bb),'b -*');shading interp;
                ylabel('Well 8R')
                xlabel('t')
                
                subplot(7,3,20)
                plot(t_1(1:bb),Na_9(1:bb),'b -*');shading interp;
                ylabel('Well 9R')
                xlabel('t')
                
                subplot(7,3,21)
                plot(t_1(1:bb),Na_10(1:bb),'b -*');shading interp;
                ylabel('Well 10R')
                xlabel('t')
                
                drawnow;
              
              
    end         
end 
 

        
save('schroedinger.mat','wf1')

toc
