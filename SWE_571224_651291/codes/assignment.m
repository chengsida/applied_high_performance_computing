function waterequation()
   clear all; close all; clc;
   global g N_x N_y  delta_x delta_y;
   % simulation parameters
   x_min=0.00;
   x_max=100.00;
   y_min=0.00;
   y_max=100.00;
   t_min=0.00;
   t_max=100.00;
   g=9.81;
   delta_x=1;
   delta_y=1;
   delta_t=0.1;
   x=x_min:delta_x:x_max;
   y=y_min:delta_y:y_max;
   t=t_min:delta_t:t_max;
   N_x=length(x);
   N_y=length(y);
   N_t=length(t);
 
   
  %allocate arrays
   h=zeros(N_x,N_y);
   vx=zeros(N_x,N_y);
   vy=zeros(N_x,N_y);

   
   %set initial condition
 
   for i=1:N_x
       for j=1:N_y
           h (i,j)=1+0.5*exp((-1/25)*(((x(i)-30).^2)+((y(j)-30).^2)));
       end
   end
   
   vx(:,:)=0;
   vy(:,:)=0;
   

 
 figure('WindowStyle', 'docked');
    Solution        = surf(x, y, h(:,:));
    axis('equal', [x_min x_max y_min y_max 0 5]);
    axis off;
    grid off;
    
    view([45 25]);
    color('winter');
    drawnow;
 
   
    for l=1:N_t-1

        [kx1,ky1,kh1]          = f(vx(:,:),vy(:,:),h(:,:));
        [kx2,ky2,kh2]          = f(vx(:,:) + (delta_t/2)*kx1(:,:),vy(:,:)+(delta_t/2)*ky1(:,:),h(:,:)+(delta_t/2)*kh1(:,:));
        [kx3,ky3,kh3]          = f(vx(:,:) + (delta_t/2)*kx2(:,:),vy(:,:)+(delta_t/2)*ky2(:,:),h(:,:)+(delta_t/2)*kh2(:,:));
        [kx4,ky4,kh4]          = f(vx(:,:) + delta_t*kx3(:,:),vy(:,:)+delta_t*ky3(:,:),h(:,:)+delta_t*kh3(:,:));
        vx(:,:)  = vx(:,:) + delta_t  *(kx1(:,:)/6 + kx2(:,:)/3 + kx3(:,:)/3 + kx4(:,:)/6);
        vy(:,:)  = vy(:,:) + delta_t  *(ky1(:,:)/6 + ky2(:,:)/3 + ky3(:,:)/3 + ky4(:,:)/6);
        h(:,:)   = h(:,:) + delta_t  *(kh1(:,:)/6 + kh2(:,:)/3 + kh3(:,:)/3 + kh4(:,:)/6);
        
        
        % Plot the solution
      set(Solution, 'ZData', h(:,:));
      title(['t = ' num2str(t(l+1))]);
      drawnow;
        
    end
    
    

end


function [kx,ky,kh]=f(vx,vy,h)
    
    global  delta_x delta_y N_x N_y g;
    kx=zeros(N_x,N_y);
    ky=zeros(N_x,N_y);
    kh=zeros(N_x,N_y);
    
    for i=1:N_x;
            x_a=i-2;
            x_b=i-1;
            x_c=i+1;
            x_d=i+2;
            if x_a==-1
                x_a=N_x-1;
            else if x_a==0
                x_a=N_x;
                end
            end
            if x_b==0
                x_b=N_x;
            end
            if x_c==N_x+1;
                x_c=1;
            end
            if x_d==N_x+2;
                x_d=2;
            else if x_d==N_x+1;
                    x_d=1;
                end
            end
        for j=1:N_y;
            
            y_e=j-2;
            y_f=j-1;
            y_g=j+1;
            y_h=j+2;
            
            if y_e==-1
                y_e=N_y-1;
            else if y_e==0
                y_e=N_y;
                end
            end
            if y_f==0
                y_f=N_y;
            end
            if y_g==N_y+1;
                y_g=1;
            end
            if y_h==N_y+2;
                y_h=2;
            else if y_h==N_y+1;
                    y_h=1;
                end
            end
            kx(i,j)=(((-1)*vx(i,j))/(12*delta_x))*(vx(x_a,j)-8*vx(x_b,j)+8*vx(x_c,j)-vx(x_d,j))+(((-1)*vy(i,j))/(12*delta_y))*(vx(i,y_e)-8*vx(i,y_f)+8*vx(i,y_g)-vx(i,y_h))+((-1)*g/(12*delta_x))*(h(x_a,j)-8*h(x_b,j)+8*h(x_c,j)-h(x_d,j));
            ky(i,j)=(((-1)*vx(i,j))/(12*delta_x))*(vy(x_a,j)-8*vy(x_b,j)+8*vy(x_c,j)-vy(x_d,j))+(((-1)*vy(i,j))/(12*delta_y))*(vy(i,y_e)-8*vy(i,y_f)+8*vy(i,y_g)-vy(i,y_h))+((-1)*g/(12*delta_y))*(h(i,y_e)-8*h(i,y_f)+8*h(i,y_g)-h(i,y_h));
            kh(i,j)=(-1/(12*delta_x))*(vx(x_a,j)*h(x_a,j)-8*vx(x_b,j)*h(x_b,j)+8*vx(x_c,j)*h(x_c,j)-vx(x_d,j)*h(x_d,j))+(-1/(12*delta_y))*(vy(i,y_e)*h(i,y_e)-8*vy(i,y_f)*h(i,y_f)+8*vy(i,y_g)*h(i,y_g)-vy(i,y_h)*h(i,y_h));
        end
    end
end

    
        
    
  
   