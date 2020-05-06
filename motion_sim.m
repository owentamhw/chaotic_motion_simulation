% 2018/11/21
% PHYS4150 Project
% CHAOTIC BEHAVIOR IN A PLANE PENDULUM WITH A SINUSOIDAL DRIVING TORQUE AT THE PIVOT POINT

function motion_sim
    close all
    clear all
    clc
    
    num_methods = 1;    %need to match the number of methods used

    %====================================================================
    %Array, variable initialization  
    [g,l,theta_0,w_0,~,w_ext,~,h,data_pts,duration] = get_var();
 
    if data_pts > 500000
        disp('Datapoints exceeds 200000, simulation terminated.');
        return
    end
    [time, theta, w] = deal(zeros(data_pts, num_methods));
    theta(1,1:num_methods) = theta_0;
    w(1,1:num_methods) = w_0;
    
    %====================================================================
    %RK4 Method
    [time(:,1),theta(:,1),w(:,1)] = RK4(time(:,1),theta(:,1),w(:,1));

    %====================================================================
    %Euler Method
    %[time(:,2),theta(:,2),w(:,2)] = euler_1(time(:,2),theta(:,2),w(:,2));
    
    %====================================================================
    %Error with ode45
    [t,y] = builtin_ode45();
    abs_err = theta(end,1)-y(end,1);
    rel_err =(theta(end,1)-y(end,1))/min(theta(end,1),y(end,1));
    fprintf('Absolue & Relative Error With ode45: %.12f %.12f\n',abs_err,rel_err);
    
    %====================================================================
    %Plot Graph
    plot_graph(time(:,1), theta(:,1), w(:,1)); %RK4
    %plot_graph(time(:,2), theta(:,2), w(:,2)); %Euler
    
     %====================================================================
    %Animation
%     save_gif = 0;
%     animation(time(:,1),theta(:,1),w(:,1),'RK4',save_gif);

%     For Function Test (Global Truncation Error/ Simple Harmonic Motion)
%     =================================
%     w0 = (g/l)^0.5;
%     truedata = theta_0*cos(w0*duration)+w_0/w0*sin(w0*duration);
%     ans(1,1)= h;
%     ans(1,2)= theta(end,1);
%     ans(1,3)= truedata;
%     ans(1,4)= abs(theta(end,1)-truedata); %RK4
%     ans(1,5)= abs(theta(end,2)-truedata); %euler
%     ans(1,6)= truedata;
%     ans
%     =================================
    clear all
end

%====================================================================
% MATLAB ode45
function [t,y] = builtin_ode45()
    options=odeset('RelTol',5.7627e-09,'AbsTol',1.1291e-10);    
    [g,l,theta_0,w_0,~,~,~,dt,~,duration] = get_var();     
    tspan= 0:dt:duration; % set time interval
    init=[theta_0,w_0]; % set initial conditions
    [t,y]=ode45(@myode,tspan,init,options);
    figure()
    plot(t,y(:,1)/pi);
end
function dydt = myode(t,y)
    [g,l,theta_0,w_0,Fext_0,Fext_w,c,~,~,~] = get_var(); 
    dydt = [y(2); -g/l.*sin(y(1))-c*y(2)+ Fext_0.*sin(Fext_w*(t))];
end

%====================================================================
%Variables init
function [gravity,length,theta_0,omega_0,extForce_0,extForce_omega,coeffi_damp,step_size,data_pts,duration] = get_var()
    gravity = 9.8;              %ms-2
    length = 9.8;               %m
    theta_0 = 20*2*pi/360;      %rad
    omega_0 = 0*2*pi/360;       %rad/s
    extForce_0 = 0.1;           
    extForce_omega = 1/3;       %rad/s
    coeffi_damp = 0.1;         
    
    duration = 200;             %s
    step_size = 0.1*(1/2)^5;
    data_pts = duration/step_size + 1;
end

%====================================================================
% Equation of motion
function [h_f, h_g] = motion_fn()
    [g,l,~,~,extF,extF_w,c,~,~,~] = get_var();   %%Get initial parameters
    h_f = @(time,theta,w) w;
    h_g = @(time,theta,w) -g/l*sin(theta) - c*w + extF*sin(extF_w*(time));
end

%====================================================================
% Test Function (no damping, no external force, simple harmonic oscillator)
function [h_f, h_g] = test_function_1()
    [g,l,~,~,~,~,~,~,~,~] = get_var();
    h_f = @(time,theta,w) w;
    h_g = @(time,theta,w) -g/l*theta;
end 

%====================================================================
% Euler
function [time,theta,w] = euler_1(time,theta,w)
    [~,~,~,~,~,~,~,dt,~,duration] = get_var();
    [f,g] = motion_fn;
%     [f,g] = test_function_1;   %%For Function Test
    ind = 0:dt:duration;
    for i = 1:(length(ind)-1)
        %Two first order differential equations 
        w(i+1) = w(i) + dt*g(time(i),theta(i),w(i));
        theta(i+1) = theta(i)+f(time(i),theta(i),w(i))*dt;
        time(i+1) = dt*i;
    end
end

%====================================================================
% Runga-Kutta 4th order
function [time,theta,w] = RK4(time,theta,w)
    [~,~,~,~,~,~,~,dt,data_pts,duration] = get_var();
    ind = 0:dt:duration;
    error_1(1) = 0;
    for i = 1:(length(ind)-1)
        [theta(i+1),w(i+1)]= RK4_funct(time(i),theta(i),w(i),dt);
        if i>1 %%For error analysis
            [theta_2h,w_2h]= RK4_funct(time(i-1),theta(i-1),w(i-1),dt*2);
            abs_error_1(i)=theta_2h-theta(i+1);
            rel_error_1(i)=((theta_2h-theta(i+1))/min(theta_2h,theta(i+1)));
            if abs(abs_error_1(i))>1.e-3 || abs(rel_error_1(i))>1.e-3
                disp('Error too large');
                return
            end
        end
        time(i+1) = dt*i;
    end
    fprintf('Local Truncation Absolue & Relative Error: %.12f %.12f\n',max(abs(abs_error_1)),max(abs(rel_error_1)));
    %plot(error_1);
end

function [theta_new,w_new] = RK4_funct(time,theta,w,dt)
    [f, g] = motion_fn;
    %     [f,g] = test_function_1;   %%For Function Test
    k1 = f(time,          theta,                  w);
    l1 = g(time,          theta,                  w);
    k2 = f(time+dt/2,     theta+k1*dt/2,          w+l1*dt/2);
    l2 = g(time+dt/2,     theta+k1*dt/2,          w+l1*dt/2);
    k3 = f(time+dt/2,     theta+k2*dt/2,          w+l2*dt/2);
    l3 = g(time+dt/2,     theta+k2*dt/2,          w+l2*dt/2);
    k4 = f(time+dt,       theta+k3*dt,            w+l3*dt); 
    l4 = g(time+dt,       theta+k3*dt,            w+k3*dt);
    theta_new = theta + (k1 + 2*k2 + 2*k3 + k4)*dt/6;  
    w_new = w + (l1 + 2*l2 + 2*l3 + l4)*dt/6;
end


%====================================================================
% Plot Graph
function plot_graph(time,theta,w) 
    figure();
    
    subplot(2,1,1) 
    plot(time, theta/pi);
    xlabel("Time / s");
    ylabel("Theta \theta  [\pi rad]");
    title("Angular displacement");
    set(gca,'fontsize',10);
    
    subplot(2,1,2) 
    plot(time, w/pi);
    xlabel("Time / s");
    ylabel("Omega \omega  [rad.s^{-1}]");
    title("Angular velocity");
    set(gca,'fontsize',10);
    
    figure();
    plot(theta./pi,w);
    xlabel("Theta \theta  [\pi rad]");
    ylabel("Omega \omega  [rad.s^{-1}]");
    title("Phase space diagram");
    set(gca,'fontsize',10);
end

function animation(time,theta,w,gifname,save_gif)
    pause(1);    
    [~,l,~,~,~,~,~,~,data_pts,~] = get_var();
    frame_count = 1;
    figure()   
    pos = [0.1 0 0.30*1.5 0.54*1.5];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    
    for c = 1:100:data_pts
        %====================================================================
        % Phase Plot 
        subplot('position',[0.1 0.1 0.8 0.45]) 
        plot(theta./pi,w,'k','lineWidth',2);
        hold on
        cur_plot = plot(theta(c)/pi,w(c),'ro');
        set(cur_plot,'markersize',15,'markerFaceColor','g')

        xlabel('\theta  [ \pi rad ]')
        ylabel('\omega  [ rad.s^{-1}]')
        grid on
        set(gca,'fontsize',12)
        xlim([min(theta)/pi max(theta)/pi])
        %   set(gca,'xtick',-0.5:0.1:0.5)
        hold off

        %====================================================================
        % pendulum animation     
        subplot('position',[0.2 0.62 0.6 0.35])    
        X = l.*sin(theta); 
        Y = -l.*cos(theta);
        plot([0 X(c)],[0 Y(c)],'k','lineWidth',1);
        hold on

        cur_plot = plot(X(c),Y(c),'o');
        set(cur_plot,'markersize',15,'markerFaceColor','r')
        xlim([-10 10])
        ylim([-10 10])
        grid on
        set(gca,'fontsize',14) %fontsize = 14
        axis square   
        axis off
        hold off

        pause(0.001)

        %====================================================================
        % Save as gif 
        if save_gif 
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            %  On the first loop, create the file. In subsequent loops, append.
            delay = 0.2;
            if frame_count == 1
                imwrite(imind,cm,gifname,'gif','DelayTime',delay,'loopcount',inf);
            else
                imwrite(imind,cm,gifname,'gif','DelayTime',delay,'writemode','append');
            end
            frame_count = frame_count+1;
        end
	end
    hold off
end