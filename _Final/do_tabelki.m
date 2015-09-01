% G?ówny skrypt do odpalenia
clear all;
clear java;
% javaaddpath('C:/Users/Morswin/Desktop/Workspace/CodeForMATLAB/bin');
clc;
format long;

k=1;
M=1;
a=1;

% Leak=LeakOff.LeakLoader.loadCarterLinear();

%wersja no leak off
% N=20
% epsilon=1e-6;
% grid=Nowa_Siatka_2(N,epsilon,1,5);
% bench=Benchmark_Similar(a);
% wersja carter leak off
N=100
% 
epsilon=1e-7;
grid=QuadraticGrid(N,epsilon);

% epsilon=1e-3;
% grid=RegularGrid(N,epsilon);

bench=Benchmark_CC_2(k,M,a,2,1); %ratio 0.9
% bench=Benchmark_CC_2(k,M,a,1,3); %ratio 0.5

% epsilon=4*10^-4/N^2;
% t=0:1:100;
t=[0 100];
% t=0:.00001:.0001;
% t=0:1:100;


%wybór siatki
% grid=RegularGrid(N,epsilon);
%grid=RegularGrid_2(N,epsilon);
% grid=ArcTanhGrid(N,epsilon);
% grid=QuadraticGrid(N,epsilon);
%grid=MixedGrid(N,epsilon);
% grid=Nowa_Siatka(N,epsilon);
% grid=QuadraticGrid_2(N,epsilon,2.5);



% bench=BenchmarkSolution2(k,M,a);



%opcja dyskretyzacji
diffs=FiniteDifferences(grid);
% diffs=Diffrences.FD(N,grid.xi);
% diffs=Diffrences.FDloader.loadFD(N,grid.xi);
% diffs=Polynomial2nd(grid);
% diffs=DiffSpline(grid);


% bench=Benchmark_Custom(k,M,a);
% bench=Benchmark_Custom_2(k,M,a);
% bench=Benchmark_Carter(k,M,a);
% bench=Benchmark_Presure_Carter(k,M,a);
% bench=Benchmark_C_Presure(k,M,a);
% bench=BenchmarkSolution2(k,M,a);
% bench=Benchmark_CC(k,M,a);

% bench=Benchmark2(k,M,a,1,1);
% bench=Benchmark_Similar(a);

% bench=BenchmarkSolution2(k,M,a);
% bench=BenchmarkSolution1(k,M,a,grid);
% bench=Benchmark_filmowy_3(); 

ic=Initial_Condition(grid,bench);
% ic.q_0=@(t)1-min(t,1+10^-14);
% ic.q_0=@(t)ic.q_0(t);
% ic.q_0=@(t)-1+0*t;
% ic.q_0=@(t)-0.75+min(.75,t);
% ic.q_0=@(t)ic.q_0(t)-1.2*ic.q_0(t)*min(1,(t-a));
% ic.q_0=@(t)ic.q_0(t)*cos(t/4);
% ic.q_0=@(t) t*0;

% leak off handle
% Leak=Leak_handle_Presure_Carter(ic);
% Leak=Leak_handle_Carter(ic);
% Leak=Leak_handle_C_presure(ic);



%sposób radzenia sobie z L
% L_handle=dL_handle(ic);
L_handle=dL_2_handle(ic);

% warunek na x=0, sposób liczenia tam pochodnej
left_handle=left_handle_1(ic);
% left_handle=left_handle_1_long_open_pipe(ic);
% left_handle=left_handle_2(ic);

% warunek na x=1, sposób radzenia sobie tam z pochodna i w_0

%asymptotyka 2 czlony
right_handle=right_handle_1(ic);
% asymptotyka 1  cz?on
% right_handle=right_handle_2(ic);
% newton 1 punkt 2 czlony asymptotyki
% right_handle=right_handle_3(ic);
% warunek 2 na fd
% right_handle=right_handle_4(ic);
% warunek 1 na fd
% right_handle=right_handle_5(ic);


%liczenie uk?adu
%carter
% sys1=Crack_System_W_3(ic,diffs,0,1,1,1,1);

% %leak z benchmarka
sys1=Crack_System_W_3(ic,diffs,0,1,3,1,1);


% sys1=Crack_System_W_Leak_off_nonzero(ic,Leak,L_handle,left_handle,right_handle,diffs);

% profile on
tic
% [w,L,sol]=sys1.solve(t);
    [sol]=sys1.solve(t);
toc

[X,T]=meshgrid(grid.xi(1:N),sol.x);
figure
ww=bench.w(T,X);
mesh(X,T,sol.y(1:N,:)');
mesh(X,T,abs(sol.y(1:N,:)'-ww)./ww);

fprintf('wrel %e\n',max(max(abs(sol.y(1:N,:)'-ww)./ww)));
fprintf('wabs %e\n',max(max(abs(sol.y(1:N,:)'-ww))));
L=sol.y(N+1,:).^.5;
L_rel=abs((L-bench.L(sol.x))./L);
fprintf('Lrel %e\n',max(L_rel));

% profile off
% profile viewer

% figure
% plot(sol.x,sol.y.^(.5));

% figure
% L=sol.y(end,:).^.5;
% loglog(sol.x+a,L);
% % loglog(KT,L)
% hold on
% L_small=bench.L(sol.x);
% loglog(sol.x+a,L_small,'r');
% 
% LL=2/pi*sqrt(sol.x+a);
% loglog(sol.x+a,LL,'g');
% hold off;

% plot(sol.x,sol.y(end,:).^.5)

% [X,T]=meshgrid(grid.xi,sol.x);
% figure
% ww=bench.w(T,X);
% % mesh(X,T,abs(sol.y(1:N,:)'-ww)./ww);
% % mesh(X,T,sol.y(1:N,:)');
% %  plot(sol.x,sol.y(end,:))
%  plot(sol.x,sol.y)

% max(sol.y(1:N-1,end))
% 
% figure
% mesh(grid.xi(2:N),sol.x,(sol.y(1:N-1,:)'));
% figure
% % plot(sol.x,abs(sol.y(N+1,:).^(.5)));
% plot(sol.x,sol.y(end,:).^(.5));

% Pump(1)=0;
% for i=1:length(sol.x)
%     V(i)=trapz(grid.xi,sol.y(1:N,i));
%     V(i)=V(i)*sol.y(end,i).^.5;
%     
%     if(i>1)
%         Pump(i)=Pump(i-1)+quad(@(t)ic.q_0(t),sol.x(i-1),sol.x(i));
%     end
% %     V(i)-V(1)
% %     Pump(i)
% %     pause
% end
% % figure
% % plot(sol.x,V);
% figure 
% plot(sol.x,abs(abs(V)-abs(Pump+V(1)))./(Pump+V(1)));

% % % % for i=1:length(sol.x)
% % % %     dwdx(1:N,i)=diffs.calcSecondDer(sol.y(1:N,i));
% % % % end
% % % % 
% % % % figure
% % % % mesh(X,T,dwdx');
    
% figure
% plot(sol.x,Pump)
% plot(sol.x,);
% figure
% plot(sol.x,V);
% figure
% plot(sol.x,Pump);

% L=sol.y(end,:).^.5;
% figure
% plot(sol.x,abs(L-bench.L(sol.x))./bench.L(sol.x));




% [X,T]=meshgrid(grid.xi,sol.x);
% wa=bench.w(T,X);
% 
% figure
% mesh(X,T,abs((wa-sol.y(1:N,:)')./wa));
% % zlim([0 3*10^-3]);
% set(gca,'fontsize',28)
% 
% figure
% mesh(X,T,abs((wa-sol.y(1:N,:)')));
% set(gca,'fontsize',28)

% mesh(X(1:5:250),T(1:2:100),abs((wa(1:5:250,1:2:100)-wn(1:5:250,1:2:100))./wa(1:5:250,1:2:100)));
% surf(X(1:5:100,1:10:250),T(1:5:100,1:10:250),abs((wa(1:5:100,1:10:250)-wn(1:5:100,1:10:250))./wa(1:5:100,1:10:250)),'FaceLighting','phong')
% max(max(abs((wa-wn)./wa)))


% [X,T]=meshgrid(grid.xi,sol.x);
% wa=bench.w(T,X);
% 
% NN=length(sol.x);
% sn=2;
% st=2;
% 
% figure
% mesh(grid.xi([1:sn:N N]),sol.x(1:st:NN),(abs(wa(1:st:NN,[1:sn:N N])-sol.y([1:sn:N N],1:st:NN)')...
%     ./abs(wa(1:st:NN,[1:sn:N N]))))
% 
% % mesh(grid.xi([1:sn:N N]),sol.x(1:st:NN),wa(1:st:NN,[1:sn:N N]));
% 
% 
% % surf(X(1:5:100,1:5:N),T(1:5:100,1:5:N),...
% %     abs((abs(wa(1:5:100,1:5:N))-abs(wn(1:5:100,1:5:N)))./wa(1:5:100,1:5:N)));
% % surf(X([1:5:100 100],[1:5:N 97 98 99 N]),T([1:5:100 100],[1:5:N 97 98 99 N]),...
% %      abs((abs(wa([1:5:100 100],[1:5:N  97 98 99 N]))-abs(wn([1:5:100 100],[1:5:N  97 98 99 N])))./wa([1:5:100 100],[1:5:N  97 98 99 N])));
% hold on
% figureHandle = gcf;
% zlim([0 3*10^-3]);
% set(findall(figureHandle,'type','text'),'fontSize',17,'fontWeight','bold')
% hold off
% set(gca,'fontsize',18)
% ylabel('$t$','Interpreter','LaTex','FontSize',26);
% xlabel('$x$','Interpreter','LaTex','FontSize',26);
% zlabel('$\delta w$','Interpreter','LaTex','FontSize',26);
% % set(get(gca,'ZLabel'),'Rotation',eps);
% % xlabh = get(gca,'ZLabel');
% % set(xlabh,'Position',get(xlabh,'Position') - [0 -10 0])
% 
% 
% 
% 
% 
% 
% figure
% mesh(grid.xi([1:sn:N N]),sol.x(1:st:NN),(abs(wa(1:st:NN,[1:sn:N N])-sol.y([1:sn:N N],1:st:NN)')...
%     ./1))
% % surf(X(1:5:100,1:5:N),T(1:5:100,1:5:N),...
% %     abs((abs(wa(1:5:100,1:5:N))-abs(wn(1:5:100,1:5:N)))));
% % mesh(X,T,...
% %     abs((abs(wa)-abs(wn))));
% % surf(X([1:5:100 100],[1:5:N N]),T([1:5:100 100],[1:5:N N]),...
% %     abs((abs(wa([1:5:100 100],[1:5:N N]))-abs(wn([1:5:100 100],[1:5:N N])))));
% hold on
% figureHandle = gcf;
% set(findall(figureHandle,'type','text'),'fontSize',17,'fontWeight','bold')
% hold off
% set(gca,'fontsize',20)
% ylabel('$t$','Interpreter','LaTex','FontSize',26);
% xlabel('$x$','Interpreter','LaTex','FontSize',26);
% zlabel('$\Delta w$','Interpreter','LaTex','FontSize',26);
% % set(get(gca,'ZLabel'),'Rotation',eps);
% % xlabh = get(gca,'ZLabel');
% % set(xlabh,'Position',get(xlabh,'Position') - [-0.05 -13 0])



% plot(t,L)
% plot(t,abs(L-bench.L(t)))
% plot(sol.x,abs(L-bench.L(sol.x))./bench.L(sol.x),'r')

% [X,T]=meshgrid(grid.xi,sys1.L_history(1,2:end));
% mesh(X,T,sys1.f);

% 
% figure
% L=sol.y(N+1,:).^.5;
% plot(sol.x,abs(L-bench.L(sol.x))./bench.L(sol.x),'r')

