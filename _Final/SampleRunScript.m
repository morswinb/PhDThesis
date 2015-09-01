%sample (simplified )running script
clear all;
clc;

a=1;%initial time
t=[0 100]; %computation interval

epsilon=1e-3;
N=100;
%%grid choice
% grid=RegularGrid(N,epsilon);
% grid=ArcTanhGrid(N,epsilon);
grid=QuadraticGrid(N,epsilon);

%choice of a benchmark
bench=BenchmarkSimilar(a);

%object that hold all initial conditins
ic=InitialCondition(grid,bench);

%%derivative approximation choice
% diffs=FiniteDifferences(grid);
diffs=AsymFD(grid,ic);
% diffs=Polynomial2nd(grid);
% diffs=DiffSpline(grid);

% methods for dealing with BC
left=leftBC1(ic);
right=rightBC1(ic);

%setting up the sysem
sys1=CrackSystemW(ic,left,right,diffs);

%computation
tic
    [sol]=sys1.solve(t);
toc

%optional results ploting
[X,T]=meshgrid(grid.xi(1:N),sol.x);
figure
ww=bench.w(T,X);
mesh(X,T,sol.y(1:N,:)');
mesh(X,T,abs(sol.y(1:N,:)'-ww)./ww);
figure
L=sol.y(N+1,:).^.5;
plot(sol.x,abs(L-bench.L(sol.x))./bench.L(sol.x),'r')

