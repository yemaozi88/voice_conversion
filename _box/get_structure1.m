%
% 2008-10-21 set F1-F2 structure for CyberGlove structure
%

clear all;
cd 'H:\u-tokyo\HMTS-win\PCA\dgv_B\'

%% load basis of PCA
dirlist = dir('ex1');
dirlength = length(dirlist);
X0 = [];
i = 1;
while i < dirlength + 1
	%% except ".", "..", "DS_Store"
  	if length(dirlist(i).name) > 3 
        filename = dirlist(i).name;
        A = load_dgv(['ex1\' filename]);
        X0 = [X0; A];
    end
    i = i + 1;
end
[EVec1, EVal1, u1] = PCA(X0);
Y0 = PCA_Trans(X0, EVec1, u1, 2);


%% distribution of basic 28 gestures
hold on
fp = fopen('test.txt', 'wt');

M = [];
a = 1;
while a < 29
    b = num2str(a);
    if a < 10
        filename1 = ['0' b '_1.dgv'];
        filename2 = ['0' b '_2.dgv'];
    else
        filename1 = [b '_1.dgv'];
        filename2 = [b '_2.dgv'];
    end
    
    X1 = load_dgv(['ex1\' filename1]);
    X2 = load_dgv(['ex1\' filename2]);
    
    Y1 = PCA_Trans(X1, EVec1, u1, 2);
    Y2 = PCA_Trans(X2, EVec1, u1, 2);
    
    %plot(Y1(:,1), Y1(:,2), 'b.');
    %plot(Y2(:,1), Y2(:,2), 'b.');
    
    Y = [Y1;Y2];
    m = mean(Y);
    M = [M; m];
    plot(m(:,1),m(:,2),'+');
    fprintf(fp, '%2d %f, %f\n', a, m(:,1), m(:,2));

    a = a + 1;
end
clear filename1;
clear filename2;
clear X1;
clear X2;
clear Y;
clear m;
fclose(fp);


%% the structure for DataGlove
% F1-F2
% a = [668.050565, 1083.927686];
% i = [277.051006, 2143.158583];
% u = [285.692970, 1410.754763];
% e = [474.069421, 1807.148577];
% o = [456.793823, 777.843775];
% scep
a = [0.3089, 0.0485];
i = [-0.0398, -0.0488];
u = [-0.0044, -0.1997];
e = [0.1300, -0.2714];
o = [0.0397, 0.0057];

A = [a;i;u;e;o];

%plot(A(:,1)/8, A(:,2)/8, 'm*');

t1 = M(1,:); % origin, a
t2 = M(10,:); % target, o

vo = t2 - t1;
vi = o - a; % basis

% ratio
k = sqrt(vo(1)^2+vo(2)^2)/sqrt(vi(1)^2+vi(2)^2);

% rotation deg
rot = (1/k)*inv([vi(1), -vi(2);vi(2), vi(1)])*vo';
rot = rot';
rad = acos(rot(1));
R = [cos(rad),-sin(rad);sin(rad), cos(rad)];

ii = repmat(a, 5,1); % origin
B = (A - ii) * R * k;

tt = repmat(t1, 5,1);
T = tt + B;
plot(T(1,1), T(1,2), 'r*');
plot(T(2,1), T(2,2), 'g*');
plot(T(3,1), T(3,2), 'c*');
plot(T(4,1), T(4,2), 'm*');
plot(T(5,1), T(5,2), 'k*');

hold off