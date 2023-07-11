clear
clc

J_ex = 1; % 交換相互作用係数

N = 10; % 1辺のスピンの数
spin = zeros(1,3,N,N,N);
X = zeros(N,N,N);
Y = zeros(N,N,N);
Z = zeros(N,N,N);
% 初期スピン配列を生成
for i = 1:N
    for j = 1:N
        for k = 1:N
            r = 1;
            theta = pi*rand;
            phi = 2*pi*rand;
            spin(:,:,i,j,k) = [r*sin(theta)*cos(phi) r*sin(theta)*sin(phi) r*cos(theta)];
            X(i,j,k) = i;
            Y(i,j,k) = j;
            Z(i,j,k) = k;
        end
    end
end

Nstep = 1000; % 繰り返し計算を行う回数。これを大きくするとスピンがよりそろうが、時間がかかる。 
for step = 1:Nstep
    % あ
    i = randi(N);
    j = randi(N);
    k = randi(N);
    r = 1;
    theta = pi*rand;
    phi = 2*pi*rand;
    updated_spin = [r*sin(theta)*cos(phi) r*sin(theta)*sin(phi) r*cos(theta)];

    % 更新しない場合のまわりとの交換相互作用の総和を計算
    w_0 = 0;
    if i+1 <= N 
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i+1,j,k));
    end
    if i-1 >= 1
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i-1,j,k));
    end
    if j+1 <= N
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i,j+1,k));
    end
    if j-1 >= 1
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i,j-1,k));
    end
    if k+1 <= N
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i,j,k+1));
    end
    if k-1 >= 1
        w_0 = w_0 + -2*J_ex*dot(spin(:,:,i,j,k),spin(:,:,i,j,k-1));
    end

    % 更新した場合のまわりとの交換相互作用の総和を計算
    w_1 = 0;
    if i+1 <= N 
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i+1,j,k));
    end
    if i-1 >= 1
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i-1,j,k));
    end
    if j+1 <= N
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i,j+1,k));
    end
    if j-1 >= 1
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i,j-1,k));
    end
    if k+1 <= N
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i,j,k+1));
    end
    if k-1 >= 1
        w_1 = w_1 + -2*J_ex*dot(updated_spin,spin(:,:,i,j,k-1));
    end
    
    if w_1 < w_0 % 更新したほうがエネルギー的に安定ならば
        % スピンを更新
        spin(:,:,i,j,k) = updated_spin;
    end

    % プロット
    U = squeeze(spin(1,1,:,:,:));
    V = squeeze(spin(1,2,:,:,:));
    W = squeeze(spin(1,3,:,:,:));
    quiver3(X,Y,Z,U,V,W,0,Alignment='center')
    drawnow
    %exportgraphics(gcf,'result3d.gif','Append',true);
end
