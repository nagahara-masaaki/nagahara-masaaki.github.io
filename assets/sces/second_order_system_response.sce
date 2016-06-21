//// 2次系のインパルス応答と過渡応答

// 初期設定
s = poly(0,'s'); // Laplace 変換の記号 s の定義

// 2次系の構成
zeta = 0.5; // 減衰係数
wn = 1; // 自然角周波数
K = 1; // ゲイン
Gs = K*wn^2/(s^2 + 2*zeta*wn*s + wn^2); // 1次系の伝達関数

// インパルス応答
T=1/wn;
t = 0:T/100:T*20; // シミュレーションの時間
gt = csim('impulse', t, Gs); // インパルス応答

// インパルス応答の描画
figure;
plot2d(t,gt,style=5)
title('2nd-order system (zeta=0.5,wn=1,K=1)')
xlabel('time (sec)')
ylabel('amplitude')

// ステップ応答
t = 0:T/100:T*20; // シミュレーションの時間
yt = csim('step', t, Gs); //ステップ応答

// ステップ応答の描画
plot2d(t,yt,style=2);
A=gca();
P=A.children.children;
P(1).thickness=3;
P(2).thickness=3;
xgrid();
legend('impulse response','step response',1);
]
// 減衰係数の違いによるステップ応答の変化
figure;
//plot(t,ones(t),'--m')
title('2nd-order system step response (zeta=0.1,0.5,1,2)')
xlabel('time (sec)')
ylabel('amplitude')
n=1;
for zeta = [0.1,0.5,1,2];
    Gs = K*wn^2/(s^2 + 2*zeta*wn*s + wn^2); // 1次系の伝達関数
    yt = csim('step', t, Gs); //ステップ応答
    plot2d(t,yt,style=n+1)
    n=n+1;
end

A=gca();
P=A.children.children;
P(1).thickness=3;
P(2).thickness=3;
P(3).thickness=3;
P(4).thickness=3;
xgrid();
legend('zeta=0.1','zeta=0.5','zeta=1','zeta=2',1);

