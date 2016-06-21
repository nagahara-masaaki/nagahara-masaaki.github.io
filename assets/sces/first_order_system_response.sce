//// 1次系のインパルス応答と過渡応答

// 初期設定
s = poly(0,'s'); // Laplace 変換の記号 s の定義

// 1次系の構成
T = 10; // 時定数
K = 1; // ゲイン
Gs = K/(T*s+1); // 1次系の伝達関数

// インパルス応答
t = 0:T/100:T*5; // シミュレーションの時間
gt = csim('impulse', t, Gs); // インパルス応答

// インパルス応答の描画
figure;
plot2d(t,gt,style=5)
title('1st-order system')
xlabel('time (sec)')
ylabel('amplitude')

// ステップ応答
t = 0:T/100:T*5; // シミュレーションの時間
yt = csim('step', t, Gs);//ステップ応答

// ステップ応答の描画
plot2d(t,yt,style=2);

A=gca();
P=A.children.children;
P(1).thickness=3;
P(2).thickness=3;
xgrid();
legend('impulse response','step response',1);
