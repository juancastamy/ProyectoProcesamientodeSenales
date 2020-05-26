%% Filtros
% Recuperación de datos 

[x, Fs] = audioread('A Beautiful Lie (Acoustic).wav');

ch1_x = x(:,1);
ch2_x = x(:,2);
% Definición del Periodo de Muestreo
Ts = 1/Fs;
% Reproducción de la entrada
centered_x = x - repmat(mean(x), size(x,1), 1);
maxx = max(abs(centered_x(:)));
scaled_centered_x = centered_x ./ maxx;
player_x = audioplayer (scaled_centered_x, Fs);
%play (player_x);
%stop (player_x);
%% Diseño de Filtro 1 
% Frecuencia y período de muestreo
Fl_1 = 500;
Fh_1 = 650;
order = 8;
[b, a] = butter ((order/2), [Fl_1 Fh_1]/(Fs/2));
% Se observa la función de transferencia
Hf = tf(b, a, Ts);
% Se observa el filtro diseñado
fvtool(b,a,'Fs',Fs);
% Aplicación de Ecuación de Diferencias
N = length(x);
ch1_y_1 = zeros(N, 1);  
ch2_y_1 = zeros(N, 1);
% Definición de variables recurrentes
% Recurrencias de x para el primer canal del filto 1
xn1_11 = 0;
xn2_11 = 0;
xn3_11 = 0;
xn4_11 = 0;
xn5_11 = 0;
xn6_11 = 0;
xn7_11 = 0;
xn8_11 = 0;
% Recurrencias de x para el primer canal del filto 2
xn1_12 = 0;
xn2_12 = 0;
xn3_12 = 0;
xn4_12 = 0;
xn5_12 = 0;
xn6_12 = 0;
xn7_12 = 0;
xn8_12 = 0;
% Recurrencias de y para el primer canal del filto 1
yn1_11 = 0;
yn2_11 = 0;
yn3_11 = 0;
yn4_11 = 0;
yn5_11 = 0;
yn6_11 = 0;
yn7_11 = 0;
yn8_11 = 0;
% Recurrencias de y para el primer canal del filto 2
yn1_12 = 0;
yn2_12 = 0;
yn3_12 = 0;
yn4_12 = 0;
yn5_12 = 0;
yn6_12 = 0;
yn7_12 = 0;
yn8_12 = 0;
% Ciclo for para muestrear los datos de la ecuación de diferencias
%*************************************************************************
%------------------------ PRIMER CANAL -----------------------------------
%*************************************************************************
for n = 1:N
    ch1_y_1(n) = b(1)*ch1_x(n) + b(2)*xn1_11 + b(3)*xn2_11 + b(4)*xn3_11 + ....
        b(5)*xn4_11 +b(6)*xn5_11 + b(7)*xn6_11 + b(8)*xn7_11 + b(9)*xn8_11 - ....
        a(2)*yn1_11 - a(3)*yn2_11 - a(4)*yn3_11 - a(5)*yn4_11 - a(6)*yn5_11 - ...
        a(7)*yn6_11 - a(8)*yn7_11 - a(9)*yn8_11;
    % Actualización de variables de x
    xn9 = xn8_11;
    xn8_11 = xn7_11;
    xn7_11 = xn6_11;
    xn6_11 = xn5_11;
    xn5_11 = xn4_11;
    xn4_11 = xn3_11;
    xn3_11 = xn2_11;
    xn2_11 = xn1_11;
    xn1_11 = ch1_x(n);
    % Actualización de variables de y
    yn8_12 = yn7_12;
    yn7_12 = yn6_12;
    yn6_12 = yn5_12;
    yn5_12 = yn4_12;
    yn4_12 = yn3_12;
    yn3_12 = yn2_12;
    yn2_12 = yn1_12;
    yn1_12 = ch2_y_1(n);
end
%*************************************************************************
%------------------------ SEGUNDO CANAL ----------------------------------
%*************************************************************************
for n = 1:N
    ch2_y_1(n) = b(1)*ch1_x(n) + b(2)*xn1_12 + b(3)*xn2_12 + b(4)*xn3_12 + ....
        b(5)*xn4_12 +b(6)*xn5_12 + b(7)*xn6_12 + b(8)*xn7_12 + b(9)*xn8_12 - ....
        a(2)*yn1_12 - a(3)*yn2_12 - a(4)*yn3_12 - a(5)*yn4_12 - a(6)*yn5_12 - ...
        a(7)*yn6_12 - a(8)*yn7_12 - a(9)*yn8_12;
    % Actualización de variables de x
    xn9 = xn8_12;
    xn8_12 = xn7_12;
    xn7_12 = xn6_12;
    xn6_12 = xn5_12;
    xn5_12 = xn4_12;
    xn4_12 = xn3_12;
    xn3_12 = xn2_12;
    xn2_12 = xn1_12;
    xn1_12 = ch2_x(n);
    % Actualización de variables de y
    yn8_12 = yn7_12;
    yn7_12 = yn6_12;
    yn6_12 = yn5_12;
    yn5_12 = yn4_12;
    yn4_12 = yn3_12;
    yn3_12 = yn2_12;
    yn2_12 = yn1_12;
    yn1_12 = ch2_y_1(n);
end
% Concatenación de ambos canales para obtener una señal stereo
y_1 = [ch1_y_1, ch2_y_1];
% Reproducción de la salida
centered_y_1 = y_1 - repmat(mean(y_1), size(y_1,1), 1);
maxy_1 = max(abs(centered_y_1(:)));
scaled_centered_y_1 = centered_y_1 ./ maxy_1;
player_y_1 = audioplayer (scaled_centered_y_1, Fs);
%play (player_y_1);
%stop (player_y_1);
%% Graficas de filtrado para Primera Señal
%Canal 1 de x
ch1_x = x(:,1);
yX1_x1 = fft(ch1_x);
yX1_bi_x1 = abs(yX1_x1/length(ch1_x));
yX1_uni_x1 = yX1_bi_x1(1:length(ch1_x)/2+1);
yX1_uni_x1(2:end-1) = 2*yX1_uni_x1(2:end-1);
axis_x_2 = Fs*(0:(length(ch1_x)/2))/length(ch1_x);
%Canal 2 de x
ch2_x = x(:,2);
yX1_x2 = fft(ch2_x);
yX1_bi_x2 = abs(yX1_x2/length(ch2_x));
yX1_uni_x2 = yX1_bi_x2(1:length(ch2_x)/2+1);
yX1_uni_x2(2:end-1) = 2*yX1_uni_x2(2:end-1);
%Canal 1 de y_1
ch1_y_1 = y_1(:,1);
yX1_y1_1 = fft(ch1_y_1);
yX1_bi_y1_1 = abs(yX1_y1_1/length(ch1_y_1));
yX1_uni_y1_1 = yX1_bi_y1_1(1:length(ch1_y_1)/2+1);
yX1_uni_y1_1(2:end-1) = 2*yX1_uni_y1_1(2:end-1);
%Canal 2 de y_1
ch2_y_1 = y_1(:,2);
yX1_y2_1 = fft(ch2_y_1);
yX1_bi_y2_1 = abs(yX1_y2_1/length(ch2_y_1));
yX1_uni_y2_1 = yX1_bi_y2_1(1:length(ch2_y_1)/2+1);
yX1_uni_y2_1(2:end-1) = 2*yX1_uni_y2_1(2:end-1);
%Gráfica
figure(2);
clf;
subplot(2, 2, 1);
plot(axis_x_2(1:44100), yX1_uni_x1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 1 de x.');
subplot(2, 2, 2);
plot(axis_x_2(1:44100), yX1_uni_x2(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 2 de x.');
subplot(2, 2, 3);
plot(axis_x_2(1:44100), yX1_uni_y1_1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 1 de y_1.');
subplot(2, 2, 4);
plot(axis_x_2(1:44100), yX1_uni_y2_1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 2 de y_1.');
%% Diseño de Filtro 2 
% Frecuencia y período de muestreo
Fl_2 = 73;
Fh_2 = 150;
order = 8;
[b2, a2] = butter ((order/2), [Fl_2 Fh_2]/(Fs/2));
% Se observa la función de transferencia
Hf2 = tf(b2, a2, Ts);
% Se observa el filtro diseñado
fvtool(b2,a2,'Fs',Fs);
% Aplicación de Ecuación de Diferencias
N2 = length(x);
ch1_y_2 = zeros(N2, 1);  
ch2_y_2 = zeros(N2, 1);
%Definición de variables recurrentes
% Recurrencias de x para el primer canal del filto 2
xn1_21 = 0;
xn2_21 = 0;
xn3_21 = 0;
xn4_21 = 0;
xn5_21 = 0;
xn6_21 = 0;
xn7_21 = 0;
xn8_21 = 0;
% Recurrencias de x para el segundo canal del filto 2
xn1_22 = 0;
xn2_22 = 0;
xn3_22 = 0;
xn4_22 = 0;
xn5_22 = 0;
xn6_22 = 0;
xn7_22 = 0;
xn8_22 = 0;
% Recurrencias de y para el primer canal del filto 2
yn1_21 = 0;
yn2_21 = 0;
yn3_21 = 0;
yn4_21 = 0;
yn5_21 = 0;
yn6_21 = 0;
yn7_21 = 0;
yn8_21 = 0;
% Recurrencias de y para el segundo canal del filto 2
yn1_22 = 0;
yn2_22 = 0;
yn3_22 = 0;
yn4_22 = 0;
yn5_22 = 0;
yn6_22 = 0;
yn7_22 = 0;
yn8_22 = 0;
% Ciclos for para muestrear los datos de la ecuación de diferencias
%*************************************************************************
%------------------------ PRIMER CANAL -----------------------------------
%*************************************************************************
for n = 1:N2
    % Ecuación de Diferencias
    ch1_y_2(n) = (b2(1)*ch1_x(n) + b2(2)*xn1_21 + b2(3)*xn2_21 + b2(4)*xn3_21 + ....
        b2(5)*xn4_21 +b2(6)*xn5_21 + b2(7)*xn6_21 + b2(8)*xn7_21 + b2(9)*xn8_21 - ....
        a2(2)*yn1_21 - a2(3)*yn2_21 - a2(4)*yn3_21 - a2(5)*yn4_21 - a2(6)*yn5_21 - ...
        a2(7)*yn6_21 - a2(8)*yn7_21 - a2(9)*yn8_21);
    % Actualización de variables de x
    xn9 = xn8_21;
    xn8_21 = xn7_21;
    xn7_21 = xn6_21;
    xn6_21 = xn5_21;
    xn5_21 = xn4_21;
    xn4_21 = xn3_21;
    xn3_21 = xn2_21;
    xn2_21 = xn1_21;
    xn1_21 = ch1_x(n);
    % Actualización de variables de y
    yn8_21 = yn7_21;
    yn7_21 = yn6_21;
    yn6_21 = yn5_21;
    yn5_21 = yn4_21;
    yn4_21 = yn3_21;
    yn3_21 = yn2_21;
    yn2_21 = yn1_21;
    yn1_21 = ch1_y_2(n);
end
%*************************************************************************
%------------------------ SEGUNDO CANAL ----------------------------------
%*************************************************************************
for n = 1:N2
    ch2_y_2(n) = (b2(1)*ch1_x(n) + b2(2)*xn1_22 + b2(3)*xn2_22 + b2(4)*xn3_22 + ....
        b2(5)*xn4_22 +b2(6)*xn5_22 + b2(7)*xn6_22 + b2(8)*xn7_22 + b2(9)*xn8_22 - ....
        a2(2)*yn1_22 - a2(3)*yn2_22 - a2(4)*yn3_22 - a2(5)*yn4_22 - a2(6)*yn5_22 - ...
        a2(7)*yn6_22 - a2(8)*yn7_22 - a2(9)*yn8_22);
    % Actualización de variables de x
    xn9 = xn8_22;
    xn8_22 = xn7_22;
    xn7_22 = xn6_22;
    xn6_22 = xn5_22;
    xn5_22 = xn4_22;
    xn4_22 = xn3_22;
    xn3_22 = xn2_22;
    xn2_22 = xn1_22;
    xn1_22 = ch2_x(n);
    % Actualización de variables de y
    yn8_21 = yn7_21;
    yn7_21 = yn6_21;
    yn6_21 = yn5_21;
    yn5_21 = yn4_21;
    yn4_21 = yn3_21;
    yn3_21 = yn2_21;
    yn2_21 = yn1_21;
    yn1_21 = ch1_y_2(n);
end
% Concatenación de ambos canales para crear una señal de salida stereo
y_2 = [ch1_y_2, ch2_y_2];
% Reproducción de la salida
centered_y_2 = y_2 - repmat(mean(y_2), size(y_2,1), 1);
maxy_2 = max(abs(centered_y_2(:)));
scaled_centered_y_2 = centered_y_2 ./ maxy_2;
player_y_2 = audioplayer (scaled_centered_y_2, Fs);

%play (player_y_2);
stop (player_y_2);
%% Graficas de filtrado para Segunda Señal
%Canal 1 de x
ch1_x = x(:,1);
yX1_x1 = fft(ch1_x);
yX1_bi_x1 = abs(yX1_x1/length(ch1_x));
yX1_uni_x1 = yX1_bi_x1(1:length(ch1_x)/2+1);
yX1_uni_x1(2:end-1) = 2*yX1_uni_x1(2:end-1);
axis_x_2 = Fs*(0:(length(ch1_x)/2))/length(ch1_x);
%Canal 2 de x
ch2_x = x(:,2);
yX1_x2 = fft(ch2_x);
yX1_bi_x2 = abs(yX1_x2/length(ch2_x));
yX1_uni_x2 = yX1_bi_x2(1:length(ch2_x)/2+1);
yX1_uni_x2(2:end-1) = 2*yX1_uni_x2(2:end-1);
%Canal 1 de y_1
ch1_y_2 = y_2(:,1);
yX1_y1_1 = fft(ch1_y_2);
yX1_bi_y1_1 = abs(yX1_y1_1/length(ch1_y_2));
yX1_uni_y1_1 = yX1_bi_y1_1(1:length(ch1_y_2)/2+1);
yX1_uni_y1_1(2:end-1) = 2*yX1_uni_y1_1(2:end-1);
%Canal 2 de y_1
ch2_y_2 = y_2(:,2);
yX1_y2_1 = fft(ch2_y_2);
yX1_bi_y2_1 = abs(yX1_y2_1/length(ch2_y_2));
yX1_uni_y2_1 = yX1_bi_y2_1(1:length(ch2_y_2)/2+1);
yX1_uni_y2_1(2:end-1) = 2*yX1_uni_y2_1(2:end-1);
%Gráfica
figure(4);
clf;
subplot(2, 2, 1);
plot(axis_x_2(1:44100), yX1_uni_x1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 1 de x.');
subplot(2, 2, 2);
plot(axis_x_2(1:44100), yX1_uni_x2(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 2 de x.');
subplot(2, 2, 3);
plot(axis_x_2(1:44100), yX1_uni_y1_1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 1 de y_2.');
subplot(2, 2, 4);
plot(axis_x_2(1:44100), yX1_uni_y2_1(1:44100));
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro unilateral del canal 2 de y_2.');