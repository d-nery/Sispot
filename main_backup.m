Vth = 13800 / sqrt(3); % Fonte de tensao do equivalente de Thevenin do nao da subestao
Zth = (Vth^2) / (((1/3) * 10E8) * (cos(80 * pi / 180) - 1i * sin(80 * pi / 180))); % Impedancia equivalente do Thevenin do no da subestacao
Ith = Vth / Zth; % Fonte de corrente do equivalente de Norton do no da subestacao
Yth = 1 / Zth; % Admit^ancia do equivalente de Norton do no da subestacao
Rmax = 10; % Valor maximo da resist^encia de falta

% Dados de topologia
% <no pai ><no filho >< comprimento [m]>< resistencia do trecho [ ohms /m]>< reat^ancia do trecho [ ohms /m]>
TOP000;

% Dados das cargas
% <no onde a carga esta conectada >< pot^encia [kW] da carga >< fator de pot^encia ( indutivo )>
CAR000;

% Dados das tensoes medidas no no da subestacao para curtos - circuitos em pontos desconhecidos
% <numero do evento de curto - circuito >< parte real da tensao [V]>< parte imaginaria da tensao [V]>
VOL000;

resultado = [];
resposta = [];


%posicao_falta = 2;


ultima_posicao = size(topologia, 1) + 2;

for contcasos = 1:10
    E10Ameas = Emedido(contcasos, 2) + 1i * Emedido(contcasos, 3); % Leitura da fase A ( equivalente monofasico ) do primeiro caso de simulacao
    E10Bmeas = Emedido(contcasos, 4) + 1i * Emedido(contcasos, 5); % Leitura da fase A ( equivalente monofasico ) do primeiro caso de simulacao
    E10Cmeas = Emedido(contcasos, 6) + 1i * Emedido(contcasos, 7); % Leitura da fase A ( equivalente monofasico ) do primeiro caso de simulacao

    % Calculo das imped^ancias das cargas
    for aux = 1:size(cargas, 1)
        Zcarga(aux) = (13800^2) / (1000 * (cargas(aux, 2) - 1i * cargas(aux, 2) .* tan(acos(cargas(aux, 3)))));
    end

    resultadof = Inf;


    for posicao_falta = 1:size(topologia, 1) %varia a posicao da falta para todos os pontos possiveis

        matriz_Q = zeros(ultima_posicao, ultima_posicao);
        matriz_Q(ultima_posicao, ultima_posicao) = 1;

        for i = 1:(ultima_posicao - 2)
            if i == posicao_falta
                matriz_Q(i, (topologia(i, 1)/10)) = 1;
                matriz_Q(i, ultima_posicao) = -1;
                matriz_Q((ultima_posicao - 1), (topologia(i, 2)/10)) = 1;
                matriz_Q((ultima_posicao - 1), ultima_posicao) = -1;
            else
                matriz_Q(i, (topologia(i, 1)/10)) = 1;
                matriz_Q(i, (topologia(i, 2)/10)) = -1;
            end
        end

        YPR = zeros(ultima_posicao, ultima_posicao);
        for  i = 1:size(topologia, 1)
            YPR(i, i) = inv(topologia(i, 3) * (topologia(i, 4) + 1i * topologia(i, 5)));
        end

        matrix_cargas = zeros(ultima_posicao, ultima_posicao);
        matrix_cargas(1,1) = Yth;
        for  i = 2:(ultima_posicao-1)
            matrix_cargas (i, i) = inv(Zcarga(i-1));
        end


        funcao = Inf;
        

        for  x = 1:1:topologia(posicao_falta, 3) - 1 % Laco for que testa todos as distancias com passo de 1[m]

            for rf = 0.1:0.1:Rmax % Laco for que testa resistencias de 0.1 a Rmax [ ohms ] com passo de 0.1 [ ohm ]

                YPR(posicao_falta, posicao_falta) = inv(x * (topologia(posicao_falta, 4) + 1i * topologia(posicao_falta, 5)));
                YPR(ultima_posicao - 1, ultima_posicao - 1) = inv((topologia(posicao_falta, 3) - x) * (topologia(posicao_falta, 4) + 1i * topologia(posicao_falta, 5)));
                YPR(ultima_posicao, ultima_posicao) = inv(rf);
                YNOS = transpose(matriz_Q) * YPR * matriz_Q;

                YNOS = YNOS + matrix_cargas;

                % 3) Calculo das tensoes nodais
                E = YNOS \ [Ith; 0; 0; 0; 0];
                E10calc = E(1);
                % 4) Calculo da funcao a ser otimizada
                % 5) Armazenamento da distancia e da resistencia de falta para o menor valor encontrado para a funcao a ser otimizada
                funcao_old = abs(E10Ameas - E10calc) / abs(E10Ameas);

                if funcao_old < funcao
                    xcalc = x;
                    rcalc = rf;
                    funcao = funcao_old;
                end
            end
        end
        linha = [contcasos topologia(posicao_falta, 1) topologia(posicao_falta, 2) xcalc rcalc funcao];
        resultado(end+1, :) = linha;

        if funcao < resultadof 
            resultadof = funcao;
            variavel = [contcasos topologia(posicao_falta, 1) topologia(posicao_falta, 2) xcalc rcalc resultadof];
        end

        fprintf("%02.f,%03.f,%03.f,%03.f,%2.1f,%2.3f\n", contcasos, topologia(posicao_falta, 1), topologia(posicao_falta, 2), xcalc, rcalc, funcao);
    end

    resposta(end+1, :) = variavel; 

end

csvwrite('OUT000CERTO.csv',resultado);
csvwrite('RELBIGOLA.csv',resposta);