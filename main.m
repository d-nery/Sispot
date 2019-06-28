clear;

A = 1;
B = 2;
C = 3;

V = 13800 / sqrt(3);

Vth = [V; V * exp(i * deg2rad(-120)); V * exp(i * deg2rad(120))];
Zth = (V^2) / ((10E8/3) * exp(-i * deg2rad(80))); % Impedancia equivalente do Thevenin do no da subestacao
Ith = Vth ./ Zth;
Yth = 1 / Zth; % Admitancia do equivalente de Norton do no da subestacao

Rmax = 10; % Valor maximo da resistencia de falta

TOP014;
CAR014;
VOL014;

resultado = [];
resposta = [];

trechos = size(topologia, 1);
tamanho = 3 * (trechos + 2); % Trechos + falta-gnd + falta-pai

Zcarga = (13800 ^ 2) ./ (1000 * (cargas(:, 2) - i * cargas(:, 2) .* tan(acos(cargas(:, 3))))) ./ 3;

%Loop principal para cada um dos 10 casos
for contcasos = 1:10
    E10meas = [
        complex(Emedido(contcasos, 2), Emedido(contcasos, 3)); % Leitura da fase A do caso de simulacao
        complex(Emedido(contcasos, 4), Emedido(contcasos, 5)); % Leitura da fase B do caso de simulacao
        complex(Emedido(contcasos, 6), Emedido(contcasos, 7)); % Leitura da fase C do caso de simulacao
    ];

    resultadof = Inf;

    for posicao_falta = 1:trechos
        % Varia a posicao da falta e para todos os pontos possiveis de falta cria uma matriz Q
        % Inclui o ramo pai-falta no lugar de pai-filho original
        % E inclui ramo falta-filho na penultima linha da matriz Q
        Q = zeros(tamanho);
        Q(end-2:end, end-2:end) = eye(3);

        for i = 1:trechos
            aux = 3 * i;
            no_pai   = topologia(i, 1) / 10;
            no_filho = topologia(i, 2) / 10;

            np = 3*no_pai;
            nf = 3*no_filho;

            if i == posicao_falta
                Q(aux-2:aux,   np-2 :np)  =  eye(3);
                Q(aux-2:aux,   end-2:end) = -eye(3);
                Q(end-5:end-3, nf-2 :nf)  =  eye(3);
                Q(end-5:end-3, end-2:end) = -eye(3);
            else
                Q(aux-2:aux, np-2:np) =  eye(3);
                Q(aux-2:aux, nf-2:nf) = -eye(3);
            end
        end

        % Constroi YPR. Sem levar em consideracao a falta. YPR eh posteriormente ajustada para incluir a posicao onde ocorre a falta
        YPR = zeros(tamanho);

        for i = 1:size(topologia, 1)
            % Matriz Auxiliar para casos de 1 a 5 em que eh necessario calcular a media dos valores e equilibrar a rede
            % Para casos 6 ao 10 eh igual a matriz topologia original
            topologia_aux = topologia(i, :);

            if (contcasos < 6)
                % Media das resistencias proprias
                topologia_aux([4 12 20]) = mean(topologia(i, [4 12 20]));
                % Media das reatancias proprias
                topologia_aux([5 13 21]) = mean(topologia(i, [5 13 21]));
                % Media das resistencias mutuas
                topologia_aux([6 8 10 14 16 18]) = mean(topologia(i, [6 8 10 14 16 18]));
                % Media das reatancias mutuas
                topologia_aux([7 9 11 15 17 19]) = mean(topologia(i, [7 9 11 15 17 19]));
            end

            aux = 3 * i;

            distancia = topologia_aux(3);
            YPR(aux - 2, aux - 2) =  inv(distancia * complex(topologia_aux(4), topologia_aux(5)));
            YPR(aux - 2, aux - 1) = -inv(distancia * complex(topologia_aux(6), topologia_aux(7)));
            YPR(aux - 2, aux - 0) = -inv(distancia * complex(topologia_aux(8), topologia_aux(9)));

            YPR(aux - 1, aux - 2) =  YPR((aux - 2), (aux - 1));
            YPR(aux - 1, aux - 1) =  inv(distancia * complex(topologia_aux(12), topologia_aux(13)));
            YPR(aux - 1, aux - 0) = -inv(distancia * complex(topologia_aux(14), topologia_aux(15)));

            YPR(aux - 0, aux - 2) =  YPR((aux - 2), (aux - 0));
            YPR(aux - 0, aux - 1) =  YPR((aux - 1), (aux - 0));
            YPR(aux - 0, aux - 0) =  inv(distancia * complex(topologia_aux(20), topologia_aux(21)));
        end

        % Constroi uma matriz com as cargas individuais em cada No
        matrix_cargas = zeros(tamanho);
        matrix_cargas(1, 1) = Yth;
        matrix_cargas(2, 2) = Yth;
        matrix_cargas(3, 3) = Yth;

        for i = 2:trechos
            aux = 3 * i;
            matrix_cargas(aux - 2, aux - 2) = inv(Zcarga(i - 1));
            matrix_cargas(aux - 1, aux - 1) = inv(Zcarga(i - 1));
            matrix_cargas(aux - 0, aux - 0) = inv(Zcarga(i - 1));
        end

        funcao = Inf;

        for x = 1:1:topologia(posicao_falta, 3) - 1 % Laco for que testa todos as distancias com passo de 1[m]
            for rf = 0.1:0.1:Rmax % Laco for que testa resistencias de 0.1 a Rmax [ ohms ] com passo de 0.1 [ ohm ]
                % Matriz Auxiliar para casos de 1 a 5 em que eh necessario calcular a media dos valores e equilibrar a rede
                % Para casos 6 ao 10 eh igual a matriz topologia original
                topologia_aux = topologia(posicao_falta, :);

                if (contcasos < 6)
                    % Media das resistencias proprias
                    topologia_aux([4 12 20]) = mean(topologia(posicao_falta, [4 12 20]));
                    % Media das reatancias proprias
                    topologia_aux([5 13 21]) = mean(topologia(posicao_falta, [5 13 21]));
                    % Media das resistencias mutuas
                    topologia_aux([6 8 10 14 16 18]) = mean(topologia(posicao_falta, [6 8 10 14 16 18]));
                    % Media das reatancias mutuas
                    topologia_aux([7 9 11 15 17 19]) = mean(topologia(posicao_falta, [7 9 11 15 17 19]));
                end

                % Atualiza primitiva na posicao da pai-falta
                aux = 3 * posicao_falta;
                YPR(aux - 2, aux - 2) =  inv(x * complex(topologia_aux(4), topologia_aux(5)));
                YPR(aux - 2, aux - 1) = -inv(x * complex(topologia_aux(6), topologia_aux(7)));
                YPR(aux - 2, aux - 0) = -inv(x * complex(topologia_aux(8), topologia_aux(9)));

                YPR(aux - 1, aux - 2) =  YPR((aux - 2), (aux - 1));
                YPR(aux - 1, aux - 1) =  inv(x * complex(topologia_aux(12), topologia_aux(13)));
                YPR(aux - 1, aux - 0) = -inv(x * complex(topologia_aux(14), topologia_aux(15)));

                YPR(aux - 0, aux - 2) =  YPR((aux - 2), (aux - 0));
                YPR(aux - 0, aux - 1) =  YPR((aux - 1), (aux - 0));
                YPR(aux - 0, aux - 0) =  inv(x * complex(topologia_aux(20), topologia_aux(21)));

                % Atualiza primitiva na posicao da filho-falta
                aux = tamanho - 3;
                dist = topologia_aux(3);
                YPR(aux - 2, aux - 2) =  inv((dist - x) * complex(topologia_aux(4), topologia_aux(5)));
                YPR(aux - 2, aux - 1) = -inv((dist - x) * complex(topologia_aux(6), topologia_aux(7)));
                YPR(aux - 2, aux - 0) = -inv((dist - x) * complex(topologia_aux(8), topologia_aux(9)));

                YPR(aux - 1, aux - 2) =  YPR(aux - 2, aux - 1);
                YPR(aux - 1, aux - 1) =  inv((dist - x) * complex(topologia_aux(12), topologia_aux(13)));
                YPR(aux - 1, aux - 0) = -inv((dist - x) * complex(topologia_aux(14), topologia_aux(15)));

                YPR(aux - 0, aux - 2) =  YPR(aux - 2, aux - 0);
                YPR(aux - 0, aux - 1) =  YPR(aux - 1, aux - 0);
                YPR(aux - 0, aux - 0) =  inv((dist - x) * complex(topologia_aux(20), topologia_aux(21)));

                % Atualiza primitiva da falta com sua resistencia para o GND
                matrix_cargas(end - 2, end - 2) = inv(rf);
                matrix_cargas(end - 1, end - 1) = inv(rf);
                matrix_cargas(end - 0, end - 0) = inv(rf);

                YNOS = Q.' * YPR * Q;

                YNOS = YNOS + matrix_cargas;

                % 3) Calculo das tensoes nodais
                I = zeros(tamanho, 1);
                I(1:3) = Ith;

                E = YNOS \ I;
                E10calc = E(1:3);

                % 4) Calculo da funcao a ser otimizada
                funcao_old = sum(abs(E10calc - E10meas) ./ abs(E10meas));

                % 5) Armazenamento da distancia e da resistencia de falta para o menor valor encontrado para a funcao a ser otimizada
                if funcao_old < funcao
                    xcalc = x;
                    rcalc = rf;
                    funcao = funcao_old;
                end
            end
        end

        % Salva menor caso para cada posicao de falta entre os ramos.
        resultado = [resultado; contcasos topologia_aux(1) topologia_aux(2) xcalc rcalc funcao];

        % Salva menor caso de todos
        if funcao < resultadof
            resultadof = funcao;
            variavel = [contcasos topologia_aux(1) topologia_aux(2) xcalc rcalc resultadof];
        end

        fprintf("%02.f,%03.f,%03.f,%03.f,%2.1f,%2.3f\n", contcasos, topologia_aux(1), topologia_aux(2), xcalc, rcalc, funcao);
    end

    % Salva menor caso de todos
    resposta = [resposta; variavel];
end

csvwrite('OUT081_2.csv', resultado);
csvwrite('REL081_2.csv', resposta);
