%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETUDE D'UN CODAGE CONVOLUTIF
% COMPARAISON CHAINE CODEE / CHAINE NON CODEE
% RABEFANIRAKA Yan, Mai 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETRES GENERAUX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fe = 12000;       %Fréquence d'échantillonnage
Te = 1/Fe;        %Période d'échantillonnage
Rb = 3000;        %Débit binaire souhaité
N = 1000;         %Nombre de bits générés

M = 2;         %Ordre de la modulation          
Rs = Rb / log2(M);         %Débit symbole   
Ns = 1 / (Rs * Te);         %Facteur de suréchantillonnage

%tableau des valeurs de SNR par bit souhaité à l'entrée du récpeteur en dB
tab_Eb_N0_dB = [0: 6]; 
%Passage au SNR en linéaire
tab_Eb_N0 = 10 .^ (tab_Eb_N0_dB / 10);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETRES DU CODE CONVOLUTIF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Générateurs du code convolutif (taux 1/2, contrainte K=3)
g_bin = [1 0 1; 1 1 1]; % Représentation binaire
g_oct = [5 7];          % Représentation octale
g_poly = {'1 + D^2', '1 + D + D^2'}; % Représentation polynomiale

treillis = poly2trellis(3, g_oct); % Création du treillis pour Viterbi

%codage convolutif: {0,1}^N -> {0,1}^{2N}

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUCLE SUR LES NIVEAUX DE Eb/N0 A TESTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indice_bruit = 1:length(tab_Eb_N0_dB)

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VALEUR DE Eb/N0 TESTEE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eb_N0_dB = tab_Eb_N0_dB(indice_bruit)
    Eb_N0 = tab_Eb_N0(indice_bruit);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALISATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nb_erreurs = 0;   %Variable permettant de compter le nombre d'erreurs cumulées
    nb_cumul = 0;     %Variables permettant de compter le nombre de cumuls réalisés
    TES_BPSK = 0;
    TEB_BPSK = 0;
    
    TES_BPSK_conv = 0;
    TEB_BPSK_conv_dur = 0;
    TEB_BPSK_conv_souple = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BOUCLE POUR PRECISION TES ET TEBS MESURES :COMPTAGE NOMBRE ERREURS
    % (voir annexe texte TP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(nb_erreurs < 500)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %GENERATION DE L'INFORMATION BINAIRE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits = randi([0, 1], 1, N);

        % Codage convolutif (taux 1/2, K=3)
        bits_conv = convenc(bits, treillis);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Mapping BPSK
        symboles_BPSK = pskmod(bits, M, 0, InputType="bit", PlotConstellation=false);
        a_k = real(symboles_BPSK);
        b_k = zeros(size(a_k));

        %Mapping BPSK avec codage convolutif
        symboles_BPSK_conv = pskmod(bits_conv, M, 0, InputType="bit", PlotConstellation=false);
        a_k_conv = real(symboles_BPSK_conv);
        b_k_conv = zeros(size(a_k_conv));
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SURECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        somme_Diracs_ponderes_BPSK = kron(a_k, [1 zeros(1, Ns-1)]);
        somme_Diracs_ponderes_BPSK_conv = kron(a_k_conv, [1 zeros(1, Ns-1)]);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %FILTRAGE DE MISE EN FORME 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = ones(1, Ns);

        Signal_emis_BPSK = filter(h, 1, somme_Diracs_ponderes_BPSK);
        Signal_emis_BPSK_conv = filter(h, 1, somme_Diracs_ponderes_BPSK_conv);
 
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CANAL DE PROPAGATION AWGN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %POUR MODULATION BPSK
        P_signal = (1 / (length(Signal_emis_BPSK))) * sum(norm(Signal_emis_BPSK).^2);
        P_bruit = (P_signal * Ns) / (2 * log2(M) * Eb_N0);
        Bruit = sqrt(P_bruit) * randn(1, length(Signal_emis_BPSK));
        Signal_recu_BPSK = Signal_emis_BPSK + Bruit;

        %POUR MODULATION BPSK avec codage convolutif
        P_signal_conv = (1 / (length(Signal_emis_BPSK_conv))) * sum(norm(Signal_emis_BPSK_conv).^2);
        P_bruit_conv = (P_signal_conv * Ns) / (2 * log2(M) * Eb_N0);
        Bruit_conv = sqrt(P_bruit_conv) * randn(1, length(Signal_emis_BPSK_conv));
        Signal_recu_BPSK_conv = Signal_emis_BPSK_conv + Bruit_conv;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %FILTRAGE DE RECEPTION ADAPTE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hr = ones(1, Ns);

        Signal_recu_filtre_BPSK = filter(hr, 1, Signal_recu_BPSK);
        Signal_recu_filtre_BPSK_conv = filter(hr, 1, Signal_recu_BPSK_conv);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n0 = 4;     % A COMPLETER
        Signal_echantillonne_BPSK = Signal_recu_filtre_BPSK(n0:Ns:end);
        Signal_echantillonne_BPSK_conv = Signal_recu_filtre_BPSK_conv(n0: Ns: end);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DECISIONS SUR LES SYMBOLES 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        symboles_recus_decimaux_BPSK = real(Signal_echantillonne_BPSK) > 0;
        symboles_recus_decimaux_BPSK = symboles_recus_decimaux_BPSK * 2 - 1;
        
        symboles_recus_decimaux_BPSK_conv = real(Signal_echantillonne_BPSK_conv) > 0;
        symboles_recus_decimaux_BPSK_conv = symboles_recus_decimaux_BPSK_conv * 2 - 1;

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CALCUL DU TAUX D'ERREUR SYMBOLE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TES_BPSK = TES_BPSK + sum(real(symboles_BPSK) ~= symboles_recus_decimaux_BPSK) / length(real(symboles_BPSK));
        TES_BPSK_conv = TES_BPSK_conv + sum(real(symboles_BPSK_conv) ~= symboles_recus_decimaux_BPSK_conv) / length(real(symboles_BPSK_conv));

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DEMAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits_recus_BPSK = pskdemod(symboles_recus_decimaux_BPSK, M);

        bits_recus_BPSK_conv = pskdemod(symboles_recus_decimaux_BPSK_conv, M);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DECODAGE VITERBI DUR (distance de Hamming)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits_decodes_BPSK_conv_dur = vitdec(bits_recus_BPSK_conv, treillis, 5*3, 'trunc', 'hard');

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DECODAGE VITERBI SOUPLE (distance euclidienne)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = 3; % nombre de niveaux (3 bits)
        dist_souple = round(((1 - Signal_echantillonne_BPSK_conv) / 2) * (2^L - 1));
        % Saturation dans l'intervalle [0, 2^L-1]
        dist_souple(dist_souple < 0) = 0;
        dist_souple(dist_souple > 2^L-1) = 2^L-1;
        bits_decodes_BPSK_conv_souple = vitdec(dist_souple, treillis, 5*3, 'trunc', 'soft', L);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CALCUL DU TAUX D'ERREUR BINAIRE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TEB_BPSK = TEB_BPSK + sum(bits ~= bits_recus_BPSK) / length(bits);
        TEB_BPSK_conv_dur = TEB_BPSK_conv_dur + sum(bits ~= bits_decodes_BPSK_conv_dur(1:N)) / length(bits);
        TEB_BPSK_conv_souple = TEB_BPSK_conv_souple + sum(bits ~= bits_decodes_BPSK_conv_souple(1:N)) / length(bits);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CUMUL DU NOMBRE D'ERREURS ET NOMBRE DE CUMUL REALISES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb_erreurs = nb_erreurs + sum(bits ~= bits_decodes_BPSK_conv_dur(1:N));
        nb_cumul = nb_cumul + 1;

    end  %fin boucle sur comptage nombre d'erreurs

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CALCUL DU TAUX D'ERREUR SYMBOLE ET DU TAUX D'ERREUR BINAIRE POUR LA
    %VALEUR TESTEE DE Eb/N0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TES_simule_BPSK(indice_bruit) = TES_BPSK / nb_cumul;
    TEB_simule_BPSK(indice_bruit) = TEB_BPSK / nb_cumul;
    TES_simule_BPSK_conv(indice_bruit) = TES_BPSK_conv / nb_cumul;

    TEB_simule_BPSK_conv_dur(indice_bruit) = TEB_BPSK_conv_dur / nb_cumul;
    TEB_simule_BPSK_conv_souple(indice_bruit) = TEB_BPSK_conv_souple / nb_cumul;

end  %fin boucle sur les valeurs testées de Eb/N0

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TES ET TEB THEORIQUES CHAINES IMPLANTEES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TES_THEO_BPSK = qfunc(sqrt(2 * tab_Eb_N0));
TEB_THEO_BPSK = TES_THEO_BPSK;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACES DES TES ET TEB OBTENUS EN FONCTION DE Eb/N0
%COMPARAISON AVEC LES TES et TEBs THEORIQUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(tab_Eb_N0_dB, TES_THEO_BPSK,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TES_simule_BPSK,'b-o')
hold on
semilogy(tab_Eb_N0_dB, TES_simule_BPSK_conv,'g-s')
legend('TES théorique BPSK','TES simulé BPSK', 'TES simulé BPSK avec codage convolutif')
xlabel('E_b/N_0 (dB)')
ylabel('TES')

figure
semilogy(tab_Eb_N0_dB, TEB_THEO_BPSK,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TEB_simule_BPSK,'b-o')
hold on
semilogy(tab_Eb_N0_dB, TEB_simule_BPSK_conv_dur,'g-s')
hold on
semilogy(tab_Eb_N0_dB, TEB_simule_BPSK_conv_souple,'k-d')
legend('TEB théorique BPSK','TEB simulé BPSK', 'TEB simulé BPSK avec codage convolutif en dur', 'TEB simulé BPSK avec codage convolutif en souple')
% legend('TEB théorique BPSK', 'TEB simulé BPSK avec codage convolutif en souple')
xlabel('E_b/N_0 (dB)')
ylabel('TEB')