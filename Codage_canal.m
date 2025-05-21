%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETUDE D'UN CODAGE DE TYPE BLOC ET D'UN CODAGE CONVOLUTIF
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

col_1 = [1; 1; 1; 0];
col_2 = [0; 1; 1; 1];
col_3 = [1; 1; 0; 1];
G = [eye(4) col_1 col_2 col_3];

%codage bloc: {0,1}^4 -> {0,1}^7


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
    
    TES_BPSK_hamming = 0;
    TEB_BPSK_hamming = 0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BOUCLE POUR PRECISION TES ET TEBS MESURES :COMPTAGE NOMBRE ERREURS
    % (voir annexe texte TP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(nb_erreurs < 200)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %GENERATION DE L'INFORMATION BINAIRE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits = randi([0, 1], 1, N);
        bits_hamming = mod(reshape(bits, [N/4, 4]) * G, 2); % mot de codes
        bits_hamming = reshape(bits_hamming', length(bits_hamming), [])';
        bits_hamming_line = reshape(bits_hamming, numel(bits_hamming), [])';
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Mapping BPSK
        symboles_BPSK = pskmod(bits, M, 0, InputType="bit", PlotConstellation=false);
        a_k = real(symboles_BPSK);
        b_k = zeros(size(a_k));

        %Mapping BPSK avec Hamming
        symboles_BPSK_hamming = pskmod(bits_hamming_line, M, 0, InputType="bit", PlotConstellation=false);
        a_k_hamming = real(symboles_BPSK_hamming);
        b_k_hamming = zeros(size(a_k));
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SURECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        somme_Diracs_ponderes_BPSK = kron(a_k, [1 zeros(1, Ns-1)]);
        somme_Diracs_ponderes_BPSK_hamming = kron(a_k_hamming, [1 zeros(1, Ns-1)]);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %FILTRAGE DE MISE EN FORME 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Génération de la réponse impulsionnelle du filtre de mise en forme
        h = ones(1, Ns);

        %Filtrage de mise en forme
        Signal_emis_BPSK = filter(h, 1, somme_Diracs_ponderes_BPSK);
        Signal_emis_BPSK_hamming = filter(h, 1, somme_Diracs_ponderes_BPSK_hamming);
 
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CANAL DE PROPAGATION AWGN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %POUR MODULATION BPSK
        %Calcul de la puissance du signal émis en BPSK
        P_signal = (1 / (length(Signal_emis_BPSK))) * sum(norm(Signal_emis_BPSK).^2);

        %Calcul de la puissance du bruit à ajouter au signal pour obtenir la valeur
        %souhaité pour le SNR par bit à l'entrée du récepteur (Eb/N0)
        P_bruit = (P_signal * Ns) / (2 * log2(M) * Eb_N0);

        %Génération du bruit gaussien à la bonne puissance en utilisant la fonction
        %randn de Matlab (bruit réel ici car en BPSK l'enveloppe complexe est réelle)
        Bruit =sqrt(P_bruit) * randn(1, length(Signal_emis_BPSK));

        %Ajout du bruit canal au signal émis => signal à l'entrée du récepteur
        % Signal_recu_BPSK = Signal_emis_BPSK;
        Signal_recu_BPSK = Signal_emis_BPSK + Bruit;


        %POUR MODULATION BPSK avec Hamming
        P_signal_hamming = (1 / (length(Signal_emis_BPSK_hamming))) * sum(norm(Signal_emis_BPSK_hamming).^2);

        P_bruit_hamming = (P_signal_hamming * Ns) / (2 * log2(M) * Eb_N0);

        Bruit_hamming =sqrt(P_bruit_hamming) * randn(1, length(Signal_emis_BPSK_hamming));

        Signal_recu_BPSK_hamming = Signal_emis_BPSK_hamming + Bruit_hamming;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %FILTRAGE DE RECEPTION ADAPTE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Réponse impulsionnelle du filtre de réception
        hr = ones(1, Ns);

        %Filtrage de réception
        Signal_recu_filtre_BPSK = filter(hr, 1, Signal_recu_BPSK);
        Signal_recu_filtre_BPSK_hamming = filter(hr, 1, Signal_recu_BPSK_hamming);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ECHANTILLONNAGE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Choix de n0 
        n0 = 4;     % A COMPLETER
        %Echantillonnage à n0+mNs
        Signal_echantillonne_BPSK = Signal_recu_filtre_BPSK(n0:Ns:end);
        Signal_echantillonne_BPSK_hamming = Signal_recu_filtre_BPSK_hamming(n0: Ns: end);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DECISIONS SUR LES SYMBOLES 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        symboles_recus_decimaux_BPSK = real(Signal_echantillonne_BPSK) > 0;
        symboles_recus_decimaux_BPSK = symboles_recus_decimaux_BPSK * 2 - 1;
        
        symboles_recus_decimaux_BPSK_hamming = real(Signal_echantillonne_BPSK_hamming) > 0;
        symboles_recus_decimaux_BPSK_hamming = symboles_recus_decimaux_BPSK_hamming * 2 - 1;
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CALCUL DU TAUX D'ERREUR SYMBOLE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TES_BPSK = TES_BPSK + sum(real(symboles_BPSK) ~= symboles_recus_decimaux_BPSK) / length(real(symboles_BPSK));
        TES_BPSK_hamming = TES_BPSK_hamming + sum(real(symboles_BPSK_hamming) ~= symboles_recus_decimaux_BPSK_hamming) / length(real(symboles_BPSK_hamming));

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DEMAPPING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bits_recus_BPSK = pskdemod(symboles_recus_decimaux_BPSK, M);
        bits_recus_BPSK_hamming = pskdemod(symboles_recus_decimaux_BPSK_hamming, M);
        bits_recus_BPSK_hamming = reshape(bits_recus_BPSK_hamming, [length(bits_recus_BPSK_hamming)/7, 7]);

        %cas DUR - A REFAIRE CAR IMPLEMENTATION NON CORRECTE
        % likelihood_matrix = int2bit(bits_recus_BPSK_hamming * bits_hamming, 7);
        % likelihood_matrix = reshape(likelihood_matrix', [250 250 7]);
        % likelihood_matrix = sum(likelihood_matrix, 3);
        % [mins, min_matrix_id] = min(likelihood_matrix, [], 2);
        % bits_decides_BPSK_hamming_dur = bits_hamming(min_matrix_id);

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CALCUL DU TAUX D'ERREUR BINAIRE CUMULE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TEB_BPSK = TEB_BPSK + sum(bits ~= bits_recus_BPSK) / length(bits);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CUMUL DU NOMBRE D'ERREURS ET NOMBRE DE CUMUL REALISES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb_erreurs=nb_erreurs + sum(bits ~= bits_recus_BPSK);
        nb_cumul = nb_cumul + 1;

    end  %fin boucle sur comptage nombre d'erreurs

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CALCUL DU TAUX D'ERREUR SYMBOLE ET DU TAUX D'ERREUR BINAIRE POUR LA
    %VALEUR TESTEE DE Eb/N0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TES_simule_BPSK(indice_bruit) = TES_BPSK / nb_cumul;
    TEB_simule_BPSK(indice_bruit) = TEB_BPSK / nb_cumul;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DIAGRAMME DE L'OEIL EN SORTIE DU FILTRE DE RECEPTION AVEC BRUIT
    %TRACE POUR CHAQUE VALEUR DE Eb/N0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MODULATION BPSK
    % oeil=reshape(Signal_recu_BPSK, Ns, length(Signal_recu_filtre_BPSK) / Ns);
    % figure
    % plot(oeil)
    % title(['Tracé du diagramme de l"oeil en sortie du filtre de réception (BPSK) pour E_b/N_0 = ' num2str(Eb_N0_dB) 'dB'])


end  %fin boucle sur les valeurs testées de Eb/N0

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TRACE DES CONSTELLATIONS APRES ECHANTILLONNAGE POUR CHAQUE VALEUR DE Eb/N0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MODULATION BPSK
    figure
    plot(a_k, b_k, 'b*');
    xlabel('a_k')
    ylabel('b_k')
    title(['Tracé de la constellation en sortie du filtre de réception (BPSK) pour E_b/N_0 = ' num2str(Eb_N0_dB) 'dB'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TES ET TEB THEORIQUES CHAINES IMPLANTEES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TES_THEO_BPSK = qfunc(sqrt(2 * tab_Eb_N0));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCUL DU TES ET TEB THEORIQUE DE LA CHAINE IMPLANTEE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEB_THEO_BPSK = TES_THEO_BPSK;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRACES DES TES ET TEB OBTENUS EN FONCTION DE Eb/N0
%COMPARAISON AVEC LES TES et TEBs THEORIQUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(tab_Eb_N0_dB, TES_THEO_BPSK,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TES_simule_BPSK,'b-o')
legend('TES théorique BPSK','TES simulé BPSK')
xlabel('E_b/N_0 (dB)')
ylabel('TES')

figure
semilogy(tab_Eb_N0_dB, TEB_THEO_BPSK,'r-x')
hold on
semilogy(tab_Eb_N0_dB, TEB_simule_BPSK,'b-o')
legend('TEB théorique BPSK','TEB simulé BPSK')
xlabel('E_b/N_0 (dB)')
ylabel('TEB')