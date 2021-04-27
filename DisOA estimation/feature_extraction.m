clear all
close all
clc

% Positions of source signal for training and testing for each distance

%train 50 points 0.5 meter
x=[3.996200909 3.98486137 3.966153701 3.940362191 3.907878778 3.86919709 3.82490495 3.775675437 3.722256659 3.665460389 3.606149722 3.545225966 3.48361494 3.422252907 3.362072347 3.303987786 3.248881897 3.19759209 3.150897784 3.109508561 3.074053388 3.045071054 3.023001986 3.008181552 3.00083497 3.001073881 3.008894654 3.024178442 3.046692987 3.07609615 3.11194111 3.153683153 3.200687953 3.252241207 3.307559493 3.365802176 3.426084177 3.48748943 3.549084798 3.609934256 3.669113113 3.725722065 3.778900863 3.827841383 3.871799906 3.910108423 3.942184783 3.967541543 3.985793372 3.996662909];
y=[3.438480429 3.37789573 3.319166576 3.263185429 3.210803003 3.16281532 3.11995161 3.082863267 3.052113879 3.02817073 3.01139767 3.002049589 3.000268542 3.006081596 3.019400414 3.040022596 3.067634763 3.101817308 3.142050782 3.187723782 3.238142245 3.292539994 3.350090379 3.409918846 3.471116219 3.53275252 3.593891102 3.653602881 3.710980454 3.765151892 3.815293986 3.860644758 3.900515041 3.934298953 3.9614831 3.981654382 3.99450627 3.999843461 3.99758485 3.987764758 3.970532417 3.946149694 3.914987118 3.877518248 3.834312474 3.786026365 3.733393697 3.677214293 3.618341876 3.557671094];
label = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% test 25 points 0.5 meter
% x=[3.985998123 3.944776703 3.878644449 3.791305263 3.687650796 3.573486475 3.455206361 3.339435026 3.232656535 3.140851285 3.069161059 3.021601048 3.00083497 3.008025881 3.042771035 3.103124444 3.185705863 3.285890113 3.398066125 3.515951198 3.632942885 3.742488772 3.838453467 3.915462228 3.969201985];
% y=[3.382501812 3.271584405 3.173459985 3.093624258 3.036548623 3.005429744 3.002010512 3.02648243 3.077474886 3.152131921 3.246272179 3.3546231 3.471116219 3.589227046 3.702340489 3.804121346 3.888869124 3.951837312 3.989499219 3.999745495 3.982002271 3.937263302 3.868034306 3.778192626 3.67277007];
% label = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


% train 50 points 1 meter
% x=[4.492401819 4.469722739 4.432307401 4.380724382 4.315757555 4.238394181 4.149809901 4.051350873 3.944513318 3.830920778 3.712299445 3.590451932 3.46722988 3.344505814 3.224144694 3.107975571 2.997763794 2.89518418 2.801795567 2.719017122 2.648106776 2.590142108 2.546003971 2.516363104 2.50166994 2.502147762 2.517789308 2.548356884 2.593385974 2.6521923 2.723882219 2.807366306 2.901375905 3.004482414 3.115118987 3.231604352 3.352168354 3.47497886 3.598169596 3.719868511 3.838226225 3.95144413 4.057801727 4.155682766 4.243599812 4.320216845 4.384369566 4.435083086 4.471586744 4.493325818];
% y=[3.376960858 3.255791464 3.138333151 3.026370859 2.921606007 2.82563064 2.739903234 2.665726535 2.604227758 2.556341461 2.522795341 2.504099178 2.500537085 2.512163193 2.538800827 2.580045193 2.635269526 2.703634616 2.784101564 2.875447565 2.976284491 3.085079987 3.200180759 3.319837692 3.442232438 3.56550504 3.687782204 3.807205761 3.921960908 4.030303785 4.130587972 4.221289516 4.301030083 4.368597906 4.4229662 4.463308765 4.48901254 4.499686922 4.495169699 4.475529517 4.441064833 4.392299388 4.329974237 4.255036496 4.168624947 4.072052731 3.966787394 3.854428586 3.736683753 3.615342188];
% label = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

% test 25 points 1 meters
% x=[4.471996247 4.389553407 4.257288898 4.082610527 3.875301592 3.646972951 3.410412721 3.178870052 2.96531307 2.78170257 2.638322118 2.543202097 2.50166994 2.516051761 2.58554207 2.706248888 2.871411726 3.071780227 3.296132249 3.531902396 3.765885769 3.984977543 4.176906935 4.330924456 4.43840397];
% y=[3.265003624 3.043168809 2.846919971 2.687248516 2.573097246 2.510859488 2.504021025 2.55296486 2.654949773 2.804263842 2.992544358 3.2092462 3.442232438 3.678454092 3.904680978 4.108242691 4.277738248 4.403674624 4.478998437 4.499490989 4.464004542 4.374526605 4.236068612 4.056385252 3.845540141];
% label = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

% train 50 points 1.5 meters
% x=[4.988602728 4.954584109 4.898461102 4.821086573 4.723636333 4.607591271 4.474714851 4.32702631 4.166769977 3.996381166 3.818449167 3.635677899 3.45084482 3.26675872 3.086217041 2.911963357 2.746645691 2.59277627 2.452693351 2.328525683 2.222160164 2.135213163 2.069005957 2.024544656 2.00250491 2.003221642 2.026683962 2.072535326 2.140078961 2.22828845 2.335823329 2.461049459 2.602063858 2.75672362 2.92267848 3.097406528 3.278252531 3.46246829 3.647254394 3.829802767 4.007339338 4.177166196 4.33670259 4.483524149 4.615399718 4.730325268 4.826554349 4.902624629 4.957380116 4.989988726];
% y=[3.315441287 3.133687196 2.957499727 2.789556288 2.63240901 2.48844596 2.359854852 2.248589802 2.156341637 2.084512191 2.034193011 2.006148766 2.000805627 2.018244789 2.058201241 2.120067789 2.202904289 2.305451925 2.426152347 2.563171347 2.714426736 2.877619981 3.050271138 3.229756539 3.413348657 3.59825756 3.781673306 3.960808642 4.132941363 4.295455677 4.445881958 4.581934274 4.701545124 4.802896858 4.8844493 4.944963147 4.98351881 4.999530383 4.992754549 4.963294275 4.91159725 4.838449081 4.744961355 4.632554744 4.502937421 4.358079096 4.200181091 4.031642879 3.855025629 3.673013282];
% label = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

%test 25 points 1.5 meters
% x=[4.95799437 4.83433011 4.635933347 4.37391579 4.062952388 3.720459426 3.365619082 3.018305077 2.697969605 2.422553855 2.207483178 2.064803145 2.00250491 2.024077642 2.128313105 2.309373331 2.557117589 2.85767034 3.194198374 3.547853594 3.898828654 4.227466315 4.515360402 4.746386684 4.907605955];
% y=[3.147505436 2.814753214 2.520379956 2.280872775 2.109645869 2.016289233 2.006031537 2.07944729 2.232424659 2.456395763 2.738816538 3.063869301 3.413348657 3.767681139 4.107021467 4.412364037 4.666607371 4.855511936 4.968497656 4.999236484 4.946006814 4.811789907 4.604102918 4.334577878 4.018310211];
% label =[2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

% train 50 points 2 meter
% x=[5.484803637 5.439445478 5.364614802 5.261448764 5.13151511 4.976788362 4.799619801 4.602701747 4.389026637 4.161841555 3.924598889 3.680903865 3.434459759 3.189011627 2.948289387 2.715951142 2.495527588 2.29036836 2.103591134 1.938034244 1.796213552 1.680284217 1.592007943 1.532726208 1.50333988 1.504295523 1.535578616 1.596713768 1.686771948 1.8043846 1.947764438 2.114732612 2.30275181 2.508964827 2.730237974 2.963208703 3.204336708 3.44995772 3.696339192 3.939737023 4.17645245 4.402888261 4.615603454 4.811365532 4.987199624 5.14043369 5.268739132 5.370166172 5.443173488 5.486651635];
% y=[3.253921716 3.011582928 2.776666302 2.552741718 2.343212014 2.15126128 1.979806469 1.83145307 1.708455516 1.612682921 1.545590682 1.508198355 1.501074169 1.524326385 1.577601654 1.660090386 1.770539051 1.907269233 2.068203129 2.25089513 2.452568981 2.670159974 2.900361518 3.139675385 3.384464876 3.63101008 3.875564408 4.114411522 4.343921817 4.560607569 4.761175944 4.942579032 5.102060165 5.237195811 5.345932399 5.426617529 5.47802508 5.499373845 5.490339399 5.451059033 5.382129667 5.284598775 5.159948473 5.010072992 4.837249894 4.644105462 4.433574787 4.208857172 3.973367506 3.730684375];
% label = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

% test 25 points 2 meters
% x=[5.44399 5.27911 5.01458 4.66522 4.25060 3.79395 3.32083 2.85774 2.43063 2.06341 1.77664 1.58640 1.50334 1.53210 1.67108 1.91250 2.24282 2.64356 3.09226 3.56380 4.03177 4.46996 4.85381 5.16185 5.37681];
% y=[3.03001 2.58634 2.19384 1.87450 1.64619 1.52172 1.50804 1.60593 1.80990 2.10853 2.48509 2.91849 3.38446 3.85691 4.30936 4.71649 5.05548 5.30735 5.45800 5.49898 5.42801 5.24905 4.97214 4.61277 4.19108];
% label=[3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];


% Reading signals from the LibriSpeech dataset
% all signals belong to one speech 
[sig1,fs1] = audioread('1.flac'); 
[sig2,fs2] = audioread('2.flac'); 
[sig3,fs3] = audioread('3.flac'); 
[sig4,fs4] = audioread('4.flac'); 
[sig5,fs5] = audioread('5.flac'); 
[sig6,fs6] = audioread('6.flac'); 
[sig7,fs7] = audioread('7.flac'); 
[sig8,fs8] = audioread('8.flac'); 
[sig9,fs9] = audioread('9.flac'); 
[sig10,fs10] = audioread('10.flac'); 
% concatenate to obtain longer duration of signal
sig = [sig1 ;sig2; sig3; sig4; sig5 ;sig6; sig7; sig8; sig9; sig10];

I = length(x); % length of x and y are equal
data_arr = zeros(I,6);
for i = 1:I
    
[h11,h12] = rir_gen_dist(fs1,x(i),y(i));  % RIR Generation for two mics

withoutNoiseS11 = filter(h11,1,sig); 
withoutNoiseS12 = filter(h12,1,sig);
withoutNoiseS1 = [withoutNoiseS11,withoutNoiseS12]; 

snr = 20; % Signal-to-Noise Ratio in dB
withNoiseS11 = awgn(withoutNoiseS11,snr); %noise addition
withNoiseS12 = awgn(withoutNoiseS12,snr);

withNoiseS1 = [withNoiseS11,withNoiseS12]; %Recorded signal for noisy room

    k = i; %for defining the filename

    myfilename1 = sprintf('train_SNR_0.5m_room1_S11_%g.flac', k); % define file name 
    myfilename2 = sprintf('train_SNR_0.5m_room1_S12_%g.flac', k); % define file name 
    
    %uncomment for quiet room
    %audiowrite(myfilename1,withoutNoiseS1(:,1),fs1);  
    %audiowrite(myfilename2,withoutNoiseS1(:,2),fs1);
    
    audiowrite(myfilename1,withNoiseS1(:,1),fs1);
    audiowrite(myfilename2,withNoiseS1(:,2),fs1);
    
    data_arr(i,1) = cdr_dereverc(myfilename1,myfilename2); % DIFF extraction
    
    [RT1, DRR1] = irStats(myfilename1);
    [RT2, DRR2] = irStats(myfilename2); 
    data_arr(i,2) = DRR1; % DRR extraction for mic 1
    data_arr(i,3) = DRR2; % DRR extraction for mic 2
    
    [x2,fs_in] = audioread(myfilename1);
    [y2,fs_in] = audioread(myfilename2);
    [a,fc]=mscohere(x2,y2,hann(512),128,1024); % MSC extraction
    data_arr(i,4) = mean(a);
    
    data_arr(i,5)= magspectra(myfilename1,myfilename2); % BSMD STD extraction
    
    data_arr(i,6) = label(1,i); % labelling
    
    
end
writematrix(data_arr,'train_SNR_0.5m_room1.csv') %writing data to a csv file