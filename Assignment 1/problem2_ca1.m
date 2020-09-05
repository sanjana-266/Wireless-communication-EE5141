% Assignment 1 EE5141 %
% Simulation of Clarkes model %
% Note: All the plot generating codes have been commented out, please
% uncomment them one by one to view the plots
% Defining the parameters- maximum doppler frequency, frequency spacing%

fd=10;
N=1024;
df=2*fd/(N-1);
S=zeros(N/2,1);
time_duration = 1/df;
noise_i_0=zeros([N/2,1]);% in-phase component of noise
noise_q_0=zeros([N/2,1]);% quadrature component of noise

for iterations=1:50 % generating 50 sample functions
    
    for i=1:N/2
        S(i)=1.5/(pi*fd*sqrt(1-(i*df/fd)^2));% generating the spectral filter
        noise_i_0(i)=randn+1j*randn;% generating in-phase component of noise
        noise_q_0(i)=randn+1j*randn;% generating quadrature component of noise
    end
    
    % for negative freq, the components are complex conjugate( conjugate symmetry )
    S1 = [conj(flip(S)); S];
    noise_i = [conj(flip(noise_i_0)); noise_i_0];
    noise_q = [conj(flip(noise_q_0)); noise_q_0];
    
    output_i= noise_i.*sqrt(S1);
    time_signal_i = ifft(output_i,N);
    t1=abs(time_signal_i).*abs(time_signal_i);% computing the in-phase time series
    
    output_q= noise_q.*sqrt(S1);
    time_signal_q = ifft(output_q,N);
    t2=abs(time_signal_q).*abs(time_signal_q);% computing the quadrature time series
    
    t(:,iterations) = sqrt(t1+t2);
    normalised_t(:,iterations) = t(:,iterations)/rms(t(:,iterations));% generating the Rayleigh fading signal
    
    threshold = [0.01,0.1,1];% for calculating the level crossing rate
    
    for i=1:length(threshold)
        N_v(i)=0;% counter to store the number of times the envelope crosses the threshold
        time(i)=0;% counter to store the total time the envelope is below the threshold
        for k=2:N
            if normalised_t(k,iterations) >threshold(i) && normalised_t(k-1,iterations)< threshold(i)
                N_v(i)=N_v(i)+1;
            end 
        end
        
        for k=1:N
            if normalised_t(k,iterations) < threshold(i)
                time(i)=time(i)+1;
            end            
        end
        level_crossing_rate(i,iterations)=N_v(i);
        tau(i,iterations)=time(i)/N;
    end
end

level_crossing_rate_avg=sum(level_crossing_rate,2)./50;% taking an average over the 50 iterations
tau_avg=sum(tau,2)./50;

rayleigh_dB(:) = 20*log10(normalised_t(:,1));
rayleigh_dB(:) = rayleigh_dB(:)./rms(rayleigh_dB);
time = [1:1:1024]*time_duration/1024;

% Plot 1- normalised Rayleigh fading waveform envelope versus time 

% plot(time,rayleigh_dB);
% xlabel 'Time(in seconds)'
% ylabel 'Received envelope v (in dB)'
% title 'Clarkes model- f_d = 10Hz'

% calculation of the theoretical parameters using the formula
   
level_crossing_rate_theoretical(1:3)=time_duration*sqrt(2*pi)*fd.*threshold(:).*exp(-threshold(:).*threshold(:));
tau_theoretical(1:3)=(exp(threshold(:).*threshold(:))-1)./(threshold(:).*fd.*sqrt(2*pi));

% Plot 2- comparing Theoretical and Measured LCR

% stem(threshold,level_crossing_rate_avg,'r','LineWidth',2);
% hold on
% stem(threshold,level_crossing_rate_theoretical,'b','LineWidth',2);
% xlabel 'Threshold of normalised received envelope'
% ylabel 'Level crossing rate(N_v)'
% title 'Clarkes model- f_d = 100Hz'
% legend ('Measured LCR','Theoretical LCR');

% Plot 3- comparing Theoretical and Measured duration of fade

% stem(threshold,tau_avg,'r','LineWidth',2);
% hold on
% stem(threshold,tau_theoretical,'b','LineWidth',2);
% xlabel 'Threshold of normalised received envelope'
% ylabel 'Average duration of fades(tau)'
% title 'Clarkes model- f_d = 10Hz'
% legend ('Measured ','Theoretical');