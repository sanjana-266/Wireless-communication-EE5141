% Assignment 1 EE5141 %
% Simulation of Jakes model %
% Note: All the plot generating codes have been commented out, please
% uncomment them one by one to view the plots
% Defining the parameters- maximum doppler frequency, number of oscillators %

fd=1;
M=20;
N=4*M+2;

time_duration=1;
no_of_samples=1000;% samples per second
P=time_duration*no_of_samples;

z_r=zeros(P,1);% real part of the fading waveform
z_i=zeros(P,1);% imaginary part of the fading waveform

beta=zeros(M,1);
phi=zeros(M,1);

for i=1:M % defining all parameters of the M oscillators
    beta(i)=pi*i/(M+1);
    phi(i)=pi*(-1+2*rand);
    f(i)=fd*cos(2*pi*i/N);
end
 
for i=1:P % computation of the fading waveform for each time instant
    sum1=0;
    sum2=0;
    for j=1:M
      sum1=sum1 + cos(beta(j))*cos((2*pi*f(j)*i/1000)+phi(j)); 
      sum2=sum2 + sin(beta(j))*cos((2*pi*f(j)*i/1000)+phi(j));
    end
    z_r(i)=2*(2*sum1 + sqrt(2)*cos((2*pi*fd*i/1000)+phi(1)))/sqrt(N);
    z_i(i)=2*(2*sum2)/sqrt(N);
end
z = z_r + 1j*z_i;
m = [1:1:P];

% Plot 1- fading waveform envelope versus time 

% plot(m/1000,20*log10(abs(z)),'r','LineWidth',1);
% xlabel 'Time(in seconds)'
% ylabel 'Received envelope v (in dB)'
% title 'Jakes model- f_d = 1Hz'

rho = z/rms(z);% normalised received envelope

% Plot 2- normalised fading waveform envelope versus time 

% plot(m/1000,20*log10(abs(rho)),'r','LineWidth',2);
% xlabel 'Time(in seconds)'
% ylabel 'Normalised received envelope v (in dB)'
% title 'Jakes model- f_d = 100Hz'

threshold = [-22: 3: 5];% for calculating the level crossing rate

for i=1:length(threshold)
    N_v(i)=0;% counter to store the number of times the envelope crosses the threshold
    time(i)=0;% counter to store the total time the envelope is below the threshold
    for k=1:P-1
         if 20*log10(abs(rho(k+1))) >=threshold(i) && 20*log10(abs(rho(k)))<= threshold(i)
             N_v(i)=N_v(i)+1;
         end
    end
    
    for k=1:P     
         if 20*log10(abs(rho(k))) < threshold(i)
            time(i)=time(i)+1;
         end   
    end
    
   level_crossing_rate(i)=N_v(i)/fd;% normalised LCR
   tau(i)=time(i)*fd/P;% divided by P(total number of samples) to find average fading duration
   rho1(i)=10^(threshold(i)/20);% rho as a number
   % calculation of the theoretical parameters using the formula
   level_crossing_rate_theoretical(i)=time_duration*sqrt(2*pi)*rho1(i)*exp(-rho1(i)*rho1(i));
   tau_theoretical(i)=(exp(rho1(i)*rho1(i))-1)/(rho1(i)*sqrt(2*pi));
end

r=[1:1:length(threshold)];
% Plot 3- comparing Theoretical and Measured LCR

% stem(threshold,level_crossing_rate,'r','LineWidth',2);
% hold on
% stem(threshold,level_crossing_rate_theoretical,'b','LineWidth',2);
% xlabel 'Threshold of normalised received envelope(in dB)'
% ylabel 'Normalised level crossing rate(N_v/fd)'
% title 'Jakes model- f_d = 100Hz'
% legend ('Measured LCR','Theoretical LCR');

% Plot 4- comparing Theoretical and Measured duration of fade

% stem(threshold,tau,'r','LineWidth',2);
% hold on
% stem(threshold,tau_theoretical,'b','LineWidth',2);
% xlabel 'Threshold of normalised received envelope(in dB)'
% ylabel 'Normalized duration of fades(tau*fd)'
% title 'Jakes model- f_d = 10Hz'
% legend ('Measured ','Theoretical');

% modifying Jakes model to obtain Ricean fading

sigma = sqrt(var(z));% standard dev of the Jakes model fading
k = 10;
fd = 100;
phi=zeros(1000,1);
phi(:)=pi*(-1+2*rand);% generating a random phi

for i=1:1000 % generating the Ricean and Rayleigh components

    ricean(i) = (sqrt(k/(k+1))*sigma*exp(1j*phi(i)))+sqrt(1/(k+1))*(sqrt(sigma/2)*(randn+1j*randn));
    rayleigh(i) = sqrt(1/(k+1))*(sqrt(sigma/2)*(randn+1j*randn));
    
end

% Plot 5- Ricean and Rayleigh waveforms

% subplot(2,1,1);
% plot([1:1:1000]/1000,20*log10(ricean));
% ylabel 'Ricean Fading envelope(in dB)'
% xlabel 'Time(in s)'
% title 'Modified Jakes model to generate Ricean Fading(k=10)- f_d = 100Hz'
% subplot(2,1,2);
% plot([1:1:1000]/1000,20*log10(rayleigh));
% ylabel 'Rayleigh Fading envelope(in dB)'
% xlabel 'Time(in s)'
% title 'Modified Jakes model to generate Rayleigh Fading- f_d = 100Hz'
