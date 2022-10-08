bands = [[0 170]; [170 310]; [310 600]; [600 1000]; [1000 3000]; [3000 6000]; [6000 12000]; [12000 14000]; [14000 16000]];

while true
  file = uigetfile('*.*'); 
  [y,oldFs] = audioread(file);

  duration = length(y)/oldFs;
  torg = linspace(0,duration,length(y));

  gains = zeros(1,length(bands));
  for i=1:length(bands)
    fprintf('Enter gain (dB) of frequency band (');
    fprintf('%d',bands(i,1));
    fprintf('-');
    fprintf('%d): ',bands(i,2));
    gains(i) = input('');
  end

  choice = 0;
  while choice ~= 1 && choice ~= 2
    choice = input('Enter type of filter (1.FIR, 2.IIR): ');
  end

  Fs = 0;
  while Fs <= 0
    Fs = input('Enter output sample rate: ');
  end
  
  highest = bands(length(bands),2) + (bands(length(bands),2) - bands(length(bands),1))/8;
  if oldFs < 2 * highest
    ratio = (2*highest)/oldFs;
    aa  = log2(ratio);
    pow = ceil(log2(ratio));
    workFs = oldFs * (2^pow);
    Y = resample(y,workFs,oldFs);
    twork = linspace(0,duration,length(Y));
  else
    workFs = oldFs;
    Y = y;
    twork = torg;
  end
  
  output = zeros(length(Y),1);
  Fn = workFs/2;
  for i=1:length(bands)
    Rp = 1;  % dB
    Rs = 80; % dB
    fac = (bands(i,2) - bands(i,1)) / 10;
    if choice == 1 % FIR
      if i==1 % LOWPASS
        Wc = [bands(1,2) (bands(1,2)+fac)];
        mags = [1 0];
        devs = [(10^(-Rp/20)) 10^(-Rs/20)];
      else % BANDPASS
        if bands(i,1)-fac<0
          Wc = [0 bands(i,1) bands(i,2) (bands(i,2)+fac)];
        else
          Wc = [(bands(i,1)-fac) bands(i,1) bands(i,2) (bands(i,2)+fac)];
        end
        mags = [0 1 0];
        devs = [10^(-Rs/20) (10^(-Rp/20)) 10^(-Rs/20)];
      end
      width = bands(i,2) - bands(i,1);
      [N,Wc,beta,ftype] = kaiserord(Wc,mags,devs,workFs);
      N = ceil(800*log10(N/130));
      num = fir1(N,Wc,ftype,kaiser(N+1,beta),'noscale');
      den = 1;
      sys = tf(num,1);
      [z,p,k] = zpkdata(sys);
    else % IIR
      if i==1 % LOWPASS
        Wp = bands(1,2)/Fn;
        Ws = (bands(1,2)+fac)/Fn;
        ftype = 'low';
      else % BANDPASS
        Wp = [bands(i,1) bands(i,2)]/Fn;
        if bands(i,1)-fac<0
          Ws = [0 (bands(i,2)+fac)]/Fn;
        else
          Ws = [(bands(i,1)-fac) (bands(i,2)+fac)]/Fn;
        end
        ftype = 'bandpass';
      end
      [N,Wn] = buttord(Wp,Ws,Rp,Rs);
      [z,p,k] = butter(N,Wn,ftype);
      [fsos,g] = zp2sos(z,p,k);
      sos = zp2sos(z,p,k);
    end
    
    % GAIN
    fprintf('Gain of frequency band (');
    fprintf('%d',bands(i,1));
    fprintf('-');
    fprintf('%d): %d\n',bands(i,2),k);
    
    figure('Name',sprintf('Band (%d-%d)', bands(i,1),bands(i,2)));
    % IMPULSE RESPONSE
    subplot(2,3,1);
    if choice == 1 % FIR
        impz(num,1);
    else % IIR
        impz(sos);
    end
    title(sprintf('Impulse Response (%d %d)', bands(i,1),bands(i,2)));
  
    % STEP RESPONSE
    subplot(2,3,4);
    if choice == 1 % FIR
        stepz(num,1);
    else % IIR
        stepz(sos);
    end
    title(sprintf('Step Response (%d %d)', bands(i,1),bands(i,2)));
    
    % ORDER 
    fprintf('Order of filter: %d.\n',N);
    
    % POLES/ZEROS
    subplot(2,3,3);
    if choice == 1 % FIR
      pzmap(sys);
    else % IIR
      zplane(z,p);
    end
    title(sprintf('Poles/Zeros (%d %d)', bands(i,1),bands(i,2)));
    
    % APPLY FILTER
    if choice == 1 % FIR
      temp = filtfilt(num,1,Y);
    else % IIR
      temp = sosfilt(sos,Y);
    end
    
    % DRAW TIME
    subplot(2,3,6);
    plot(twork,temp);
    title(sprintf('Time domain (%d %d)', bands(i,1),bands(i,2)));
    
    % DRAW FREQUENCY
    tempF = fftshift(fft(temp));
    Fvec = linspace(-workFs/2, workFs/2, length(tempF));
    Fmag = abs(tempF);
    subplot(2,3,2);
    plot(Fvec, Fmag);
    title(sprintf('Magnitude (frequency) (%d %d)', bands(i,1),bands(i,2)));
    Fphase = angle(tempF);
    subplot(2,3,5);
    plot(Fvec, Fphase);
    title(sprintf('Phase (frequency) (%d %d)', bands(i,1),bands(i,2)));
    
    a = axes;
    title(sprintf('Band (%d-%d)', bands(i,1),bands(i,2)));
    set(a,'Visible','off');
    set(a,'Position',[0,0,1,0.94]);
    
    % PHASE
    figure('Name',sprintf('Band (%d-%d) Phase', bands(i,1),bands(i,2)));
    if choice == 1 % FIR
        phasez(num,1,2^16,workFs);
    else % IIR
        phasez(fsos,2^16,workFs);
    end

    % ADD BAND TO OUTPUT
    Amp = 10^(gains(i)/20);
    temp = Amp * temp;
    output = output + temp;
  end

  % RESAMPLE TO OUTPUT FS
  output = resample(output,Fs,workFs);
  tnew = linspace(0,duration,length(output));

  % COMPARE TIME
  figure('Name','Compare Time and Frequency','NumberTitle','off');
  subplot(3,1,1);
  plot(torg,y);
  hold on;
  plot(tnew,output);
  title('Compare (time)');
  hold off;
  
  % COMPARE FREQUENCY
  yF = fftshift(fft(y));
  outputF = fftshift(fft(output));
  
  Fvec1 = linspace(-oldFs/2, oldFs/2, length(yF));
  Fvec2 = linspace(-Fs/2, Fs/2, length(outputF));
  
  Fmag1 = abs(yF);
  Fmag2 = abs(outputF);
  subplot(3,1,2);
  plot(Fvec1,Fmag1);
  hold on;
  plot(Fvec2,Fmag2);
  title('Compare Magnitude (frequency)');
  hold off;
  
  Fphase1 = angle(yF);
  Fphase2 = angle(outputF);
  subplot(3,1,3);
  plot(Fvec1,Fphase1);
  hold on;
  plot(Fvec2,Fphase2);
  title('Compare Phase (frequency)');
  hold off;
  
  a = axes;
  title('Compare Original vs Output signal');
  set(a,'Visible','off');
  set(a,'Position',[0,0,1,0.94]);
  
  % PLAY
  sound(output,Fs);
  
  % WRITE FILE
  audiowrite('output.wav',output,Fs);
  
  choice = input('Enter "Y"/"y" to continue, otherwise to exit: ', 's');
  if choice ~= 'Y' && choice ~='y'
    break
  end
end