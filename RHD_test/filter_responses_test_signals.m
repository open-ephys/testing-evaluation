% create test signals for evaluating a RHD based open-ephys headstage
% assumes a nidaq Aout card
%
%
%
%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%

%% config
samplerate=50000;

freq_min=-1; % from 0.1
freq_max=4;  % to 10000Hz
freq_N=50;   % N freq. points

secs_sweep=10; % s per sweep
triangle_N=10;

freq_points= logspace(freq_min,freq_max,freq_N);

amp_points=[1 3 5]; % 1-5mV

%% set up nidaq
s = daq.createSession('ni');

s.Rate=samplerate;
s.addAnalogOutputChannel('Dev1',0,'Voltage') & output test signal

%% run test
%start at 0
train=zeros(samplerate*1,1);
s.queueOutputData(train);
s.startForeground();
pause(5); % 5 sec delay after each

for amp=amp_points
    
    for freq=freq_points
        train=sin( linspace(0, pi*freq*2*secs_sweep ,samplerate*secs_sweep) )*amp;
        train(end)=0;
        s.queueOutputData(train');
        pause(5);
        fprintf('start %f\n',freq);
        plot(train);
        drawnow;
        s.startBackground();
        pause(secs_sweep+5); % 5 sec delay after each
    end;
    
    % impulse responses
    for i=1:10
        train=zeros(samplerate*1,1);
        
        train(1:round(samplerate/1000))=amp;
        s.queueOutputData(train);
        s.startForeground();
        pause(5); % 5 sec delay after each
    end;
    
    triangle_N=1;
    triangle_freq=10;
    % 1Hz triangle wave
    train=abs(mod( linspace(0,triangle_N*triangle_freq,samplerate*triangle_N)+triangle_N/4 ,1)-.5)*2*amp;
    
    s.queueOutputData(train');
    s.startForeground();
    
    pause(5); % 5 sec delay after each
    
    square_N=10;
    square_freq=.5;
    % 1Hz square wave
    train=(mod( linspace(0,square_N*square_freq,samplerate*square_N) ,1)>.5)*amp;
    s.queueOutputData(train');
    s.startForeground();
    
    pause(5); % 5 sec delay after each
    
end;


%% test input range with triangle wave
triangle_N=1;
triangle_freq=2;
amp=6;
% 1Hz triangle wave
train=abs(mod( linspace(0,triangle_N*triangle_freq,samplerate*triangle_N)+triangle_N/2 ,1)-.5)*2*amp;
train=[train, -train];
train(end)=0;
plot(train);
s.queueOutputData(train');
s.startForeground();
