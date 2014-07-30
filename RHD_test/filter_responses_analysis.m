% This code requires manual configuration of the paths and the boundaries
% between different test blocks
%
% quantifies freq. response by measuring abs peak amp per frequency block
%
% measures phase distortion by measuring the deviation between the hilbert
% transform and a straight line.
% However, the analysis is frequency-agnostic, so interference patterns
% above nyquist could result in very low reported phase distortion if the
% resulting interference pattern is neatly sinusodial
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


base_dir='/media/data_ephys3/oe_test'

exps(1).path='2014-03-28_14-01-06_0_10k_filter';    % 0.1Hz low cut
exps(1).filter=[0.1 10000];

exps(2).path='2014-03-28_15-10-09_1_75k_filter';    % 1Hz low cut
exps(2).filter=[1 75000];

exps(3).path='2014-03-28_16-09-07_324_5982_filter'; % spike band filter
exps(3).filter=[300 6000];

exps(4).path='2014-07-17_15-58-57_100uV';           % filter wide open but low stim. amp
exps(4).filter=[0.1 10000];



chs=[1:10 17]; % all channels that were used

use_channel=3;

freq_min=-1; % from 0.1
freq_max=4;  % to 10000Hz
freq_N=50;   % N freq. points

secs_sweep=10; % s per sweep
triangle_N=10;

freq_points= logspace(freq_min,freq_max,freq_N);

set(0,'DefaultFigureWindowStyle','docked');

%%

for exp_num=1:4
    
    fname=fullfile(base_dir,exps(exp_num).path);
    for ch=use_channel
        [data, timestamps, info] =load_open_ephys_data(fullfile(fname,['100_CH',num2str(ch),'.continuous']));
    end;
    
    data=data.*info.header.bitVolts;
    data_hilb=angle(hilbert(data));
    data_phase=angle(data_hilb);
    %% define block onsets <--- change this manually
    clear blocks;
    switch exp_num
        case 1
            blocks.onsets(1:50)=linspace(13.5,1001.5,50);
            blocks.amp(1:50)=1000;%uV
            blocks.ampnum(1:50)=1;
            
            blocks.onsets(51:100)=linspace(13.5,1001.5,50)+1091;
            blocks.amp(51:100)=3000;%uV
            blocks.ampnum(51:100)=2;
            
            blocks.onsets(101:150)=linspace(13.5,1001.5,50)+2182;
            blocks.amp(101:150)=5000;%uV
            blocks.ampnum(101:150)=3;
            
        case 2
            blocks.onsets(1:50)=linspace(13.5,1001.5,50)+2.9;
            blocks.amp(1:50)=1000;%uV
            blocks.ampnum(1:50)=1;
            
            blocks.onsets(51:100)=linspace(13.5,1001.5,50)+1093;
            blocks.amp(51:100)=3000;%uV
            blocks.ampnum(51:100)=2;
            
            blocks.onsets(101:150)=linspace(13.5,1001.5,50)+2184;
            blocks.amp(101:150)=5000;%uV
            blocks.ampnum(101:150)=3;
            
        case 3
            blocks.onsets(1:50)=linspace(13.5,1001.5,50)+0;
            blocks.amp(1:50)=1000;%uV
            blocks.ampnum(1:50)=1;
            
            blocks.onsets(51:100)=linspace(13.5,1001.5,50)+1090.5;
            blocks.amp(51:100)=3000;%uV
            blocks.ampnum(51:100)=2;
            
            blocks.onsets(101:150)=linspace(13.5,1001.5,50)+2181.5;
            blocks.amp(101:150)=5000;%uV
            blocks.ampnum(101:150)=3;
        case 4                          % only ran one experiment at low amplitude here
            blocks.onsets(1:50)=linspace(23,766,50)+0;
            blocks.amp(1:50)=100; %uV
            blocks.ampnum(1:50)=1;
    end;
    
    blocks.freq=[freq_points,freq_points,freq_points];
    
    
    %quantify peak attenuation
    for i=1:numel(blocks.onsets)
        t=blocks.onsets(i)+[4 10]'; % cut out first few seconds to let lc settle
        from=find_halfspace_iterative(timestamps,t(1));
        to=find_halfspace_iterative(timestamps,t(2));
        blocks.peak_attenuation(i)= quantile(abs(data(from:to)),.99);
    end;
    
    % plot block boundaries
    subs=5;
    
    figure(1+exp_num*10);
    clf; hold on;
    plot(timestamps(1:subs:end),data(1:subs:end));
    
    %  figure(2+exp_num*10);
    %  clf; hold on;
    % plot(timestamps(1:subs:end),data_phase(1:subs:end).*100,'r');
    
    plotvert(blocks.onsets,'k',.5);
    for i=1:numel(blocks.onsets)
        t=blocks.onsets(i);
        plot(t+[3 10],[0 0],'r','LineWidth',2 );
        plot(t+[3 10],[1 1].*  blocks.peak_attenuation(i),'r','LineWidth',2 );
    end;
    
    %% quantify phase distortion
    %
    % just fit a straight line to the hilbert transform and look at max.
    % deviation normalized to pi, this should give a nice estimate of max. phase
    % distortion in percent
    %
    phase_quant_points=100; % N points over which to quantify phase distortion
    for i=1:numel(blocks.onsets)
        t=blocks.onsets(i)+[0 10]'; % cut out first few seconds to let lc settle
        from=find_halfspace_iterative(timestamps,t(1));
        to=find_halfspace_iterative(timestamps,t(2));
        
        datablock=data(from:to);
        %smooth HF stuff to clean up phase in relevant freq. range
        freq=blocks.freq(i);
        sd=ceil( 100/freq );
        f=normpdf([-sd*4:sd*4],0,sd); f=f./sum(f);
        datablock_sm=conv(datablock,f,'same');
        
        datablock_hilb=angle(hilbert(datablock_sm));
        
        crossings= find( [0; diff(sign(datablock_hilb))] .* (abs(datablock_hilb)>1) ) ;
        crossings(crossings<50000)=[]; % remove early crossings so that we only get clean waveforms
        ncrosses=min(numel(crossings),20);
        if ncrosses>1
            e=zeros(ncrosses,phase_quant_points);
            for j=1:ncrosses-1
                ii=crossings(j):crossings(j+1)-1;
                %x=datablock_hilb(ii)-linspace(-pi,pi,numel(ii))'; % raw phase error
                
                if numel(ii)<5;
                    % just straight line
                    x=datablock_hilb(ii)-linspace(datablock_hilb(ii(1)),datablock_hilb(ii(end)),numel(ii))'; % raw phase error
                else
                    % do linear fit, accounts for crossing detection jitter -
                    % should be mroe accurate representation of true phase
                    % distortion
                    r=[ones(size(ii));ii];
                    b=regress(datablock_hilb(ii),r');
                    x=datablock_hilb(ii)-(b'*r)';
                    plot(ii,b'*r,'r');
                end;
                
                if numel(x)>2
                    e(j,:)=interp1(x,linspace(1,numel(x),phase_quant_points)); % interpolate to same number of points
                else
                    e(j,:)=zeros(1,phase_quant_points);
                end;
            end;
            if size(e,1)>1
                phase_error=mean(e);
            else
                phase_error=e;
            end;
        else
            phase_error=zeros(1,phase_quant_points);
        end;
        
        phase_error=(phase_error/pi)*100; % in percent of maximum phase error
        
        clf; hold on
        
        plot(crossings,0,'ks');
        plot(datablock.*.01,'g');
        plot(datablock_sm.*.01,'r');
        plot(datablock_hilb,'k');
        plot(linspace(0,100000,phase_quant_points),phase_error.*1,'LineWidth',2);
        
        drawnow;
        % pause(1);
        blocks.phase_error(i,:)=phase_error;
        
    end;
    
    
    %% plot freq. response and max. phase distortion
    figure(3+exp_num*10);
    clf; subplot(2,1,1);
    for i=1:size(blocks.onsets,1)
        c=[1-i/4,0,i/3]';
        signamp=mean(blocks.amp(blocks.ampnum==i));
        semilogx(freq_points, blocks.peak_attenuation(blocks.ampnum==i)./signamp,'o-','color',c);
        %plot(log(freq_points), blocks.peak_attenuation(blocks.ampnum==i)./signamp,'o-','color',c);
        text(10,i/10,[num2str(signamp),'mV'],'color',c);
        hold on;
    end;
    grid on;
    plot(exps(exp_num).filter(1).*[1 1],[0 1],'r');
    plot(exps(exp_num).filter(2).*[1 1],[0 1],'r');
    xlabel('Frequency (Hz)');
    ylabel('Attenuation log ratio');
    title(['experiment: ',exps(exp_num).path],'interpreter','none');
    
    %plot phase distortion
    subplot(2,1,2);
    for i=1:size(blocks.onsets,1)
        c=[1-i/4,0,i/3]';
        signamp=mean(blocks.amp(blocks.ampnum==i));
        semilogx(freq_points, max(abs(blocks.phase_error(blocks.ampnum==i,:)')) ,'o-','color',c); % max. distortion across phase
        %plot(log(freq_points), blocks.peak_attenuation(blocks.ampnum==i)./signamp,'o-','color',c);
        text(10,i/10,[num2str(signamp),'mV'],'color',c);
        hold on;
    end;
    xlabel('Frequency (Hz)');
    ylabel('Max abs. harmonic distortion (%)');
    
    grid on;
    
    
    plot2svg(['plots/freq_response',exps(exp_num).path,'.svg']);
    
    saveas(gca,['plots/freq_response',exps(exp_num).path,'.fig'])
    
    
    
    %% plot phase distortion
    
    figure(4+exp_num*10); clf
    
    plotrange=1;
    
    for amp=1:3
        
        
        subplot(3,2,(amp*2)-1);
        hold on;
        d=blocks.phase_error([1:50]+(amp-1)*50,:)';
        for j=1:size(d,2)
            yoffs=j.*.1;
            plot(linspace(-pi,pi,phase_quant_points),d(:,j)+yoffs);
            text(3.3,yoffs,num2str(freq_points(j)))
        end;
        axis tight;
        xlabel('phase (rad)');
        ylabel('phase error (rad)');
        title(['experiment: ',exps(exp_num).path],'interpreter','none');
        
        subplot(3,2,amp*2);
        hold on;
        d=blocks.phase_error([1:50]+(amp-1)*50,:)';
        d(:,end+1)=linspace(-plotrange,plotrange,phase_quant_points);
        
        I=zeros(size(d,1),size(d,2),3);
        
        
        I(:,:,1)=(d>0).*d;
        I(:,:,1)=I(:,:,1)-min(min(I(:,:,1)));
        I(:,:,1)=I(:,:,1)/max(max(I(:,:,1)));
        
        I(:,:,2)=.01;
        I(:,:,2)=(abs(d)<.001).*.1;
        
        I(:,:,3)=(d<0).*-d;
        I(:,:,3)=I(:,:,3)-min(min(I(:,:,3)));
        I(:,:,3)=I(:,:,3)/max(max(I(:,:,3)));
        
        imagesc(sqrt(I));
        
        axis tight;
        colormap gray;
        
        title(['phase error in radians for amp=',num2str( blocks.amp(amp*50)),' scale on right end range=',num2str(plotrange)]);
        
        xlabel('frequency');
        ylabel('phase (rad)');
        
        set(gca,'xTick',[5:10:50])
        set(gca,'xTickLabel', mat2cell(round(freq_points([1:10:50]))) )
        
    end;
    
    pause(1);
    plot2svg(['plots/phase_error',exps(exp_num).path,'.svg']);
    saveas(gca,['plots/phase_error',exps(exp_num).path,'.fig'])
    
end;



