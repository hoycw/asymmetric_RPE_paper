% 4. set parameters to calculate timelocked mutual information
tvec = -1:0.05:2.25; % 50ms steps
twin = 0.2; % window length

Nbins = 2.^[2:1:6];

% assign memory
mutinf = NaN([size(datamtl.trial,1) length(tvec) length(Nbins)]);

for tr = 1:size(datamtl.trial,1)

    for t = 1:length(tvec)

        c1 = squeeze(datamtl.trial(tr,1,nearest(datamtl.time,tvec(t)-twin):nearest(datamtl.time,tvec(t)+twin)));
        c2 = squeeze(datapfc.trial(tr,1,nearest(datamtl.time,tvec(t)-twin):nearest(datamtl.time,tvec(t)+twin)));

        for n = 1:length(Nbins)

            c1bin = rh_discretize(c1, Nbins(n), 'UniCount'); % Uni Count Method
            c2bin = rh_discretize(c2, Nbins(n), 'UniCount');

            mutinf(tr,t,n) = MutualInformation(c1bin,c2bin); clear c1bin c2bin  

        end
        clear c1 c2
    end
end
