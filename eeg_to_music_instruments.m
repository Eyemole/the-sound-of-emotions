
NOTE_PATH = fullfile(pwd, 'data', 'music_data'); % Path to note files 
DATA_PATH = fullfile(pwd, 'data', 'eeg_data'); % Path to EEG data
OUTPUT_PATH = fullfile(pwd, 'output'); % Where to output the file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the filenames below for processing new files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG_FILENAME = 'sad_ica.mat'; % Name of the EEG file
OUTPUT_FILENAME = 'the_sound_of_sadness.wav'; % Name of the output file 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play around with the parameters below to change the sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SCALINGS = [1 6 2 4]; % How long each instrument's sample lasts (relatively)
VOLUME = [1 4 25 100]; % How loud each instrument is, relatively 
                       % note - it makes sense to make instruments that
                       % correspond to higher frequencies louder since the
                       % EEG power spectrum usually takes a 1/f shape (so
                       % higher frequencies have lower amplitude and they
                       % will be drowned out otherwise) 
NUM_TONES = [4 4 1 2]; % The number of different tones played at each time interval for each instrument 
EEG_RANGE = [4 8; 8 12; 12 30; 30 60]; % Frequency ranges (# of frequencies must be <= # of instruments)

PLAYBACK_RATE = 1; % How much to speed up / slow down the file (if < 1, file will be slowed down; if > 1, will be sped up)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
note_tones = get_note_tone_data(NOTE_PATH);
toneData = [];
load(fullfile(DATA_PATH, EEG_FILENAME));

for i = 1:min(256, size(data,1)) % 256 is the channel limit for .wav files in audiowrite 
    [toneData1 Fs] = eegToTones(data(i,:), note_tones(:,:,1:100000), SCALINGS, EEG_RANGE);
    toneData(1:size(toneData1',1),i) = toneData1';
end

savename = fullfile(OUTPUT_PATH, OUTPUT_FILENAME);
audiowrite(savename,toneData,Fs*PLAYBACK_RATE);

function [eegData Fs] = eegToTones(eeg, notes, scalings, volume_factor, num_tones, freq_range, lowf,highf,win,shift,duration)
%  eegToTones(eeg,numTones,lowf,highf,win,shift,duration,filename)
%
%  eeg:         1 x nSamples of EEG (or any) data
%  notes:       nInstruments x nNotes x nSamples files with the vector
%        representations of .wav files 
%  note_tones:  nInstruments x nNotes files with frequencies associated with each note 
%  scalings:    1 x nInstruments vector of the durations for each instrument
%        relative to (duration)
%  freq_ranges: nEEGFrequencies x 2 vector specifying the boundaries of
%        each frequency
%  numTones:    number of maximum amplitude tones to play for each window
%  lowf:        lowest frequency in eeg to map to audible range
%  highf:       highest frequency in eeg to map to audible range
%  win:         window size in eeg to apply FFT to
%  shift:       samples of eeg to shift window, coarseness of spectrogram
%  duration:    seconds (or fraction of) to play tones for each window
%
%  The first two arguments are required.
%
%  This code was based off of the code by Chuck Anderson, Dept. of Computer Science,
%  Colorado State University
%  see http://www.cs.colostate.edu/eeg  for more information
%
%  This code may be copied, distributed, and modified, as long as
%  the author is credited. Also, if you publish this or any derived
%  code on the web, include a link to the above URL.

if nargin < 3, scalings = ones(1, size(notes, 1)); end;
if nargin < 4, volume_factor = logspace(0, 2, size(notes,1)); end;
if nargin < 5, num_tones = ones(1, size(notes,1))*4; end;   
if nargin < 6, freq_range = [ 4 8; 8 12; 12 30; 30 60]; end; % Theta - 4-8 Hz, Alpha - 8-12 Hz, Beta - 12-30 Hz, Gamma - 30-60 Hz
if nargin < 7, lowf = min(min(freq_range)); end
if nargin < 8, highf = max(max(freq_range)); end
if nargin < 9, win = 512; end
if nargin < 10, shift = 64; end
if nargin < 11, duration = shift / 256; end

Fs = 44100; % Most samples are at 44kHz
seconds = duration;

%%%%%%%%%%%%%%%%%
% Rescale all the instrument samples according to scalings
%%%%%%%%%%%%%%%%
for i = 1:size(notes,1)
    scaling = scalings(i);
    for j = 1:size(notes,2)
        rescaledVector = notes(i,j,1:Fs*seconds*scaling);
        notes(i,j,:) = 0;
        notes(i,j,1:size(rescaledVector,3)) = rescaledVector;
    end
end

[signal,f] = myspecgram(eeg,lowf,highf,win,shift);  

MEMORY_LIMIT = 2*1000000;
maxDuration = min(MEMORY_LIMIT/size(notes,1), size(signal,2)*Fs*seconds);
f = f';

frange = f([1 end]);

maxsignal = max(signal(:));
scaled = signal / maxsignal;
tones = [];

for col = 1:size(signal,2)
    
    for eegf = 1:size(freq_range,1)
        
        nwinners = num_tones(eegf);
        toneRows = zeros(1,nwinners);
        [junk,order]=sort(-abs(scaled(:,col)));
        currfreqs = find((f(order) >= freq_range(eegf,1)) & (f(order) <= freq_range(eegf,2)));
        if (length(currfreqs) > 1) 
            order = order(currfreqs);
            winners = order(1:nwinners);

            for i = 1:nwinners
                toneRows(i) = winners(i);
                toneCols = round((col-1)*Fs*seconds+1 : min(maxDuration, (col + scalings(eegf) - 1)*Fs*seconds));
                numNotes = size(find(notes(eegf,:,1)), 2);
                scaledTone = floor((f(toneRows(i)) - freq_range(eegf,1))/(freq_range(eegf,2)-freq_range(eegf,1))*(numNotes-1) + 1); 
                tones(toneRows(i), toneCols) = notes(eegf,scaledTone,1:size(toneCols,2)).*(scaled(toneRows(i), col)*volume_factor(eegf));

            end
        end
        
    end
            
end

melody = sum(tones);
melody = melody / size(tones,1);

eegData = melody;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,f] = myspecgram(eeg,lowf,highf,win,shift)
% [signal,f] = myspecgram(eeg,lowf,highf,win,shift)
%
% [signal,f] = myspecgram(eeg,5,20,512,20);

signal = [];
n = length(eeg);
for w = 1:shift:n-win

  [mag,f] = myfft(eeg(w:w+win));

  signal = [signal mag'];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mag,f] = myfft(data)
Fs = 256;
Fn = Fs/2;

%see http://www.mathworks.com/support/tech-notes/1700/1702.shtml
NFFT=2.^(ceil(log(length(data))/log(2)));
z = fft(data,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
z = z(1:NumUniquePts);
z = abs(z);
z = z*2; %threw out half the points
z(1) = z(1)/2; %but not DC or Nyquist freq
if ~rem(NFFT,2)
  z(length(z)) = z(length(z)) / 2;
end
z = z/length(data);
%now make frequency vector
f=(0:NumUniquePts-1)*2/NFFT;
% Multiply this by the Nyquist frequency 
% (Fn ==1/2 sample freq.)
f=f*Fn;

mag = z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function note_tones = get_note_tone_data(folder_path)
%  get_note_tone_data(folder_path)
%
%   folder_path: the path to the data folder. It is assumed that the data folder
%       contains folders for each instrument. Each instrumed folder is assumed to 
%       contain the .wav files of note samples (or other samples increasing
%       in pitch)
%
%  Returns:  
%    note_tones: a matrix where each column is the vector representation of
%        the signal for each note of each instrument, where note_tones(i, j) 
%        has the signal for note j of instrument i

    instruments = ls(folder_path);
    note_tones = [];
    for i = 3:size(instruments,1) % first 2 outputs of ls will be '.' and '..'
        instrument = strtrim(instruments(i,:));
        notes_files = ls(fullfile(folder_path, instrument, '*.wav'));
        for j = 1:size(notes_files,1)
            notename = strtrim(notes_files(j,:));
            [signal, f] = audioread(fullfile(folder_path, instrument,notename));
            note_tone = signal(:,1);
            note_tones(i-2,j,1:size(note_tone,1)) = note_tone;
        end
    end
end

