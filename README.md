# The Sound of Emotions   

## Notes:  
This code was based off of the code by Chuck Anderson, Dept. of Computer Science, Colorado State University  
(see http://www.cs.colostate.edu/eeg  for more information)   

This code may be copied, distributed, and modified, as long as the author is credited. Also, if you publish this or any derived code on the web, include a link to the above URL.  

## How to Use (With Existing Data):  

1. Place any EEG data (2 x n_samples matrix, saved as a .mat file) in the data/eeg_data folder.  
     - You can see the preprocessing for the EEG data in this repository in the emotions_as_audio.ipynb file. Unfortunately, the raw EEG file can't be added here since it's huge, but you can download it free of charge [here](http://headit.ucsd.edu/studies/3316f70e-35ff-11e3-a2a9-0050563f2612/)    
2. Place instrument samples in the data/music_data/{instrument folder name} folder. Instrument folders will be mapped to frequency bands in alphabetical order (i.e the first folder in data/music_data will correspond to the first frequency band you defined and so forth). Each instrument folder should contain at least one .wav sample - anything goes, as long as subsequent samples go up in pitch. The number of samples can be different for each instrument.  
3. Open the eeg_to_music_instruments.m script (requires Matlab, sadly), play around with the parameters and run it   
4. Your freshly baked Eyebomination will be delivered to the output folder. Enjoy?  

## How To Use (Live)   

1. Clone or download this repository.   
2. Install [Muse Direct](http://www.choosemuse.com/developer#direct) (requires Windows 10).   
3. In Muse Direct, click the 'Add' button and choose 'OSC UDP' for the destination. Select the settings exatly as shown in the picture below (except for the nickname, it really doesn't matter):   
![pic](https://user-images.githubusercontent.com/13011161/31425092-9726292e-ae12-11e7-9060-686779c54301.PNG)  
4. In a terminal, change to the the-sound-of-emotions folder and run ```python eeg_to_music_live.py```   
5. :sunglasses:

## Examples:   

You can download samples of what the results sound like from the output folder in this repository. The samples below were generated from EEG data of participants experiencing different emotions:  
- [Disgust](https://github.com/Eyemole/the-sound-of-emotions/blob/master/output/the_sound_of_disgust.wav)   
- [Fear](https://github.com/Eyemole/the-sound-of-emotions/blob/master/output/the_sound_of_fear.wav)    
- [Joy](https://github.com/Eyemole/the-sound-of-emotions/blob/master/output/the_sound_of_joy.wav)    
- [Sadness](https://github.com/Eyemole/the-sound-of-emotions/blob/master/output/the_sound_of_sadness.wav)  
