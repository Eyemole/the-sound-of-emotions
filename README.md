# The Sound of Eyemole  

## Notes:  
This code was based off of the code by Chuck Anderson, Dept. of Computer Science, Colorado State University  
(see http://www.cs.colostate.edu/eeg  for more information)   

This code may be copied, distributed, and modified, as long as the author is credited. Also, if you publish this or any derived code on the web, include a link to the above URL.  

## How to Use:  

1. Place any EEG data (2 x n_samples matrix, saved as a .mat file) in the data/eeg_data folder.  
     - You can see the preprocessing for the EEG data in this repository in the emotions_as_audio.ipynb file. Unfortunately, the raw EEG file can't be added here since it's huge, but you can download it free of charge [here](http://headit.ucsd.edu/studies/3316f70e-35ff-11e3-a2a9-0050563f2612/)    
2. Place instrument samples in the data/music_data/{instrument folder name} folder. Instrument folders will be mapped to frequency bands in alphabetical order. Each instrument folder should contain at least one .wav sample - preferably the samples will be ordered as notes on the chromatic scale, but anything goes, as long as they go up in pitch (i.e sample 2 has higher pitch than sample 1)  
3. Open the eeg_to_music_instruments.m script (requires Matlab, sadly), play around with the parameters and run it   
4. Your freshly baked Eyebomination will be delivered to the output folder. Enjoy?   
