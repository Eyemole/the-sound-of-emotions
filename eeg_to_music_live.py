import socket
import os.path
import sys
import time
import os
import argparse
import numpy as np
import pyaudio
import pythonosc_extended.dispatcher as dp
import pythonosc_extended.osc_server as osc_server
from scipy.io import wavfile
from scipy import signal as s

SAMPLE_TIME = 2
SAMPLE_RATE = 256
NUM_CHANNELS = 2
RECORDING_RATE = 44100
pwd = os.getcwd()
# NOTE_PATH = os.path.join(pwd, 'data', 'music_data')

# # How long each instrument's sample lasts (relatively)
# SCALINGS = [1, 6, 2, 4]
# VOLUME = [1, 4, 25, 100]  # How loud each instrument is, relatively
# # The number of different tones played at each time interval for each
# # instrument
# NUM_TONES = [4, 4, 1, 2]
# Frequency ranges (# of frequencies must be <= # of instruments)
# EEG_RANGE = [[4, 8], [8, 12], [12, 30], [30, 59]]
MAX_FREQ = 54
NUM_MAX_FREQS = 4

AUDIBLE_RANGE = [200, 2000]
# How much to speed up / slow down the file (if < 1, file will be slowed
# # down; if > 1, will be sped up)
# PLAYBACK_RATE = 1

# CONVERTED_LENGTH = SAMPLE_TIME * SAMPLE_RATE / RECORDING_RATE


class EEGMusicProperties(object):


	"""Class to hold different EEG Player settings.

	Attributes:
		num_channels (int): The number of channels in the incoming data
		cycle_duration (float): How long (in s) each data collection cycle lasts
		sample_rate (int): The sample frequency of the incoming data
		maximum_frequency (int): Everything above this frequency will be filtered out
		num_max_frequencies (int): The number of highest-power frequencies that will be converted to sound
		converted_range (int[]): The range of frequencies that the converted data can contain
		playback_rate (int): The rate at which converted EEG data is played back
		buffer_size (int): The array size of each
		f (float[]): An array of log-spaced frequencies for the converted data
	"""

	def __init__(self, num_channels=NUM_CHANNELS, cycle_duration=SAMPLE_TIME, sample_rate=SAMPLE_RATE, maximum_frequency=MAX_FREQ, num_max_frequencies=NUM_MAX_FREQS, converted_range=AUDIBLE_RANGE, playback_rate=RECORDING_RATE):

		self.num_channels = num_channels
		self.cycle_duration = cycle_duration
		self.sample_rate = sample_rate
		self.maximum_frequency = maximum_frequency
		self.num_max_frequencies = num_max_frequencies
		self.converted_range = converted_range
		self.playback_rate = playback_rate
		self.buffer_size = sample_rate * cycle_duration
		self.f = self.get_frequency_array()

	def get_frequency_array(self):
		"""Return the log-spaced converted frequency array based on the range of converted frequencies

		Returns:
			float[]: The returned array
		"""
		return np.logspace(np.log10(self.converted_range[0]), np.log10(
			self.converted_range[1]), num=129)[:self.maximum_frequency]


class EEGMusicPlayer(object):


	"""EEG to music converter, containing EEG and converted audio data.

	Attributes:
		data (float64[][]): The raw data from the EEG headset
		output_data (float32[][]): The converted data in float32 format
		stream (PyAudio.Stream()): The actual sound stream
		sample_count (int): The number of samples in the current cycle
		settings (EEGMusicProperties): The parameters for this session
	"""

	def __init__(self, settings):
		self.settings = settings
		self.eegcount = 0
		self.data = np.empty((self.settings.num_channels, self.settings.buffer_size))
		self.output_data = np.empty((self.settings.num_channels, self.settings.playback_rate * self.settings.cycle_duration))
		p = pyaudio.PyAudio()
		self.stream = p.open(format=pyaudio.paFloat32,
							 channels=self.settings.num_channels,
							 rate=self.settings.playback_rate,
							 output=True,
							 stream_callback=lambda in_data, frame_count, time_info, status_flag: (np.array(self.output_data), pyaudio.paContinue))

	def setup_server(self):
		"""Return a UDP OSC server listening to the incoming EEG data.

		By default, the server listens to the local IP on port 5000. Currently, only the raw data from the frontal channels
		(AF7 and AF8) is being processed.

			Returns:
				ThreadingOSCUDPServer: the server listening to EEG data
		"""

		parser = argparse.ArgumentParser()
		ip = socket.gethostbyname(socket.gethostname())
		parser.add_argument("--ip",
							default="127.0.0.1",
							help="The ip to listen on")
		parser.add_argument("--port",
							type=int,
							default=5000,
							help="The port to listen on")
		args = parser.parse_args()

		dispatcher = dp.Dispatcher()
		dispatcher.map("/debug", print)
		dispatcher.map("/muse/eeg", lambda addr, args, ch1, ch2, ch3, ch4, ch5,
					   ch6: self.eeg_handler(addr, args, ch1, ch2, ch3, ch4, ch5, ch6), "EEG")

		server = osc_server.ThreadingOSCUDPServer(
			(args.ip, args.port), dispatcher)
		server.socket.setblocking(0)
		print("Serving on {}".format(server.server_address))
		return server

	def eeg_handler(self, unused_addr, args, ch1, ch2, ch3, ch4, ch5, ch6):
		""" Handle the incoming raw EEG data by adding it to the buffer and converting to audio when the current cycle is complete.

		Only the data from the frontal channels (AF7 and AF8) is currently being used.

		Args:
			unused_addr (str): The address of the incoming OSC message - "/muse/eeg" for this handler. This parameter is not used
			args: Other (unused) args of the incoming OSC message
			ch1 .. ch6 (float64): Raw EEG data from the Muse electrodes. Channels are in following order: TP9, AF7, AF8, TP10, DLR/REF, MUX
		"""

		self.data[0, self.eegcount] = ch2
		self.data[1, self.eegcount] = ch3

		self.eegcount += 1

		if self.eegcount >= self.settings.buffer_size - 1:

			self.eegcount = 0
			old_data = self.data
			for i in range(self.settings.num_channels):
				self.output_data[i] = self.eeg_to_tones(old_data[i])

			self.stream.start_stream()

			while self.stream.is_active():
				time.sleep(self.settings.buffer_size /
						   self.settings.playback_rate)

			self.stream.stop_stream()
			old_data = []

	def eeg_to_tones(self, data):
		""" Returns the converted audio signal by rescaling EEG signal into audible range and restricting the number of frequencies.

		Only the top {settings.num_max_frequencies} frequencies are kept in the audio signal (according to their spectral power over this recording cycle). The volume of the
		resulting signal is modulated by spectral power.

		Note: The conversion to float32 is due to PyAudio not supporting 64-bit float, which is the format that Muse Direct sends data in over OSC/

		Args:
			data(float64[]): The raw EEG data that will be converted into audio

		Returns:
			float32[]: The transformed EEG data
		"""

		_, _, signal = s.spectrogram(data, self.settings.sample_rate)

		maxsignal = np.max(signal)
		signal = np.mean(signal / maxsignal,
						 axis=1)[:self.settings.maximum_frequency]

		endindex = self.settings.maximum_frequency - 1
		startindex = endindex - self.settings.num_max_frequencies

		highest_freqs = self.settings.f[np.argpartition(
			signal, self.settings.num_max_frequencies)][startindex:endindex]
		signal = signal[np.argpartition(signal, self.settings.num_max_frequencies)][
										startindex:endindex]

		print(highest_freqs)
		values = np.sum([np.sin(2 * np.pi * np.arange(self.settings.playback_rate * self.settings.cycle_duration)
								* highest_freqs[i] / self.settings.playback_rate) * signal[i] for i in range(self.settings.num_max_frequencies)], axis=0).astype(np.float32)
		return values


if __name__ == "__main__":

	settings = EEGMusicProperties()
	eeg = EEGMusicPlayer(settings)
	server = eeg.setup_server()
	while 1:
		server.handle_request()
