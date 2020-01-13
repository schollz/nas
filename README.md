# Not Another Synthesizer!

This is a notebook tracking my ideas and progress in developing a specific kind of synthesizer.

### Why *another* synthesizer?

I have an OP-1. I love it. It is beautiful in its simplicity and sounds amazing. I want to make something that also sounds amazing and is also extreme in its simplicity (and its okay if its ugly as can be). For me, simplicity also means inexpensive, so I'm going to see if I can design something out of simple components.

I have a list of things I want to exclude which helps guide what I'd like to include.

- [x] No sampling! I want realtime synthesis. I want to generate sounds in as realtime as possible so that I can take advantage of analog inputs that modulate the waveforms.
- [x] Not just a keyboard! I want variable analog inputs. Basically, I want to utilize potentiometers to modulate the waveforms in realtime.
- [x] Not just sin() and cos()! I want extremely customizable synthesizer algorithms. I want to find very simple and beautiful sounding synth algorithms and see if I can tweak them to produce neat variations.
- [x] No ADSR! Just simple continous waves. The volume knob can be used to manually make a ADSR effect (or maybe a theremin style sensor).


### 1. POS-1: grain synthesizer (Arduino based)

My first attempt which yielded something neat is the [grain synth built from an Arduino](https://code.google.com/archive/p/tinkerit/wikis/Auduino.wiki). It sounds pretty good and uses only an Arduino, pots, and a speaker. The sound is good, but sort of low quality.

### 2. POS-2: mozzi synthesizer (Arduino based)

Another Arduino option is to use the [Mozzi synthesizer](https://sensorium.github.io/Mozzi/). I'm not too impressed with the sounds on this and the extensibility. I think I'd like to try with Raspberry Pi audio.

### 3. POS-3: padsynth in C (Linux based)

Moving to a Linux based synth opens a lot of doors in terms of quality.

I've found that the [PADsynth algorithm](https://zynaddsubfx.sourceforge.io/doc/PADsynth/PADsynth.htm) by Nasca Paul is quite beautiful. It can be made polyphonic, and can be made continuous (each is a repeatable unit). Still I'm trying to figure out:

- What kind of sounds can I get?

There are a lot of things that can be done, especially in generating formants. For example, the simple chord uses the [following formant distribution](https://www.wolframalpha.com/input/?i=Plot%5Bexp%28-%28%28x*30.3-600.0%29%2F150.0%29%5E2.0%29%2Bexp%28-%28%28x*30.3-900.0%29%2F250.0%29%5E2.0%29%2Bexp%28-%28%28x*30.3-2200.0%29%2F200.0%29%5E2.0%29%2Bexp%28-%28%28x*30.3-2600.0%29%2F250.0%29%5E2.0%29%2Bexp%28-%28%28x*30.3%29%2F3000.0%29%5E2.0%29*0.1%2Cx%3D%5B0%2C256%5D%5D).

- Can samples be swapped in the middle?
- Should I pre-generate samples?
- How should the sound the outputted in realtime?

For the last question there are a number of possibilities:

- Use JACK realtime server (I tried this and it was painful to switch between pulseaudio and JACK and I also don't find JACK intuitive at all).
- Use ALSA (API seems to be pretty good)
- Use a Raspberry Pi with a MCP4725  ([sampling rate may be too low](https://engineer.john-whittington.co.uk/2015/03/raspberry-pi-dac-mcp4725-with-wiringpi/
))

#### ALSA route

- Try taking the [barebones ALSA example](https://www.alsa-project.org/alsa-doc/alsa-lib/_2test_2pcm__min_8c-example.html) and generate random int and try to play noise after converting [int to unsigned char](https://stackoverflow.com/a/4630167/8133281).
- Try overriding the min buffer for ALSA using [to start ALSA immediately](https://stackoverflow.com/a/25961708/8133281)
- Try [this example for a sine wave](https://albertlockett.wordpress.com/2013/11/06/creating-digital-audio-with-alsa/), and try modulating frequency "on the fly"
- To modulate "on the fly" keep "previous", "current", and "future" in memory
