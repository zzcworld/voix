                               VoiX 1.0.0 Beta 5
                    (http://vocaleliminator.sourceforge.net)
                                  Readme file

                               by Zhang Zhichao
                            (zzcworld@hotmail.com)
================================================================================
* Introduction

  VoiX is an Open Source Vocal Eliminator. There are so many vocal elimination
softwares, but few are open-source. So after trying 2 or 3 this kind of
softwares, I decide to write my own vocal eliminator.

  VoiX is distributed under the GPL license published by the "Free Software
Foundation". So you can modify and/or redistribute it under the terms of GPL.

  VoiX is developed based on VirtualDub's Center Cut filter by Avery Lee, and
written in ANSI C. If you want a Winamp Plugin based on the same algorithm, you
can download "dsp_centercut" at moitah.net.
================================================================================
* Usage

Usage: voix.exe <infile> [outfile] [options]

Options: <lowpass> <highpass>

<lowpass>: Cutoff frequency (in Hz) for low-pass filter
<highpass>: Cutoff frequency (in Hz) for high-pass filter

Example: voix.exe input.wav output.wav 200 8000
================================================================================
* Updates

  1.0.0 Beta 5
    Better realization of "Center Cut" algorithm, Wave file reader improved:
  24 Bit PCM wave and 32 Bit IEEE float wave are supported.

  1.0.0 Beta 4
    Using "Center Cut" algorithm instead of traditional phase substraction;

  1.0.0 Beta 3
    Wave reader improved.

  1.0.0 Beta 4
    Parameter reader supported. A high-pass filter added. Settings improved
    (Pan, Amplitude...).

  1.0.0 Beta 1
    First version: A channel mixer with a low-pass filter.
================================================================================
* TO-DO List

  * When processing some songs, VoiX remains some noise in the song.
    Processing speed needs improving.
================================================================================
* Contact

  If you have any suggestion, please contact me:
  Zhang Zhichao(zzcworld@hotmail.com)
================================================================================
* Special Thanks
  
================================================================================
END OF FILE
