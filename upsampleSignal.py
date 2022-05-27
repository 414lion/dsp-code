# Standard library
import sys
import argparse
import numpy as np
from scipy import fft, interpolate, signal
import matplotlib.pyplot as plt
import wave
import time
import math
import scipy.io.wavfile as wavfile

# Import everything from `lib/dsp.py`.
from lib.dsp import *

# Un-comment the following line to use interactive matplotlib widget.
# %matplotlib widget
def main():
    args = parse_sys_args()
    args.func(args)

def workflow_FIR(args):
    if args.customization is not None:
        custom = args.customization
        if int(custom[0]) == 0:
            fs, x = read_audio(args.INPUT)
        else:
            fs, x = read_audio(args.INPUT)
            fs = float(custom[0])
        if int(custom[1]) == 0:
            up_factor = 6
        else:
            up_factor = int(custom[1])

        fs_up = fs * up_factor
        x_filter,b,a = upsample_FIR_filter(x, fs, up_factor, float(custom[2]), float(custom[3]), float(custom[4]), custom[5])
        out = remove_delay(x_filter, b, a)
        write_audio(args.OUTPUT, out, fs_up)

    else:
        fs, x = read_audio(args.INPUT)
        up_factor = 6

        fs_up = fs * up_factor
        x_filter,b,a = upsample_FIR_filter(x, fs, up_factor)
        out = remove_delay(x_filter, b, a)
        write_audio(args.OUTPUT, out, fs_up)

def workflow_IIR(args):
    if args.customization is not None:
        custom = args.customization
        if int(custom[0]) == 0:
            fs, x = read_audio(args.INPUT)
        else:
            fs, x = read_audio(args.INPUT)
            fs = float(custom[0])
        if int(custom[1]) == 0:
            up_factor = 6
        else:
            up_factor = int(custom[1])

        fs_up = fs * up_factor
        x_filter,b,a = upsample_IIR_filter(x, fs, up_factor, float(custom[2]), float(custom[3]), float(custom[4]), float(custom[5]), custom[6])
        out = remove_delay(x_filter, b, a)
        write_audio(args.OUTPUT, out, fs_up)

    else:
        fs, x = read_audio(args.INPUT)
        up_factor = 6

        fs_up = fs * up_factor
        x_filter,b,a = upsample_IIR_filter(x, fs, up_factor)
        out = remove_delay(x_filter, b, a)
        write_audio(args.OUTPUT, out, fs_up)

def workflow_interp(args):
    if args.customization is not None:
        custom = args.customization
        if int(custom[0]) == 0:
            fs, x = read_audio(args.INPUT)
        else:
            fs, x = read_audio(args.INPUT)
            fs = float(custom[0])
        if int(custom[1]) == 0:
            up_factor = 6
        else:
            up_factor = int(custom[1])

        fs_up = fs * up_factor
        x_interp = upsample_interp(x, fs, len(x), up_factor, custom[2])
        write_audio(args.OUTPUT, x_interp, fs_up)

    else:
        fs, x = read_audio(args.INPUT)
        up_factor = 6

        fs_up = fs * up_factor
        x_interp = upsample_interp(x, fs, len(x), up_factor)
        write_audio(args.OUTPUT, x_interp, fs_up)

def remove_delay(x, b, a):
    #Function to obtain delay
    w, gd = signal.group_delay((b, a))
    delay = int(gd[0])
    out = x[delay:]
    return out

#Upsampling function based on digital filter - FIR
def upsample_FIR_filter(x, fs, up_factor, f_p = None, f_st = None, A_s_target = None, window = None):
    """
    Realization of signal up-sampling based on digital filter.
    Parameters
    ----------
    x : array
        Low sampling rate signal.
    fs : int,float
         Original low sampling frequency(Default=8000Hz).
    up_factor : int
                Upsampling factor(Default=6).
    f_p : int,float
          Filter cutoff frequency,default is empty,f_p is in Hz.
    f_st : int,float
          Filter cutoff frequency,default is empty,f_st is in Hz.
    A_s_target : int,float
                 Stop band target attenuation,default is 76dB,A_s_target is in db.
    window :
           Filter usage window,default is 'kaiser'.
    Returns
    -------
    x_filter : ndarray
               Upsampled signal realized after filtering.
    b ,a : float
            Filter coefficient.
    Notes
    -----
    Window types:
    - 'boxcar','triang','blackman','hamming','hann'
    -  If it is not for these windows, use signal.kaiserord() to calculate the filter coefficients.
    """

    x_len = len(x)                   #Original signal length
    x_len_up = x_len * up_factor     #Signal length after upsampling
    fs_up = fs * up_factor           #Sampling frequency after upsampling

    # Insert zeros between samples.
    x_insert_0 = np.zeros(x_len_up)
    x_insert_0[::up_factor] = x

    #When fst is not set,set fst.
    if f_st == 0 or f_st is None:
        X = abs(np.fft.fft(x))
        for i in range (int(len(x)/2),int(len(x))):
            if X[i] > 0.001:
                f_st = fs*i/len(x)
                break

    #When fp is not set,set fp.
    if f_p == 0 or f_p is None:
        f_p = fs - f_st                  #The passband cut-off frequency is the sampling frequency minus the stopband cut-off frequency

    #When As is not set,set As.
    if A_s_target == 0 or A_s_target is None:
        A_s_target = round(-20*math.log(0.001/up_factor, 10)) 

    f_tran = f_st - f_p        #set the transition zone

    #When f_p > 3400,set f_p,f_st
    if f_p > 3400:
        f_p = 3400
        f_st = fs - f_p

    f_tran = f_st - f_p

    f_c = (f_p + f_st)/2       #Cut off frequency of filter
    windows = ['boxcar','triang','blackman','hamming','hann'] # windows
    a = 1/up_factor            #Gain I 
    # Compute the order 
    for i in windows:
        # choose the window
        if window == windows:
            if window == 'boxcar':
                window_tran = 1.8*(fs_up/2)  
            if window == 'triang':
                window_tran = 6.1*(fs_up/2)
            if window == 'blackman': 
                window_tran = 11*(fs_up/2)
            if window == 'hamming':
                window_tran = 6.6*(fs_up/2)
            if window == 'hann':
                window_tran = 6.2*(fs_up/2)
            N = int(np.ceil(window_tran/f_tran))
            b = signal.firwin(N, f_c,window = window, fs=fs_up)
        else:
            (N, beta) = signal.kaiserord(A_s_target, f_tran/(0.5*fs_up)) 
            window = ('kaiser', beta)
            b = signal.firwin(N, f_c,window = window, fs=fs_up)
    bands= [
            ('pass', 0, f_p), 
            ('tran', f_p, f_st), 
            ('stop', f_st, fs_up/2)
           ]
    R_p, A_s = analyze_filter(bands, b=b, a=a, show_plot=False, fs=fs_up, tick_format=tick_format_append_hz, amp_in_dB=True)

    while(A_s < A_s_target):
        N = N+1
        b = signal.firwin(N, f_c,window = window, fs=fs_up)
        R_p, A_s = analyze_filter(bands, b=b, a=a, show_plot=False, fs=fs_up, tick_format=tick_format_append_hz, amp_in_dB=True)
    x_filter = signal.lfilter(b=b, a=a, x=x_insert_0)
    return x_filter,b,a

#Upsampling function based on digital filter - IIR
def upsample_IIR_filter(x, fs, up_factor, f_p = None, f_st = None, R_p_target = None, A_s_target = None, ftype = None):
    """
    Realization of signal up-sampling based on digital filter.
    Parameters
    ----------
    x : array
        Low sampling rate signal.
    fs : int,float
         Original low sampling frequency(Default=8000Hz).
    up_factor : int
                Upsampling factor(Default=6).
    f_p : int,float
          Filter cutoff frequency,default is empty,f_p is in Hz.
    f_st : int,float
          Filter cutoff frequency,default is empty,f_st is in Hz.
    R_p_target : int,float
                 Passband ripple target,default is 0.1dB,R_p_target is in db.
    A_s_target : int,float
                 Stop band target attenuation,default is 76dB,A_s_target is in db.
    ftype :
           Filter usage dtype,default is 'butter'.
    Returns
    -------
    x_filter : ndarray
               Upsampled signal realized after filtering.
    b ,a : float
            Filter coefficient.
    Notes
    -----
    ftypes:
    - 'butter','cheby1','cheby2','ellip'
    -  If it is not for these ftypes, use 'cheby2'.
    """

    x_len = len(x)                   #Original signal length
    x_len_up = x_len * up_factor     #Signal length after upsampling
    fs_up = fs * up_factor           #Sampling frequency after upsampling
    
    # Insert zeros between samples.
    x_insert_0 = np.zeros(x_len_up)
    x_insert_0[::up_factor] = x 

    #When fst is not set,set fst.
    if f_st == 0 or f_st is None:
        X = np.fft.fft(x)
        for i in range (int(len(x)/2), int(len(x))):
            if np.abs(X[i]) > 0.001:
                f_st = i/len(x)*fs
                break

    #When fp is not set,set fp.
    if f_p == 0 or f_p is None:
        f_p = fs - f_st                  #The passband cut-off frequency is the sampling frequency minus the stopband cut-off frequency

    #When As is not set,set As.
    if A_s_target == 0 or A_s_target is None:
        A_s_target = round(-20*math.log(0.001/up_factor, 10)) 

    f_tran = f_st - f_p                  #set the transition zone

    #When ftran < 200,set f_p,f_st.
    if f_p > 3400:
        f_p = 3400
        f_st = 4600

    f_tran = f_st - f_p

    #When Rp is not set,set Rp.
    if R_p_target == 0 or R_p_target is None:
        R_p_target = 0.1

    f_c = (f_p + f_st)/2       #Cut off frequency of filter
    ftypes = ['butter','cheby1','cheby2','ellip']
    if ftype not in ftypes:
        ftype = 'cheby2'
    else:
        ftype = ftype
    # Compute the order 
    b, a = signal.iirdesign(f_p, f_st, R_p_target, A_s_target, ftype=ftype, fs=fs_up)
    b = b*up_factor
    bands= [
            ('pass', 0, f_p), 
            ('tran', f_p, f_st), 
            ('stop', f_st, fs_up/2)
            ]
    R_p, A_s = analyze_filter(bands, b=b, a=a, show_plot=False, fs=fs_up, tick_format=tick_format_append_hz, amp_in_dB=True)

    while(A_s < A_s_target):
        A_s_target = A_s_target + 1
        b, a = signal.iirdesign(f_p, f_st,R_p_target, A_s_target,ftype=ftype, fs=fs_up)
        b = b*up_factor
        R_p, A_s = analyze_filter(bands, b=b, a=a, show_plot=False, fs=fs_up, tick_format=tick_format_append_hz, amp_in_dB=True)
    x_filter = signal.lfilter(b=b, a=a, x=x_insert_0)
    return x_filter,b,a

def upsample_interp(x, fs, x_len, up_factor, kind = None):
    kinds = ['nearest','linear']
    if kind not in kinds:
        kind = 'linear'

    fs_up = fs * up_factor
    x_len_up = x_len * up_factor
    x0 = np.linspace(0, x_len, x_len)

    #生成新的xnew数组：插入点数
    xnew = np.linspace(0, x_len, x_len_up)

    #选择时域插值的方法
    f = interpolate.interp1d(x0, x, kind=kind)

    #输入新的xnew，生成时域插值的高采样频率信号
    x_interp = f(xnew)

    return  x_interp

def read_audio(in_file_name):
    return wavfile.read(in_file_name)

def write_audio(out_file_name, data, fs_up):
    out = data.astype(np.int16)
    fs_up = int(fs_up)
    wavfile.write(out_file_name, fs_up, out)

###
# Parse command line arguments.
###
def parse_sys_args():

    # create the top-level parser
    parser = argparse.ArgumentParser(description='Digital audio signal upsampling system.')
    subparsers = parser.add_subparsers(dest='subparser_name',help='sub-command help')

    # create the parser for the "a" command
    parser_a = subparsers.add_parser('FIR', help='a help')
    parser_a.add_argument('INPUT', help='path to the input file(low sampling rate signal)')
    parser_a.add_argument('OUTPUT', help='path to the output file(after upsampling rate signal)')
    help = 'User-defined design metrics.'+'If you don’t want to set the metric just input 0.'+'window choices = {boxcar, triang, blackman, hamming, hann, kaiser}'+'[Default：fs = 8000Hz,up_factor = 6,A_s_target = 76dB,window = kaiser]'
    parser_a.add_argument('-c','--customization', nargs=6, metavar=('fs','up_factor','f_p','f_st','A_s_target','window'), help = help)
    parser_a.set_defaults(func = workflow_FIR)

    # create the parser for the "b" command
    parser_b = subparsers.add_parser('IIR', help='a help')   
    parser_b.add_argument('INPUT', help='path to the input file(low sampling rate signal)')
    parser_b.add_argument('OUTPUT', help='path to the output file(after upsampling rate signal)')
    help = 'User-defined design metrics.'+'If you don’t want to set the metric just input 0.'+'ftype choices = {butter, cheby1, cheby2, ellip}'+'[Default：fs = 8000Hz,up_factor = 6,R_p_target = 0.1dB,A_s_target = 76dB,ftype = ellip]'
    parser_b.add_argument('-c','--customization', nargs=7, metavar=('fs','up_factor','f_p','f_st','R_p_target','A_s_target','ftype'), help = help)
    parser_b.set_defaults(func = workflow_IIR)

    # create the parser for the "c" command
    parser_c = subparsers.add_parser('Interp', help='a help')
    help = 'Specifies the kind of interpolation as a string or as an integer specifying the order of the spline interpolator to use.'
    # parser_c.add_argument('kind', choices=['zero','nearest','slinear','linear','quadratic','cubic'], help = help)
    parser_c.add_argument('INPUT', help='path to the input file(low sampling rate signal)')
    parser_c.add_argument('OUTPUT', help='path to the output file(after upsampling rate signal)')
    help = 'User-defined design metrics.'+'If you don’t want to set the metric just input 0.'+'kind choices = {nearest, linear}'+'[Default：fs = 8000Hz,up_factor = 6,kind = linear]'
    parser_c.add_argument('-c','--customization', nargs=3, metavar=('fs','up_factor','kind'), help = help)
    parser_c.set_defaults(func = workflow_interp)

    if len(sys.argv)==1:
        # No arguments specified.
        parser.print_help()
        parser.exit()
    else:
        args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()