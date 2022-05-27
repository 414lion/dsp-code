# `upsampleSignal.py `

##  Usage

- 显示整体帮助信息

`python upsampleSignal.py -h`

-  显示使用基于滤波器的增采样方法中的FIR滤波器的帮助信息

`python upsampleSignal.py FIR -h`

- 显示使用基于滤波器的增采样方法中的IIR滤波器的帮助信息

`python upsampleSignal.py IIR -h`

- 显示使用基于时域插值的增采样方法的帮助信息

`python upsampleSignal.py Interp -h`

-  使用基于滤波器的增采样方法中的FIR滤波器，不自行设置参数

```help
python upsampleSignal.py FIR INPUT OUTPUT
INPUT      path to the input file(low sampling rate signal).
OUTPUT     path to the output file(after upsampling rate signal).
```

- 使用基于滤波器的增采样方法中的FIR滤波器，自行设置参数

```help
python upsampleSignal.py FIR -c fs up_factor f_p f_st A_s_target window INPUT OUTPUT
fs           Original low sampling frequency(Default=8000Hz).
up_factor    Upsampling factor(Default=6).
f_p          Filter cutoff frequency,default is empty,f_p is in Hz.
f_st         Filter cutoff frequency,default is empty,f_st is in Hz.
A_s_target   Stop band target attenuation,default is 76dB,A_s_target is in db.
window       Filter usage window,Filter usage window,default is 'kaiser'.
```

- 使用基于滤波器的增采样方法中的IIR滤波器，不自行设置参数

```help
python upsampleSignal.py IIR INPUT OUTPUT
INPUT         path to the input file(low sampling rate signal)
OUTPUT        path to the output file(after upsampling rate signal)
```

- 使用基于滤波器的增采样方法中的IIR滤波器，自行设置参数

```help
python upsampleSignal.py IIR -c fs up_factor f_p f_st R_p_target A_s_target ftype INPUT OUTPUT
fs           Original low sampling frequency(Default=8000Hz).
up_factor    Upsampling factor(Default=6)
f_p          Filter cutoff frequency,default is empty,f_p is in Hz.
f_st         Filter cutoff frequency,default is empty,f_st is in Hz.
R_p_target   Passband ripple target,default is 0.1dB,R_p_target is in db.
A_s_target   Stop band target attenuation,default is 76dB,A_s_target is in db.
ftype =      Filter usage dtype,default is 'butter'.
```

- 使用基于时域插值的增采样方法中，不自行设置参数

```help
python upsampleSignal.py Interp INPUT OUTPUT
INPUT         path to the input file(low sampling rate signal)
OUTPUT        path to the output file(after upsampling rate signal) 
```

- 使用基于时域插值的增采样方法中，自行设置参数

```help
python upsampleSignal.py Interp -c fs up_factor kind INPUT OUTPUT
fs           Original low sampling frequency(Default=8000Hz).
up_factor    Upsampling factor(Default=6)
kind         Time domain interpolation method,default is ‘linear’
```



