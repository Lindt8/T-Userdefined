# Ray-tracing through materials

This simple user-defined tally is used to identify the path lengths through distinct materials traversed by particles. 
If all source particles are neutrinos (whose reaction physics are disabled by default), this means all of these paths fall on a straight line, or "ray", traced by that source particle.


It is similar to the [T-Cross] tally dump mode in when it writes to its output file. 
This tally returns the distances traveled in each history through each distinct material, as a large list.
Minimally, just the history count, material ID number, and distance traveled from the previous point in centimeters is listed.
Full output with more columns can be printed if the relevant lines in the subroutine are uncommented.

This user-defined tally does not take any input arguments aside from the `file` parameter.
```
[t-userdefined]
file = userdefined.out
```

Example output:
```
#history   matID    dist(cm)
       1    1000  1.0078E+01
       1    5001  8.9842E+00
       1    5004  3.9062E-01
       1    5003  3.9062E-01
       1    5004  3.9062E-01
       1    5005  3.9062E-01
       1    5006  7.8124E-01
       1    5008  3.9062E-01
       1    5006  3.9062E-01
       1    5005  3.9062E-01
       1    5002  7.8124E-01
       1    5008  1.9531E+00
       1    5002  1.4062E+01
       1    5007  3.9062E-01
       1    5008  3.9062E-01
       1    5006  3.9062E-01
       1    5004  3.9062E-01
       1    5003  1.9531E+00
       1    5002  3.9062E-01
       1    5001  6.6405E+00
       2    1000  1.0078E+01
       2    5001  8.9842E+00
       2    5002  3.9062E-01
       2    5003  3.9062E-01
       2    5005  3.9062E-01
       2    5007  3.9062E-01
       2    5009  3.9062E-01
       2    5013  3.9062E-01
       2    5008  3.9062E-01
       2    5002  1.5625E+00
       2    5008  1.9531E+00
       2    5003  3.9062E-01
       2    5002  1.4453E+01
       2    5006  3.9062E-01
       2    5010  3.9062E-01
       2    5008  3.9062E-01
       2    5004  3.9062E-01
       2    5003  1.1719E+00
       2    5004  3.9062E-01
       2    5001  6.6405E+00
       3    1000  1.0078E+01
       3    5001  8.9842E+00
       3    5002  3.9062E-01
       3    5003  3.9062E-01
       3    5004  3.9062E-01
       3    5006  3.9062E-01
       3    5005  3.9062E-01
       3    5008  3.9062E-01
       3    5006  3.9062E-01
       3    5004  3.9062E-01
       3    5002  7.8124E-01
       3    5003  3.9062E-01
       3    5008  1.9531E+00
       3    5006  3.9062E-01
       3    5002  1.4062E+01
       3    5006  3.9062E-01
       3    5008  3.9062E-01
       3    5007  3.9062E-01
       3    5004  3.9062E-01
       3    5003  1.5625E+00
       3    5002  3.9062E-01
       3    5001  6.6405E+00
       4    1000  1.0078E+01
       4    5001  8.9842E+00
       4    5002  3.9062E-01
       4    5003  3.9062E-01
       4    5005  3.9062E-01
       4    5006  3.9062E-01
       4    5009  3.9062E-01
       4    5011  7.8124E-01
       4    5002  1.1719E+00
       4    5008  2.7343E+00
       4    5002  1.4453E+01
       4    5008  7.8124E-01
       4    5006  3.9062E-01
       4    5004  3.9062E-01
       4    5003  1.1719E+00
       4    5004  3.9062E-01
       4    5001  6.6405E+00
       5    1000  1.0078E+01
       5    5001  8.5936E+00
       5    5002  3.9062E-01
       5    5005  3.9062E-01
       5    5003  3.9062E-01
       5    5004  3.9062E-01
       5    5005  3.9062E-01
       5    5008  3.9062E-01
       5    5012  7.8124E-01
       5    5006  3.9062E-01
       5    5007  3.9062E-01
       5    5008  7.8124E-01
       5    5005  3.9062E-01
       5    5006  7.8124E-01
       5    5008  7.8124E-01
       5    5005  3.9062E-01
       5    5002  3.5156E+00
       5    5004  3.9062E-01
       5    5002  1.0547E+01
       5    5006  7.8124E-01
       5    5005  3.9062E-01
       5    5003  3.9062E-01
       5    5004  3.9062E-01
       5    5003  1.1719E+00
       5    5002  3.9062E-01
       5    5001  6.2499E+00
       6    1000  1.0078E+01
       6    5001  8.9842E+00
       6    5003  7.8124E-01
       6    5004  3.9062E-01
       6    5005  7.8124E-01
       6    5007  7.8124E-01
       6    5006  3.9062E-01
       6    5002  1.1719E+00
       6    5008  1.9531E+00
       6    5006  3.9062E-01
       6    5002  1.4062E+01
       6    5007  7.8124E-01
       6    5006  3.9062E-01
       6    5003  1.9531E+00
       6    5002  3.9062E-01
       6    5001  6.6405E+00
```
