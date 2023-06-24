mkdir output
move trim-trial.dat output
move trim-trial.def output

del /Q *.trial*, *.000000, *.swn, *.txt, *.url*, *.log
del /Q BOTNOW, CURNOW, INPUT, PRINT, SWANOUT1
del /Q TMP*, swan*, runid

pause
