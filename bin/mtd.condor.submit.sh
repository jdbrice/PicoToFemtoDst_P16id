######################################################################
#                               requireMtdPid
######################################################################


InitialDir = /home/jdb12/work/PicoToFemtoDst_P16id/bin/
Executable = /home/jdb12/work/PicoToFemtoDst_P16id/bin/femtoMaker.app
Arguments  = /home/jdb12/work/PicoToFemtoDst_P16id/bin/config/requireMtdPid.xml --jobIndex=$(Process)

GetEnv     = True

Queue 120