#! /bin/csh
#
#
#                       C O N V E X    F I G A R O
#
#   This is the main command file used to start up Figaro.  It
#   modifies your path environment variable to include the Figaro
#   commands, outputs the initial headlines.
#
#   Invocation:
#
#      source ${FIGARO_PROG_L}figaro  Local version
#      source ${FIGARO_PROG_S}figaro  Standard version
#
#    Directory:
#
#      May be invoked from any directory.  The default directory is
#      unchanged (unless by the local or user supplied procedures).
#
#    Environment names required:
#
#      FIGARO_PROG_L   Directory   The local Figaro directory.
#      FIGARO_PROG_S   Directory   The main Figaro directory.
#
# ------------------------------------------------------------------
#
if (! $?FIGARO_PROG_L || ! $?FIGARO_PROG_S) then
   echo "FIGARO_PROG_L or FIGARO_PROG_S environment variables undefined."
   goto QUIT
endif
set FIGARO_PROG_L=$FIGARO_PROG_L
set FIGARO_PROG_S=$FIGARO_PROG_S
#
#echo "  "
#echo \
#"       Unified VMS/Unix F I G A R O  Version 2.4.5 Patch 7 (9th July 1992)"
#echo "  "
#
#   Add Figaro directories to path variable
#
set found=0
foreach directory ($path)
  if ($directory == $FIGARO_PROG_L:h || $directory == $FIGARO_PROG_S:h) then
    set found=1
    break
  endif
end
set found=0
if (! $found) then
   set path=($FIGARO_PROG_L:h $FIGARO_PROG_S:h $path)
endif
#
#   Set up the FIGARO_PROG_U variable if the user has not already
#
if (! $?FIGARO_PROG_U) then
    setenv FIGARO_PROG_U ./
endif
#
#   Add the X binary directory to the path
#
set xwindir=`cat ~figaro/local/xwindir`
set found=0
foreach directory ($path)
  if ($directory == $xwindir) then
    set found=1
    break
  endif
end
if (! $found) then
    set path=($path $xwindir)
endif
#
# Type out latest news
#
#if ( -e ${FIGARO_PROG_L}news.txt ) then
#  cat ${FIGARO_PROG_L}news.txt
#else if ( -e ${FIGARO_PROG_S}news.txt ) then
#   cat ${FIGARO_PROG_S}news.txt
#endif
#
# Check if in batch mode
#
tty -s
if ($status) then
   setenv FIGARO_MODE BATCH
   soft null
else
   setenv FIGARO_MODE INTERACTIVE
endif
#
# Do machine-specific things (this assumes this is only run on Convex & Sun)
# (Convex has no arch command, suns do
#
which arch | grep ^no > /dev/null
if ($status) then
    setenv TV_DEVICE FIGDISP
    setenv FIGDISP 0
    alias figdisp 'figdisp \!* &'
    alias pgdisp 'pgdisp \!* &'
    if (! $?LD_LIBRARY_PATH) then
	setenv LD_LIBRARY_PATH $FIGARO_LIBS_L
    else
	echo $LD_LIBRARY_PATH | grep $FIGARO_LIBS_L > /dev/null
	if ($status) then
	    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${FIGARO_LIBS_L}
	endif
    endif
else
    setenv TV_DEVICE IVAS
    setenv IVAS /dev/ga0
    setenv IMA /dev/im0
endif
#
QUIT:
unset FIGARO_PROG_L
unset FIGARO_PROG_S
unset found
unset xwindir
alias fig ${FIGARO_DEV}fig
alias rel ${FIGARO_DEV}rel
