#
# Environment variables for Figaro
#
# Define the root pathnames
setenv fig24home     ~figaro
setenv fig24root     $fig24home/figaro
setenv fig24local    $fig24home/local
setenv fig24srcroot  $fig24home/src
#
setenv FIGARO_DEV    $fig24root/dev/
setenv FIGARO_PROG_L $fig24local/figaro/
setenv FIGARO_PROG_S $fig24root/figaro/
setenv FIGARO_LIBS_S $fig24root/libs
setenv FIGARO_LIBS_L $fig24local/libs
setenv FIGARO_INCLUDES $fig24srcroot/figaro/inc
setenv FIGARO_HELP_DIR $fig24srcroot/local/fighelp
#
# The PGPLOT font file
#
setenv PGPLOT_FONT $fig24local/pgplot/grfont.dat
#
# User variables
# SCRATCH is the directory to place the vars.dst file.
#
if (! $?SCRATCH) then
  if (-e /scr/$user) then
    setenv SCRATCH /scr/${user}/
  else
    setenv SCRATCH $HOME/
  endif
endif
# MTPCKG variable
if (! $?TAPES) then
  setenv TAPES '/dev/rmt20 /dev/rmt16 /dev/rmt12 /dev/rmt8'
endif
# aliases to compile a Figaro program and to release a program
alias fig ${FIGARO_DEV}fig
alias rel ${FIGARO_DEV}rel
alias crepar ${FIGARO_DEV}crepar
alias figaro source ~/bin/figaro.com
