# ======== include this into your .login (for csh and tcsh) ========

############### setup for VEDA ###############

setenv VEDA '/.../VEDA'
# $VEDA = whole path to VEDA directory

setenv VEDA_BIN 'BIN_gfortran'
# subdirectory within $VEDA that contains binary files

setenv VEDA_COMPILE 'gfortran -Wall -O2'
# compiler and options used to create the binary files

alias veda "python -B $VEDA/mod/Veda.py"

##############################################

