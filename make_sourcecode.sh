#!/bin/sh
D=sourcecode
mkdir $D
F="makefile values.h export_type.h main.cpp adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd examples"
for X in $F; do
svn export $X $D/$X
done
exit

