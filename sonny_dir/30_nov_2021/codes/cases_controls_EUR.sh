#!/bin/bash'

awk '{ if ( $6 == "1" ) print $1,$2 }' HAS_EUR1.fam > controls_EUR.txt

awk '{ if ($6 == "2") print$1,$2 }' HAS_EUR1.fam > cases_EUR.txt

