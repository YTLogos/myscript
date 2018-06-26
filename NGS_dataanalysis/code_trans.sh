#!/bin/bash
sed "{
s/AG/R/g
s/GA/R/g
s/CT/Y/g
s/TC/Y/g
s/GT/K/g
s/TG/K/g
s/AC/M/g
s/CA/M/g
s/CG/S/g
s/GC/S/g
s/AT/W/g
s/TA/W/g
}" GL_snp.txt > GL_snp_clean.txt

