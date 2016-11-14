d:
cd "d:\data\ms\DAISIE"
del DAISIE_1.2_tar.gz
del DAISIE_1.2.zip
R CMD build DAISIE --resave-data
R CMD INSTALL --build DAISIE_1.2.tar.gz
pause
R CMD check --timings --as-cran DAISIE_1.2.tar.gz
pause