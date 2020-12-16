grep $1 $2 | cut -d" " -f2- | tr -d "," > adj_vs_cendiff.dat
python numdiff.py

