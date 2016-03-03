
#Make result figures for suc
python analyzeResults.py -i results/varySuc/N10/ -o results/varySuc/N10/builtArrays/ -b 0 -N 30 -w 0 -c 1 -v suc -S 1 -f results/varySuc/N10/ -l 0 -H 0 -p 1 -F results/FinalFigures/

#Make result figures for mdet
python analyzeResults.py -i results/varyDet/N10/suc0.8/ -o results/varyDet/N10/suc0.8/builtArrays/ -b 0 -N 30 -w 0 -c 1 -v mdet -S 1 -f results/varyDet/N10/suc0.8/figures/ -l 1 -H 0 -p 1 -F results/FinalFigures/

#Make result figures for N
python analyzeResults.py -i results/varyN/suc0.8/ -o results/varyN/suc0.8/builtArrays/ -b 0 -N 30 -w 0 -c 1 -v N -S 1 -f results/varyN/suc0.8/figures/ -l 0 -H 0 -p 1 -F results/FinalFigures/
