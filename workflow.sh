echo "Running TCC";
S=$(python TCC.py $1 $2)
tcc=($(echo $S | tr "\t" "\n"))
if [ "${tcc[0]}" = "True" ];then echo "Temporally consistent"; exit 1; fi
echo "Running PCC";
S=$(python PCC.py $1 $2 -o temp)
S="${S##*$'\n'}"
rm -r temp
pcc=($(echo $S | tr "\t" "\n"))
if [ "${pcc[0]}" = "${tcc[2]}" ];then echo "PCC matches the number of comigrations inferred by MACHINA"; exit 1; fi
echo "Running PCCH";
python PCCH.py $1 $2 -o $3
