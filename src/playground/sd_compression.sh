while getopts 'h:' OPTION ; do
    case "$OPTION" in
	h)
	    printf "Hypergraph:\t\t $OPTARG \n"
	    h=$OPTARG;;
    esac
done

file="$h"
edge_file="$file.neighbor"

#Creating Neighbor-Hypergraph with KaHyPar
/home/theuer/Dokumente/hypergraph/release/src/application/KaHyPar --hgr="$file" --e=0.03 --k=2


g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib /home/theuer/Dokumente/hypergraph/src/playground/sd_compression.cc -o /home/theuer/Dokumente/hypergraph/src/playground/sd_compression -lsdsl -ldivsufsort -ldivsufsort64
/home/theuer/Dokumente/hypergraph/src/playground/sd_compression "$edge_file"


#Get compressing factor of the LZ-Factorization
#edge_file_size=$(stat -c '%s' "$edge_file")
#edge_file_lz_size=$(stat -c '%s' "$edge_file.lz")
#RESULT=$(awk "BEGIN {printf \"%.5f\",${edge_file_size}/${edge_file_lz_size}}")
#echo "Compressed LZ File uses $RESULT times less space than original file"