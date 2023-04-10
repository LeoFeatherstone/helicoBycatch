for file in $(ls *.fasta);
do
    iqtree1 -s ${file} -m HKY+G;
done