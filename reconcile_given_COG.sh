#!/bin/bash
COG=$1 # COG=COG2801
num_core=10
gene_tree_method=iqtree
COG_calling_method=diamond

chronogram=ugam1_ChenParamsEarth_sample.chronogram
out_tree_f="tmp/${COG}_${gene_tree_method}_${COG_calling_method}_gene_tree_filename.txt"
only_this_COG_seq=tmp/$COG.faa
gene_alignment=gene_alignments/$COG.afa
trimmed_alignment=gene_alignments/trimmed_$COG.afa

generate_gene_tree()
{
	method=$1 # raxml or iqtree
	trimmed_alignment=$2
	COG=$3
	num_core=$4

	echo using \"$method\" to make Max likelihood tree for \"$COG\"
	if [[ $method == "raxml" ]]; then
		echo raxml_gene_trees/RAxML_bipartitions.$COG.tree > $out_tree_f
		# took it from baross: /usr/local/bin/raxmlHPC-AVX
		# I have to write "-w $(pwd)/gene_trees -n ${COG}.tree" instead of "gene_trees/${COG}.tree"
		./raxmlHPC -f a -m PROTGAMMAAUTO -p 33 -x 33 -s $trimmed_alignment -w $(pwd)/raxml_gene_trees -n $COG.tree -T $num_core -N 100 # num_bootstrap set to 100
	elif [[ $method == "iqtree" ]]; then
		outfile_prefix=$(pwd)/iqtree_gene_trees/$COG
		echo ${outfile_prefix}.ufboot > $out_tree_f
		iqtree -s $trimmed_alignment -m MFP -bb 1000 -wbtl -ntmax 10 -nt AUTO -T $num_core --prefix $outfile_prefix -madd C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4 -mrate E,I,G,I+G,R
	else
		echo "ERROR!! \$method can only be \"iqtree\" or \"raxml\""
	fi
}

inhouse_scripts_processing()
{
	COG=$1
	chronogram=$2
	gene_tree_method=$3
	COG_calling_method=$4
	birthdate_f=${gene_tree_method}-${COG_calling_method}-allGeneBirthDate.txt

	recPhyloXML_file=ecceTERA_output/${COG}_symmetric_reconciliation.recPhyloXML
	chronogram_internal_nodes_f=ecceTERA_analysis/${chronogram%%_sample.chronogram}_${COG}_internal_nodes.txt # %% remove suffix
	sym_event_f=ecceTERA_analysis/${COG}_symmetric.events
	sym_event_date_f=ecceTERA_analysis/${COG}_symmetric.events_event_dates.txt

	python3 getInternalNodeDates.py $chronogram $recPhyloXML_file
	# produce file: ugam1_ChenParamsEarth_sample.chronogram_internal_nodes.txt
	mv ${chronogram}_internal_nodes.txt $chronogram_internal_nodes_f

	python3 py3-scripts/recPhyloXMLEventSummary.py -i $recPhyloXML_file -o $sym_event_f --include.transfer.departure
	sed -r -i "s|\s+|\t|g" $sym_event_f

	python3 recphyloxmlinterpreterspec.py $sym_event_f $chronogram_internal_nodes_f $recPhyloXML_file

	echo $COG >> $birthdate_f
	python3 -u getGeneBirthDate.py $sym_event_date_f >> $birthdate_f # -u for the "unbuffered" swich for python
	echo ------ >> $birthdate_f
}

# sanity check
echo $COG | grep COG
if [[ $? != 0 ]]; then
	echo The COG number $COG must be in the format of "COGXXXX"
else
	echo Starting the program!
fi


if test -f "$only_this_COG_seq"; then
	echo "$only_this_COG_seq exists, use the old file"
else
	python3 grep_seq_given_COG.py $COG $COG_calling_method $only_this_COG_seq
fi


if test -f "$trimmed_alignment"; then
	echo "$trimmed_alignment exists, use the old alignment"
else
	echo aligning the COG category: $COG
	# muscle is localed in this folder, executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
	./muscle5.1 -align $only_this_COG_seq -output $gene_alignment

	# trim out positions with mostly gaps
	~/miniconda3/bin/trimal -in $gene_alignment -out $trimmed_alignment -gt 0.15 # -automated1
fi


out_tree=$(cat $out_tree_f) # using $out_tree_f to store the filename of the gene tree
if test -f "$out_tree"; then
	echo "$out_tree exists, use the old alignment"
else
	generate_gene_tree $gene_tree_method $trimmed_alignment $COG $num_core
fi


ecceTERA_sym="ecceTERA_output/${COG}_symmetric_reconciliation.recPhyloXML"
if test -f $ecceTERA_sym; then
	echo "$out_tree exists, use the old alignment"
else
	ecceTERA species.file=$chronogram gene.file=$out_tree verbose=true print.reconciliations=1 recPhyloXML.reconciliation=true amalgamate=true
	for type in asymmetric random symmetric; do
		mv reconciliationsFile_canonical_${type}.recPhyloXML ecceTERA_output/${COG}_${type}_reconciliation.recPhyloXML
	done
fi


# This has to be executed in Jimmy's conda base environment
# XML has some version changes, only the version in my base environment worked...
inhouse_scripts_processing $COG $chronogram $gene_tree_method $COG_calling_method

# but currently, R is broken in my base environment
source activate anvio-dev
Rscript --vanilla plot-gene-events-histogram.R $COG "iqtree" "diamond"
source deactivate