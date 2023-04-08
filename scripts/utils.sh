#!/bin/bash

generate_gene_tree()
{
	method=$1 # raxml or iqtree
	trimmed_alignment=$2
	gene_name=$3
	num_core=$4
	gene_tree_f=$5

	echo using \"$method\" to make Max likelihood tree for \"$gene_name\"
	if [[ $method == "raxml" ]]; then
		echo raxml_gene_trees/RAxML_bipartitions.$gene_name.tree > $gene_tree_f
		# took it from baross: 
		# I have to write "-w $(pwd)/gene_trees -n ${gene_name}.tree" instead of "gene_trees/${gene_name}.tree"
		../raxmlHPC -f a -m PROTGAMMAAUTO -p 33 -x 33 -s $trimmed_alignment -w ../raxml_gene_trees -n $gene_name.tree -T $num_core -N 100 # num_bootstrap set to 100
	elif [[ $method == "iqtree" ]]; then
		outfile_prefix=iqtree_gene_trees/$gene_name
		echo ${outfile_prefix}.ufboot > $gene_tree_f
		iqtree -s $trimmed_alignment --nmax 3000 -m MFP -bb 1000 -wbtl -ntmax 10 -nt AUTO -T $num_core --prefix $outfile_prefix -madd C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -mrate E,I,G,I+G,R
	else
		echo "ERROR!! \$method can only be \"iqtree\" or \"raxml\""
		exit 1
	fi
}

run_ecceTERA()
{
	spcies_tree=$1
	gene_tree=$2
	gene_name=$3

	ecceTERA_sym="ecceTERA_output/${gene_name}_symmetric_reconciliation.recPhyloXML"
	if test -f $ecceTERA_sym; then
		echo "$ecceTERA_sym exists, use the old symmetric reconciliation"
	else
		ecceTERA species.file=$spcies_tree gene.file=$gene_tree verbose=true print.reconciliations=1 recPhyloXML.reconciliation=true amalgamate=true
		for type in asymmetric random symmetric; do
			old_f=reconciliationsFile_canonical_${type}.recPhyloXML
			if test -f $old_f; then
				mv $old_f ecceTERA_output/${gene_name}_${type}_reconciliation.recPhyloXML
			else
				# if for some reason ecceTERA fails, stop here
				echo I do not have ---- $old_f ----, I will stop here
				exit 52
			fi
		done
	fi
}


inhouse_scripts_processing()
{
	chronogram=$1
	gene_name=$2

	recPhyloXML_file=ecceTERA_output/${gene_name}_symmetric_reconciliation.recPhyloXML
	chronogram_internal_nodes_f=ecceTERA_analysis/${chronogram%%_sample.chronogram}_${gene_name}_internal_nodes.txt # %% remove suffix
	sym_event_f=ecceTERA_analysis/${gene_name}_symmetric.events
	sym_event_date_f=ecceTERA_analysis/${gene_name}_symmetric.events_event_dates.txt
	birthdate_f=allGeneBirthDate-$(hostname).txt

	python3 scripts/getInternalNodeDates.py $chronogram $recPhyloXML_file
	# produce file: ugam1_ChenParamsEarth_sample.chronogram_internal_nodes.txt
	mv ${chronogram}_internal_nodes.txt $chronogram_internal_nodes_f

	python3 scripts/py3-scripts/recPhyloXMLEventSummary.py -i $recPhyloXML_file -o $sym_event_f --include.transfer.departure
	sed -r -i "s|\s+|\t|g" $sym_event_f

	python3 scripts/recphyloxmlinterpreterspec.py $sym_event_f $chronogram_internal_nodes_f $recPhyloXML_file

	echo $gene_name finish time: $(date +%Y-%m-%d\ %H:%M:%S) >> $birthdate_f
	python3 -u scripts/getGeneBirthDate.py $sym_event_date_f >> $birthdate_f # -u for the "unbuffered" swich for python
	echo ------ >> $birthdate_f
	echo  >> $birthdate_f
}

