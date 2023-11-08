expFile="data/exp.txt"
nwkFile="data/Temp_Net.txt"
resD="results_reproduced"
prefix="PropaNet"
TFliFile="data/TF_list.txt"
seedFile="data/bin.txt"
n=$(head -n 1 $seedFile| awk '{print NF}')

mkdir -p ${resD}/intermediate_results
mkdir -p ${resD}/TG

# Create weighted edge template network
python2.7 network_weight.py \
    -nwk ${nwkFile} \
    -exp ${expFile} \
    -o ${resD}/intermediate_results/${prefix}.templateNetwork \
    || exit 1

# Main Propanet algorithm to extract networks and subnetworks for each time step
python2.7 TF_adding_NP_noCtrl.py \
    ${TFliFile} \
    ${resD}/intermediate_results/${prefix}.templateNetwork \
    ${expFile} \
    ${seedFile} \
    -p 10 \
    -cond ${prefix} \
    -outD ${resD}/intermediate_results \
    || exit 1

# Extract target genes of the optimal TFs
for ((i=1;i<=$[$n-1];i++));
do
    python2.7 Target_genes.py \
        ${resD}/intermediate_results/${prefix}.nwk.t$i \
        ${resD}/intermediate_results/${prefix}.DEGli.t$i \
        $TFliFile \
        ${resD}/intermediate_results/${prefix}.TF_rank.t$i.trim \
        ${gSet} \
        ${resD}/TG \
        $i
done

# Final Result : Networks are comprised of resulting TFs/TGs
python2.7 makeTGDesc.py ${resD}