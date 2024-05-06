# $1 is the name of your script
# $2 is the posterior distribution of 100 trees
# $3 is the name for the SecSSE output directory, e.g. for size "Size_results"



for i in {1..100}; do

dir=$(echo "Results")
dir+=$(echo "$i")

mkdir Results # For corHMM

# mkdir $dir # For SecSSE

# mkdir $dir/$3 # For SecSSE

Rscript $1 $2 $i 

done