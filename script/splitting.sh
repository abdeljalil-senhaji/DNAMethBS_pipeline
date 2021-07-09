#!/bin/bash 

# Splitting context:

usage()
{
cat << EOF

usage: $0 <options>
splitting context.
OPTION:
  -h      show this Help message.
  -o      Output.
  -i      Input.
EOF

}

# Get options
while getopts "ho:i:" OPTION
do
  case $OPTION in
    h)  usage; exit 1;;
    o)  output=$OPTARG;;
    i)  input=$OPTARG;;
    ?)  usage; exit;;
  esac
done

# Check that all options were passed

if [[ -z $output ]] || [[ -z $input ]] 
then
  printf "\n=========================\n ERROR: missing options\n=========================\n\n"
  usage
  exit 1
fi


#in_file = snakemake.input["isoforms"]
#out_file = snakemake.output["plot"]

# Fail on the first error
set -e

##########################################################################

file=$(echo $output|rev|cut -d "/" -f 1 |rev)
path=$(echo $output|rev|cut -d "/" -f 2- |rev)

for context in "CHH" "CG" "CHG"; do 
   
   awk "NR<=1 || \$4~/$context/" $input  > $path/$context-$file ;
done


