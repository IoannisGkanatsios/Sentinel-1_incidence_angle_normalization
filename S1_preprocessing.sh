########################################################################################
#############      How to run the script      ##################
#
# First of all, check the location where the gpt is isntalled in the system and change it
# accordingly on the gptPath argument.
#
# bash S1_preprocessing.sh graphs/s1_preprocessing_ocean.xml raw aoi/aoi.txt output/
#
#########################################################################################


gptPath=~/snap/bin/gpt
graph="$1"
source_folder="$2" 
geoRegion="$3" 
target_product="$4"


function get_date_raw () {
    echo $1 | sed -r 's/[^ ]*(_(2.......)T......){2}[^ ]*\.zip/\2/g'
}


function get_source_type () {
	path="$1"
	echo ${path##*/} | grep -o .*_ | cut -d "_" -f1-3
}

printf 'SAR image pre-processing starts......\n'

AOI_coords=("${geoRegion}")
coords=`cat $AOI_coords` 

files=($(find "${source_folder}"/*.zip ))
for ((i=0; i<${#files[*]}; i++))
do
	source_product="${files[$i]}" 
	source_product_date=$(get_date_raw "${source_product}")
	source_product_type=$(get_source_type "${source_product}")
	target_product_name="subset_${source_product_type}_${source_product_date}"
	echo $target_product_name

	${gptPath} "${graph}" \
	-Psource="${source_product}" \
	-PgeoRegion="${coords}" \
	-Ptarget="${target_product}${target_product_name}" 
done

