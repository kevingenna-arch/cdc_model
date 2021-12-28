
dynare_cleanup ()  {
rtdir=$1
files=($(find . -iname "*.mod"))
for item in ${files[*]}
do
  # treat the string to pull out the file name
  fname=$(basename "$item" .mod)
  ppl="+"
  exts_array=("_results.mat" ".log")
  for ext in ${exts_array[*]}
  do
    ffi="$fname$ext"
    find $1 -type f -iname "$ffi" -exec rm {} +
  done
  # folders
  find $1 -type d -name "$ppl$fname" -exec rm -rf {} + 
  find $1 -type d -name "$fname" -exec rm -rf {} +
done  
}
