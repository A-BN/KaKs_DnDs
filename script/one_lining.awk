#!/usr/bin/gawk -f
# Format align on one line per sequence
# usage : gawk -f ./path_to/this_script.awk ./path_to/input_align ./path_to/output_align
BEGIN {}
{
  if ($0 ~ /^>/ && NR == 1) print $0
  else if ($0 ~ /^>/) print "\n"$0
  else printf $0
}
END {}