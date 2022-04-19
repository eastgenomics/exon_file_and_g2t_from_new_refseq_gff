#!/bin/bash

# first argument is file from which to remove duplicated transcripts/exon couple
file=$1

# get a default output name if none is provided
if [ ! -z "$2" ]; then
    output_name=$2
else
    file_name_without_extension="${file%.*}"
    output_name="${file_name_without_extension}_no_dups.tsv"
fi

# get transcripts to remove:
# get the tx/exon couples, sort, uniq count them
# get duplicated couples, sort, uniq
transcripts_to_remove=$(cut -f5,6 $file | sort | uniq -c | sort | awk '($1=="2"){print $2}' | sort | uniq)

formatted_transcripts_to_remove=""

# build string for reverse grepping
for tx in $transcripts_to_remove; do
    formatted_transcripts_to_remove+="${tx}|"
done

# remove last character which is a "|" because it breaks the grepping
grep_formatted_transcripts_to_remove=${formatted_transcripts_to_remove::-1}

grep -vE "${grep_formatted_transcripts_to_remove}" $file > $output_name
