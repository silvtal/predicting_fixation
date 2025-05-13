#!/bin/bash

for f in $(ls 2024_parsed* -d1)
do

cat $f/* > all_temp
sed -i '/"distrib","size","richness","dilfactor","filename","sample","nicheN","nichedist","group_success","success","final_size","evenness","shannon","gini","group_sizes"/d' ./all_temp
sed -i '/"x"/d' ./all_temp
echo '"distrib","size","richness","dilfactor","filename","sample","nicheN","nichedist","group_success","success","final_size","evenness","shannon","gini","group_sizes"' > $f"_all.csv"
cat all_temp >> $f"_all.csv"
rm all_temp

done
