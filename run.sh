#!/bin/bash

data=("email-Enron.txt" "email-Eu-core.txt" "facebook_combined.txt" "u10m_80m.txt" "musae_facebook_edges.csv")

make

for graph in ${data[@]}
do
    echo "Starting ${graph//\.txt/}"
    ./msBFSGraft ../data/${graph} 4
done
