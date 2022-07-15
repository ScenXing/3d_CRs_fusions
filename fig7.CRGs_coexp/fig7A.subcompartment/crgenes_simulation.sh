#!/bin/bash
#-Real CR-genes
# intersectBed -a crgenes.bed -b GM12878_subcompartments.bed -wao -f 0.5|awk '$10!="."' >crgenes_in_compartment.bed
intersectBed -a mllgenes.bed -b GM12878_subcompartments.bed -wao -f 0.5|awk '$9!="."' >mllgenes_in_compartment.bed

# Make Control Dataset
#bedtools shuffle -chrom -g hg19.chrom.sizes -i crgenes.bed -seed 11 >crgenes_simu.bed
#intersectBed -a crgenes_simu.bed -b GM12878_subcompartments.bed -wao -f 0.5|awk '$10!="."' >crgenes_simu_in_compartment.bed
#rm 01_scnv5_simu_${i}.bed


