# function to calculate the percentage of common genes between all and paired dataset
percentage <- function(list1,list2){
  common <- intersect(list1$V1,list2$V1)
  tot <- union(list1$V1,list2$V1)
  percent <- (length(common)*100)/length(tot)
  return(percent)
}


### LUAD down -regulated genes

new_down_all_LUAD <- read.table("new_down_all_LUAD.txt") 
new_down_paired_LUAD <- read.table("new_down_paired_LUAD.txt")
per <- percentage(new_down_all_LUAD,new_down_paired_LUAD)
print(per)

### LUAD up-regulated genes
new_up_all_LUAD <- read.table("new_up_all_LUAD.txt")
new_up_paired_LUAD <- read.table("new_up_paired_LUAD.txt")
per <- percentage(new_up_all_LUAD,new_up_paired_LUAD)
print(per)

## down and up together
union_all <- union(new_up_all_LUAD$V1,new_down_all_LUAD$V1)
union_paired <- union(new_up_paired_LUAD$V1,new_down_paired_LUAD$V1)
common <- intersect(union_all,union_paired)
tot <- union(union_paired,union_all)
(length(common)*100)/length(tot)

## LUSC down-regulated genes
new_down_all_LUSC <- read.table("new_down_all_LUSC.txt")
new_down_paired_LUSC <- read.table("new_down_paired_LUSC.txt")
per <- percentage(new_down_all_LUSC,new_down_paired_LUSC)
print(per)

## LUSC up-regulated genes
new_up_all_LUSC <- read.table("new_up_all_LUSC.txt")
new_up_paired_LUSC <- read.table("new_up_paired_LUSC.txt")
per <- percentage(new_up_all_LUSC,new_up_paired_LUSC)
print(per)

## down and up together
union_all <- union(new_up_all_LUSC$V1,new_down_all_LUSC$V1)
union_paired <- union(new_up_paired_LUSC$V1,new_down_paired_LUSC$V1)
common <- intersect(union_all,union_paired)
tot <- union(union_paired,union_all)
(length(common)*100)/length(tot)
