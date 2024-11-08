#### Biolog data analysis for permafrost mixing study ####
#### Stacey Doherty ####
#### November 7, 2024 ####



library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(splitstackshape)
library(ggplot2)
library(devtools)
library(data.table)
library(vegan)


##### import and process tf data #####
setwd("/path_here") 
#change to directory containing the .exl data files
F_path <- "/path_to_exl_files" 


#get the names of the files
files <- list.files(F_path, pattern = ".EXL")

#make dataframe that will house all values
df <- data.frame(matrix(ncol = 4, nrow = 1))
colnames(df) <- c('plate', 'read_time', 'plate_pos', 'OD')


##This for loop will look at each file individually, pull out the 
for(file1 in files){
  txt <- readLines(file1)
# create an index for the lines that are needed
t1 <- which(grepl('user_field.1=', txt))
t2 <- which(grepl('time_start=', txt))
t3 <- which(grepl('plate_position=', txt))
#select 2nd to last line that contains the OD measurements
OD <- txt[70]
#remove the delimiters (besides the commas)
OD1 <- gsub('[\"]', '', OD)
#remove the interval number from the beginning of the string (its the same in all files)
OD1 <- sub('..','', OD1)


# filter the text with the index 'ti' 
# and remove label 
lst1 <-txt[t1]
lst1 <- str_remove(lst1, "user_field.1=")

lst2 <-txt[t2]
lst2 <- str_remove(lst2, "time_start=")

lst3 <-txt[t3]
lst3 <- str_remove(lst3, "plate_position=")


#append to master dataframe
df <- rbind(df, list(lst1, lst2, lst3, OD1), stringsAsFactors = TRUE)
#provide column names

}

#end loop


#split last column into columns via comma separated 
dfsplit <- cSplit(df, 'OD', sep = ",", type.convert = TRUE)
#remove the first empty Row
dfsplit <- dfsplit[-1,]
#convert date time to correct format
dfsplit[["read_time"]] <- as.POSIXct(dfsplit[["read_time"]],format="%m/%d/%Y %I:%M:%S %p")




#create column with incubation time elapsed
#create column with start time
dfsplit$start_time = min(dfsplit$read_time)
dfsplit <- relocate(dfsplit, start_time, .after = read_time)

#calculate incubation length of each reading and make new column of the data
dfsplit$inc_time_s = dfsplit$read_time - dfsplit$start_time
dfsplit <- relocate(dfsplit, inc_time_s, .after = start_time)
dfsplit$inc_time_s <- as.numeric(dfsplit$inc_time_s) #this is time in seconds
#convert seconds to hours
dfsplit$inc_time_h = dfsplit$inc_time_s / 3600
dfsplit <- relocate(dfsplit, inc_time_h, .after = inc_time_s)
#convert to days
dfsplit$inc_time_d = dfsplit$inc_time_s / 86400
dfsplit <- relocate(dfsplit, inc_time_d, .after = inc_time_h)
#make column called "time" 
dfsplit$time = dfsplit$inc_time_d
dfsplit <- relocate(dfsplit, time, .after = inc_time_d)

#make "blank" column with the average OD in the water wells
dfsplit$blank <- rowMeans(subset(dfsplit, select=c(OD_01, OD_05, OD_09)))
dfsplit <- relocate(dfsplit, blank, .after = plate_pos)

#rename column headers to work with growthcuver
#you will need a csv file with the old labels (e.g., "OD_01") in column 1 and new labels (e.g., "A1") in column 2
label_chg <- read.csv("/path.csv")
#rename based on column ID match
setnames(dfsplit, label_chg$old, label_chg$new)
#rearrange OD value columns so columns are ordered correctly for growthcurver
dfsplit <- dfsplit[,c("plate", "read_time", "start_time","inc_time_s","inc_time_h","inc_time_d","time","plate_pos", "blank","A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2","A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3","A4", "B4", "C4", "D4", "E4", "F4", "G4", "H4","A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5","A6", "B6", "C6", "D6", "E6", "F6", "G6", "H6","A7", "B7", "C7", "D7", "E7", "F7", "G7", "H7","A8", "B8", "C8", "D8", "E8", "F8", "G8", "H8","A9", "B9", "C9", "D9", "E9", "F9", "G9", "H9","A10", "B10", "C10", "D10", "E10", "F10", "G10", "H10","A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11","A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12")]


#split the dataframe by unique plate IDs, this returns a list of dataframes, one for each plate
plates <- split(dfsplit, dfsplit$plate)
#str(plates)
#view the first dataframe in the list 
view(plates[[1]])


##### ANALYZE DATA - calculate Average well color development (AWCD) and substrate specific (SAWCD) of each plate  ####

#make dataframe that will house AWCD of all plates
df_AWCD <- data.frame(matrix(ncol = 9, nrow = 1))
colnames(df_AWCD) <- c('plate', 'inc_time_d', 'AWCD','SAWCD_amino', 'SAWCD_amine', 'SAWCD_carb', 'SAWCD_carbacid', 'SAWCD_pheno', 'SAWCD_poly')

#loop through the plates list to generate a new df with AWCD for each plate.
for(plate1 in plates){
  d <- plate1
  #correct for blanks by subtracting the averaged blank value from all wells 
  d[,10:ncol(d)] <- d[,10:ncol(d)] - d$blank
 
  # calculate AWCD for each timepoint (row) 
  # Columns to sum
  columns_to_sum <- d[,10:ncol(d)]
  # Compute the sum of specific columns
  sum_values <- rowSums(columns_to_sum)
  # Append sum as a new column to the dataframe
  d$AWCD <- sum_values/31

  
  #calculate SAWCD for each substrate class - note: this is for Ecoplates
  #identify columns by substrate class 
  c_amino_sum <- subset(d, select=c("A4","B4","C4","D4", "E4","F4","A8","B8","C8","D8", "E8","F8","A12","B12","C12","D12","E12","F12"))
  c_amine_sum <- subset(d, select=c("G4","H4","G8","H8","G12","H12"))
  c_carb_sum <- subset(d, select=c("G1","H1","A2","B2","C2","D2","E2","G2","H2","A3","G5","H5","A6","B6","C6","D6","E6","G6","H6","A7","G9","H9","A10","B10","C10","D10","E10","G10","H10","A11"))
  c_carbacid_sum <- subset(d, select=c("B1","F2","B3","E3","F3","G3","H3","B5","F6","B7","E7","F7","G7","H7","B9","F10","B11","E11","F11","G11","H11"))
  c_pheno_sum <- subset(d, select=c("C3","D3","C7","D7","C11","D11"))
  c_poly_sum <- subset(d, select=c("C1","D1","E1","F1","C5","D5","E5","F5","C9","D9","E9","F9"))
  #compute sum of specific substrate classes
  amino_sum <- rowSums(c_amino_sum)
  amine_sum <- rowSums(c_amine_sum)
  carb_sum <- rowSums(c_carb_sum)
  carbacid_sum <- rowSums(c_carbacid_sum)
  pheno_sum <- rowSums(c_pheno_sum) 
  poly_sum <- rowSums(c_poly_sum)
  # Append sum as a new column to the dataframe 
  d$SAWCD_amino <- amino_sum/18
  d$SAWCD_amine <- amine_sum/6 
  d$SAWCD_carb <- carb_sum/30 
  d$SAWCD_carbacid <- carbacid_sum/21
  d$SAWCD_pheno <- pheno_sum/6  
  d$SAWCD_poly <- poly_sum/12
 
  #append plate metadata, timepoint, AWCD, SAWCD to new dataframe containing all plates 
  df_AWCD <- rbind(df_AWCD, list(d$plate, d$inc_time_d, d$AWCD, d$SAWCD_amino, d$SAWCD_amine, d$SAWCD_carb, d$SAWCD_carbacid, d$SAWCD_pheno, d$SAWCD_poly))
  
}

#remove the first empty Row
df_AWCD <- df_AWCD[-1,]
#shorten time to whole number
df_AWCD$inc_time_d <- round(df_AWCD$inc_time_d, digits = 0)
#make metadata columns of AWCD table
df_AWCD <- separate(df_AWCD, plate, into = c('site','mix_ratio','rep','timepoint'), sep='-', remove = FALSE, fill = "right", extra = "merge")



#create necessary columns 
df_AWCD$site_mix <- paste(df_AWCD$site, "-", df_AWCD$mix_ratio)
df_AWCD$time <- paste(df_AWCD$timepoint, "-", df_AWCD$inc_time_d)
df_AWCD$site_mix_time <- paste(df_AWCD$site, "-", df_AWCD$mix_ratio, "-", df_AWCD$timepoint, "-", df_AWCD$inc_time_d)

 
#calculate average and SE of AWCD in summary table
AWCD_summary <- df_AWCD %>%
  group_by(site_mix_time) %>%
  summarise(
    sd = sd(AWCD, na.rm = TRUE),
    AWCD_avg = mean(AWCD),
    n=n(),
    se= sd/sqrt(n))
AWCD_summary <- separate(AWCD_summary, site_mix_time, into = c('site','mix_ratio','timepoint', 'inc_time'), sep='-', remove = FALSE)
#clean up the columns
AWCD_summary$inc_time <- as.numeric(AWCD_summary$inc_time)
AWCD_summary$mix_ratio <- str_trim(AWCD_summary$mix_ratio)
AWCD_summary$mix_ratio <- factor(AWCD_summary$mix_ratio, levels = c("100AL", "10PF_90AL" ,"50PF_50AL", "90PF_10AL", "100PF", "100TZ"))
AWCD_summary$mix_ratio_time <- paste(AWCD_summary$mix_ratio, "-", AWCD_summary$timepoint)


#plot AWCD over time, averaged across replicate plates
hex_colors <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4", "#91cf60") 


p_AWCD_tf <- ggplot(AWCD_summary, aes(x=inc_time, y=AWCD_avg, color=mix_ratio, shape = timepoint)) +
  geom_line(aes(group = mix_ratio_time, linetype = timepoint)) + geom_point(size = 2) + geom_errorbar(aes(ymin=AWCD_avg-se, ymax=AWCD_avg+se), width=.2, position=position_dodge(0.05)) + facet_grid(.~site) + theme_bw() + theme(
    text = element_text(size = 26),      
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14))+ scale_color_manual(values = hex_colors) + xlab("Incubation Time (d)")+ylab("AWCD")

#change x axis to 40 days
p_AWCD_tf + scale_x_continuous(limits = c(0, 40))  # Adjust 40 to your desired maximum value  



###### ANALYZE DATA - Bar charts of SAWCD ##### 

#reformat dataframe

df_SAWCD <- filter(df_AWCD,timepoint == "tf" & inc_time_d == 39)
df_SAWCD <- pivot_longer(df_SAWCD, cols = c(8:13), names_to = "sub_class", values_to = "SAWCD")
df_SAWCD$sub_class <- as.factor(df_SAWCD$sub_class)


#function to get average and standard deviation
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


df_SAWCD_1_stats <- data_summary(df_SAWCD, varname="SAWCD", 
                                 groupnames=c("mix_ratio", "sub_class", "site"))

df_SAWCD_1_stats$sub_class <- reorder(df_SAWCD_1_stats$sub_class, df_SAWCD_1_stats$SAWCD)
df_SAWCD_1_stats$sub_class <- factor(df_SAWCD_1_stats$sub_class, levels=rev(levels(df_SAWCD_1_stats$sub_class)))

#### SAWCD barplots ####

#barplot with errorbars - substrate class x-axis

level_order <- c("100AL", "10PF_90AL", "50PF_50AL",  "90PF_10AL", "100PF", "100TZ")
df_SAWCD_1_stats$mix_ratio <- factor(df_SAWCD_1_stats$mix_ratio, levels=level_order)
df_SAWCD_1_stats$sub_class <- factor(df_SAWCD_1_stats$sub_class, levels=c("SAWCD_amine", "SAWCD_amino", "SAWCD_carb",  "SAWCD_carbacid", "SAWCD_pheno", "SAWCD_poly"))


errorbar_sawcd <- ggplot(df_SAWCD_1_stats, aes(fill=mix_ratio, y=SAWCD, x=sub_class)) + 
  geom_bar(stat = "summary", position=position_dodge())+
  geom_errorbar(aes(ymin=SAWCD-sd, ymax=SAWCD+sd), width=.2,
                position=position_dodge(.9))+
  scale_fill_manual(breaks = c("100AL", "10PF_90AL", "50PF_50AL",  "90PF_10AL", "100PF", "100TZ"),
                    limits = c("100AL", "10PF_90AL", "50PF_50AL",  "90PF_10AL", "100PF", "100TZ"), 
                    values = hex_colors,
  )+
  facet_wrap(~site)+
  
  scale_x_discrete(labels=c("SAWCD_amino" = "amino acid", "SAWCD_amine" = "amine", 
                            "SAWCD_carb" = "carbohydrate", "SAWCD_carbacid" = "carboxylic acid",
                            "SAWCD_pheno"= "phenolic compound", "SAWCD_poly" = "polymer"))+
  theme_bw()+
  xlab("Substrate Class")+
  labs(fill="Mix Ratio")

errorbar_sawcd + theme(axis.text.x = element_text(angle = 20, hjust = 1 ), text = element_text(size = 26),      
                                                  axis.title = element_text(size = 16), 
                                                  axis.text = element_text(size = 14),  
                                                  legend.title = element_text(size = 16), 
                                                  legend.text = element_text(size = 14))



##### ANALYZE DATA - calculate % richness and shannon diversity of utilized substrates of each plate at endpoint reading #### 

#average technical replicates for each plate by substrate replicate wells

#make dataframe that will house averaged well values for each substrate of all plates
df_avg_subs <- data.frame(matrix(ncol = 34, nrow = 1))
colnames(df_avg_subs) <- c('plate', 'inc_time_d', "c_water", "c_phenylethylamine", "c_putrescine", "c_glycyl_l_glutamic_acid", 
                                                   "c_l_arginine", "c_l_asparagine", "c_l_phenylalanine", "c_l_serine", 
                                                   "c_l_threonine", "c_d_cellobiose", "c_d_mannitol", "c_d_xylose", 
                                                   "c_i_erythritol", "c_n_acetyl_d_glucosamine", "c_a_d_lactose", 
                                                   "c_b_methyl_d_glucoside", "c_2_hydroxy_benzoic_acid", 
                                                   "c_4_hydroxy_benzoic_acid", "c_d_galactonic_acid_lactone", 
                                                   "c_d_galacturonic_acid", "c_d_glucosaminic_acid", "c_d_malic_acid", 
                                                   "c_d_itaconic_acid", "c_d_pyruvic_acid", "c_a_ketobutyric_acid", 
                                                   "c_hydroxybutyric_acid", "c_glycerol_phosphate", "c_glucose_l_phosphate", 
                                                   "c_glycogen", "c_tween_40", "c_tween_80", "c_cyclodextrin")


#loop through the plates list to generate a new df with averaged substrates for each plate.
for(plate1 in plates){
  d <- plate1
  #correct for blanks by subtracting the averaged blank value from all wells 
  d[,10:ncol(d)] <- d[,10:ncol(d)] - d$blank
  
  # calculate average substrates for each timepoint (row) 
  #identify columns by substrate
  c_water <- subset(d, select=c("A1","A5","A9"))
  c_phenylethylamine <- subset(d, select=c("G4","G8","G12"))
  c_putrescine <- subset(d, select=c("H4","H8","H12"))
  c_glycyl_l_glutamic_acid <- subset(d, select=c("F4","F8","F12"))
  c_l_arginine<- subset(d, select=c("A4","A8","A12"))
  c_l_asparagine <- subset(d, select=c("B4","B8","B12"))
  c_l_phenylalanine <- subset(d, select=c("C4","C8","C12"))
  c_l_serine <- subset(d, select=c("D4","D8","D12"))
  c_l_threonine <- subset(d, select=c("E4","E8","E12"))
  c_d_cellobiose <- subset(d, select=c("G1","G5","G9"))
  c_d_mannitol <- subset(d, select=c("D2","D6","D10"))
  c_d_xylose <- subset(d, select=c("B2","B6","B10"))
  c_i_erythritol <- subset(d, select=c("C2","C6","C10"))
  c_n_acetyl_d_glucosamine <- subset(d, select=c("E2","E6","E10"))
  c_a_d_lactose <- subset(d, select=c("H1","H5","H9"))
  c_b_methyl_d_glucoside <- subset(d, select=c("A2","A6","A10"))
  c_2_hydroxy_benzoic_acid <- subset(d, select=c("C3","C7","C11"))
  c_4_hydroxy_benzoic_acid <- subset(d, select=c("D3","D7","D11"))
  c_d_galactonic_acid_lactone  <- subset(d, select=c("A3","A7","A11"))
  c_d_galacturonic_acid <- subset(d, select=c("B3","B7","B11"))
  c_d_glucosaminic_acid <- subset(d, select=c("F2","F6","F10"))
  c_d_malic_acid <- subset(d, select=c("H3","H7","H11"))
  c_d_itaconic_acid <- subset(d, select=c("F3","F7","F11"))
  c_d_pyruvic_acid <- subset(d, select=c("B1","B5","B9"))
  c_a_ketobutyric_acid <- subset(d, select=c("G3","G7","G11"))
  c_hydroxybutyric_acid <- subset(d, select=c("E3","E7","E11"))
  c_glycerol_phosphate <- subset(d, select=c("H2","H6","H10"))
  c_glucose_l_phosphate<- subset(d, select=c("G2","G6","G10"))
  c_glycogen <- subset(d, select=c("F1","F5","F9"))
  c_tween_40<- subset(d, select=c("C1","C5","C9"))
  c_tween_80<- subset(d, select=c("D1","D5","D9"))
  c_cyclodextrin <- subset(d, select=c("E1","E5","E9"))
  
  #store all dataframes in a list
  subs_df_list <- list(c_water,c_phenylethylamine,c_putrescine,c_glycyl_l_glutamic_acid,c_l_arginine, c_l_asparagine,c_l_phenylalanine, c_l_serine, c_l_threonine, c_d_cellobiose,
                       c_d_mannitol, c_d_xylose, c_i_erythritol, c_n_acetyl_d_glucosamine, c_a_d_lactose, c_b_methyl_d_glucoside, c_2_hydroxy_benzoic_acid, c_4_hydroxy_benzoic_acid,
                       c_d_galactonic_acid_lactone, c_d_galacturonic_acid, c_d_glucosaminic_acid, c_d_malic_acid, c_d_itaconic_acid, c_d_pyruvic_acid, c_a_ketobutyric_acid,
                       c_hydroxybutyric_acid, c_glycerol_phosphate, c_glucose_l_phosphate, c_glycogen, c_tween_40, c_tween_80, c_cyclodextrin) 
  
  # Assign names to the list elements
  names(subs_df_list) <- c("c_water", "c_phenylethylamine", "c_putrescine", "c_glycyl_l_glutamic_acid", 
                           "c_l_arginine", "c_l_asparagine", "c_l_phenylalanine", "c_l_serine", 
                           "c_l_threonine", "c_d_cellobiose", "c_d_mannitol", "c_d_xylose", 
                           "c_i_erythritol", "c_n_acetyl_d_glucosamine", "c_a_d_lactose", 
                           "c_b_methyl_d_glucoside", "c_2_hydroxy_benzoic_acid", 
                           "c_4_hydroxy_benzoic_acid", "c_d_galactonic_acid_lactone", 
                           "c_d_galacturonic_acid", "c_d_glucosaminic_acid", "c_d_malic_acid", 
                           "c_d_itaconic_acid", "c_d_pyruvic_acid", "c_a_ketobutyric_acid", 
                           "c_hydroxybutyric_acid", "c_glycerol_phosphate", "c_glucose_l_phosphate", 
                           "c_glycogen", "c_tween_40", "c_tween_80", "c_cyclodextrin")
  
  
  # Iterate and add new columns based on dataframe names
  for (name in names(subs_df_list)) {
    # Access the dataframe by name
    sub <- subs_df_list[[name]]
    # Compute the sum of specific substrate classes
    sub_sum <- rowSums(sub) / 3
    
    # Create a unique column name for the sum based on the dataframe name
    column_name <- paste0(name)
    
    # Append sum as a new column to the dataframe d
    d[[column_name]] <- sub_sum
  }
  
  #append plate metadata, timepoint, and averaged substrate values new dataframe containing all plates 
  df_avg_subs <- rbind(df_avg_subs, list(d$plate, d$inc_time_d, d$c_water, d$c_phenylethylamine, d$c_putrescine, d$c_glycyl_l_glutamic_acid, d$c_l_arginine, d$c_l_asparagine, d$c_l_phenylalanine, d$c_l_serine, d$c_l_threonine, d$c_d_cellobiose,
                                         d$c_d_mannitol, d$c_d_xylose, d$c_i_erythritol, d$c_n_acetyl_d_glucosamine, d$c_a_d_lactose, d$c_b_methyl_d_glucoside, d$c_2_hydroxy_benzoic_acid, d$c_4_hydroxy_benzoic_acid,
                                         d$c_d_galactonic_acid_lactone, d$c_d_galacturonic_acid, d$c_d_glucosaminic_acid, d$c_d_malic_acid, d$c_d_itaconic_acid, d$c_d_pyruvic_acid, d$c_a_ketobutyric_acid,
                                         d$c_hydroxybutyric_acid, d$c_glycerol_phosphate, d$c_glucose_l_phosphate, d$c_glycogen, d$c_tween_40, d$c_tween_80, d$c_cyclodextrin))
  
}


#remove the first empty Row
df_avg_subs <- df_avg_subs[-1,]
#shorten time to whole number

df_avg_subs$inc_time_d <- round(df_avg_subs$inc_time_d, digits = 0)
#make metadata columns of AWCD table
df_avg_subs <- separate(df_avg_subs, plate, into = c('site','mix_ratio','rep','timepoint'), sep='-', remove = FALSE, fill = "right", extra = "merge")



#create necessary columns 
df_avg_subs$site_mix <- paste(df_avg_subs$site, "-", df_avg_subs$mix_ratio)
df_avg_subs$time <- paste(df_avg_subs$timepoint, "-", df_avg_subs$inc_time_d)
df_avg_subs$site_mix_time <- paste(df_avg_subs$site, "-", df_avg_subs$mix_ratio, "-", df_avg_subs$timepoint, "-", df_avg_subs$inc_time_d)
#remove water data if not interested
df_avg_subs <- df_avg_subs[, !(names(df_avg_subs) %in% "c_water")]



#make dataframe that will house alpha diversity measures of all plates
df_diversity <- data.frame(matrix(ncol = 4, nrow = 1))
colnames(df_diversity) <- c('plate', 'inc_time_d', 'richness','shannon')

# convert negative values to zero (substrate not utilized)
df_avg_subs[,7:ncol(df_avg_subs)] <- lapply(df_avg_subs[,7:ncol(df_avg_subs)], function(x) ifelse(x < 0, 0, x))
  
# calculate % richness for at each timepoint
df_avg_subs$richness <- rowSums(df_avg_subs[,7:37] > 0) 

    
# calculate shannon diversity at each time
  
df_avg_subs$shannon <- diversity(df_avg_subs[,7:37], index = "shannon")

#append plate metadata, timepoint, diversity metrics to new dataframe containing alpha diversity metrics of all plates 
df_diversity <- rbind(df_diversity, list(df_avg_subs$plate, df_avg_subs$inc_time_d, df_avg_subs$richness, df_avg_subs$shannon))
  


#remove the first empty Row
df_diversity <- df_diversity[-1,]

#shorten time to just one decimal
df_diversity$inc_time_d <- round(df_diversity$inc_time_d, digits = 0)
#make metadata columns of AWCD table
df_diversity <- separate(df_diversity, plate, into = c('site','mix_ratio','rep','timepoint'), sep='-', remove = FALSE, fill = "right", extra = "merge")



df_diversity$site_mix <- paste(df_diversity$site, "-", df_diversity$mix_ratio)
df_diversity$mix_time <- paste(df_diversity$mix_ratio, "-", df_diversity$timepoint)


#clean up data
df_diversity$site <- factor(df_diversity$site, levels = c("APT", "TK1"))
df_diversity$mix_ratio <- factor(df_diversity$mix_ratio, levels = c("100AL", "10PF_90AL" ,"50PF_50AL", "90PF_10AL", "100PF", "100TZ"))
df_diversity$mix_time <- factor(df_diversity$mix_time, levels = c("100AL - t0", "100AL - tf", "10PF_90AL - tf" ,"50PF_50AL - tf", "90PF_10AL - tf", "100PF - t0", "100PF - tf", "100TZ - t0", "100TZ - tf"))

#grab the final endpoint data  for plotting 
#max of timepoint column

df_diversity_final_tp <- subset(df_diversity, inc_time_d == 39)



#plot shannon boxplot
p_shan <- ggplot(df_diversity_final_tp, aes(x=mix_time, y=shannon, fill=mix_ratio)) + geom_boxplot(na.rm = TRUE) + facet_grid(. ~ site , scales = "free", space = "free") + theme_bw() + theme(
  text = element_text(size = 18),      
  axis.title = element_text(size = 14), 
  axis.text = element_text(size = 12),  
  legend.title = element_text(size = 14), 
  legend.text = element_text(size = 12)) + scale_fill_manual(values = hex_colors) +xlab("Mix Ratio")+ylab("Substrate Utilization Shannon Diversity Index") #+geom_jitter(width=0.2, size=0.5)
p_shan

#end  
  
  
  
  