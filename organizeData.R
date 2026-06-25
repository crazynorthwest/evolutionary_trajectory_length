fileName_phe_record <- paste0(path_analysis,"/data_phe1.dat")
fileName_mut_record <- paste0(path_analysis,"/data_mut1.dat")

data_record_phe <- cbind("repeat_number",  "generation", rbind(paste0(rep("mean_phe", trait_number),(1:trait_number))), rbind(paste0(rep("var_phe", trait_number),(1:trait_number))), "fitness_mean", "fitness_var","population_id", "FST")
data_record_mut <- cbind("repeat_number",  "mut_size", "generation", "frequence")
write.table(data_record_phe, file = fileName_phe_record, row.names = FALSE, col.names = FALSE)
write.table(data_record_mut, file = fileName_mut_record, row.names = FALSE, col.names = FALSE)
fitness_idx <- 0.95
generation_record <- NULL
average_mut_size <- NULL
fix_mut_number <- NULL
mut_fix_sum <- NULL
for (idx_1 in 1:repeat_times){
# for (idx_1 in 501:1000){
  fileName_phe <- paste0(path_result,"/data_phe",idx_1,".dat")
  fileName_mut1 <- paste0(path_result,"/data_p1_mut",idx_1,".dat")
  fileName_mut2 <- paste0(path_result,"/data_p2_mut",idx_1,".dat")
  fileName_mut_size <- paste0(path_mutation,"/",idx_1,".dat")
  data_phe <- read.csv(fileName_phe, header = TRUE, sep = " ")
  data_mut1 <- read.csv(fileName_mut1, header = TRUE, sep = " ")
  data_mut2 <- read.csv(fileName_mut2, header = TRUE, sep = " ")
  data_mut_size <- read.csv(fileName_mut_size, header = TRUE, sep = " ")
  data_record_phe <- cbind(matrix(rep(idx_1, length(data_phe$generation)), nrow = length(data_phe$generation)), data_phe)
  data_record_mut_sum <- NULL
  for (idx_2 in 1:genome_size) {
    data_record_mut1 <- cbind(matrix(rep(c(idx_1, data_mut_size[idx_2, 2], "p1"), length(data_mut1$generation)), nrow = length(data_mut1$generation), byrow = TRUE), unname(data_mut1[,c(1,idx_2+1)]))
    data_record_mut2 <- cbind(matrix(rep(c(idx_1, data_mut_size[idx_2, 2], "p2"), length(data_mut2$generation)), nrow = length(data_mut2$generation), byrow = TRUE), unname(data_mut2[,c(1,idx_2+1)]))
    data_record_mut_sum <- rbind(data_record_mut_sum, data_record_mut1, data_record_mut2)
    }
  data_record_mut_sum <- as.data.frame(data_record_mut_sum)
  colnames(data_record_mut_sum) <- c("repeat_number",  "mut_size", "population_id","generation", "frequence")
  data_record_mut_sum$mut_size <- as.numeric(data_record_mut_sum$mut_size)

  data_phe1 <- subset(data_phe, population_id == "p1")
  generation_idx1 <- data_phe1$generation[which(data_phe1$fitness_mean > fitness_idx)[1]]
  fitness_final1 <- data_phe1$fitness_mean[data_phe1$generation==4998 & data_phe1$population_id == "p1"]
  fitness_f1 <- data_phe$fitness_mean[data_phe$generation==4999 & data_phe$population_id == "p3"]
  fitness_f2 <- data_phe$fitness_mean[data_phe$generation==5000 & data_phe$population_id == "p3"]
  generation_record <- rbind(generation_record, cbind(idx_1, "p1", generation_idx1, fitness_final1, fitness_f1, fitness_f2))
  data_phe2 <- subset(data_phe, population_id == "p2")
  generation_idx2 <- data_phe2$generation[which(data_phe2$fitness_mean > fitness_idx)[1]]
  fitness_final2 <- data_phe2$fitness_mean[data_phe2$generation==4998 & data_phe2$population_id == "p2"]
  generation_record <- rbind(generation_record, c(idx_1, "p2", generation_idx2, fitness_final2, fitness_f1, fitness_f2))

  data_mut_sub1 <- subset(data_record_mut_sum,  generation ==  generation_idx1 & population_id =="p1")
  fix_mut_number <- c(fix_mut_number, sum(data_mut_sub1$frequence > 0.95))
  data_mut_sub2 <- subset(data_record_mut_sum,  generation ==  generation_idx2 & population_id =="p2")
  fix_mut_number <- c(fix_mut_number, sum(data_mut_sub2$frequence > 0.95))
  average_mut_size <- c(average_mut_size, sum(data_mut_sub1$mut_size * data_mut_sub1$frequence))
  average_mut_size <- c(average_mut_size, sum(data_mut_sub2$mut_size * data_mut_sub2$frequence))
  mut_fix1 <- data_mut_sub1$mut_size[data_mut_sub1$frequence > 0.95]
  mut_fix2 <- data_mut_sub2$mut_size[data_mut_sub2$frequence > 0.95]
  if (!is.na(generation_idx1) & length(mut_fix1) > 0){
    data_mut_sub1_2 <- subset(data_record_mut_sum, mut_size %in% mut_fix1 & frequence > 0.01 & population_id == "p1")
    appear_first1 <- aggregate( generation ~ mut_size, data = data_mut_sub1_2, FUN = min)
    colnames(appear_first1)[colnames(appear_first1) == "generation"] <- "appear_first"
    fix_first1 <- aggregate( generation ~ mut_size, data = data_mut_sub1_2[data_mut_sub1_2$frequence > 0.95,], FUN = min)
    appear_first1$fix_first <- fix_first1$generation
    appear_first1$fix_time <- appear_first1$fix_first - appear_first1$appear_first
    appear_first1$generation_record <- rep(generation_idx1, sum(data_mut_sub1$frequence > 0.95))
    appear_first1$par <- qq[idx_1]
    appear_first1$population_id <- "q1"
    mut_fix_sum <- rbind(mut_fix_sum, appear_first1)
  }
  if (!is.na(generation_idx2) & length(mut_fix2) > 0){
    data_mut_sub2_2 <- subset(data_record_mut_sum, mut_size %in% mut_fix2 & frequence > 0.01 & population_id == "p2")
    appear_first2 <- aggregate( generation ~ mut_size, data = data_mut_sub2_2, FUN = min)
    colnames(appear_first2)[colnames(appear_first2) == "generation"] <- "appear_first"
    fix_first2 <- aggregate( generation ~ mut_size, data = data_mut_sub2_2[data_mut_sub2_2$frequence > 0.95,], FUN = min)
    appear_first2$fix_first <- fix_first2$generation
    appear_first2$fix_time <- appear_first2$fix_first - appear_first2$appear_first
    appear_first2$generation_record <- rep(generation_idx2, sum(data_mut_sub2$frequence > 0.95))
    appear_first2$par <- qq[idx_1]
    appear_first2$population_id <- "q2"
    mut_fix_sum <- rbind(mut_fix_sum, appear_first2)
  }
}

generation_record <- as.data.frame(generation_record)
colnames(generation_record) <- c("repeat_number", "population_id", "generation_record", "fitness_final", "fitness_f1", "fitness_f2")
generation_record$par <- rep(qq, each = 2)
generation_record$par_str <-  paste0("q=", generation_record$par)

generation_record$mut_size_ave <- average_mut_size
generation_record$fix_mut_number <- fix_mut_number
mut_fix_sum$par_str <-  paste0("q=", mut_fix_sum$par)
write.table(generation_record, file = paste0(path_figures,"/generation_record.dat"), row.names = FALSE, col.names = TRUE)
write.table(mut_fix_sum, file = paste0(path_figures,"/mut_fix_sum.dat"), row.names = FALSE, col.names = TRUE)



fileName_test <- paste0(path_figures,"/power_law_test1.csv")
data_output <- as.data.frame(generation_record$generation_record[!is.na(generation_record$generation_record) & generation_record$par == 0.5])
colnames(data_output) <- "Steps"
data_output$Steps <- as.integer(data_output$Steps)
write.table(data_output, file = fileName_test, row.names = FALSE, col.names = TRUE)


fileName_test <- paste0(path_figures,"/power_law_test2.csv")
data_output <- as.data.frame(generation_record$generation_record[!is.na(generation_record$generation_record) & generation_record$par == 0.75])
colnames(data_output) <- "Steps"
data_output$Steps <- as.integer(data_output$Steps)
write.table(data_output, file = fileName_test, row.names = FALSE, col.names = TRUE)


fileName_test <- paste0(path_figures,"/power_law_test3.csv")
data_output <- as.data.frame(generation_record$generation_record[!is.na(generation_record$generation_record) & generation_record$par == 1])
colnames(data_output) <- "Steps"
data_output$Steps <- as.integer(data_output$Steps)
write.table(data_output, file = fileName_test, row.names = FALSE, col.names = TRUE)


fileName_test <- paste0(path_figures,"/power_law_test4.csv")
data_output <- as.data.frame(generation_record$generation_record[!is.na(generation_record$generation_record) & generation_record$par == 1.25])
colnames(data_output) <- "Steps"
data_output$Steps <- as.integer(data_output$Steps)
write.table(data_output, file = fileName_test, row.names = FALSE, col.names = TRUE)


fileName_test <- paste0(path_figures,"/power_law_test5.csv")
data_output <- as.data.frame(generation_record$generation_record[!is.na(generation_record$generation_record) & generation_record$par == 1.5])
colnames(data_output) <- "Steps"
data_output$Steps <- as.integer(data_output$Steps)
write.table(data_output, file = fileName_test, row.names = FALSE, col.names = TRUE)
