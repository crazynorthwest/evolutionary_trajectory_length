####################################################
##  
## 
## 12/02/2025
####################################################

#Please set repeat times to run the simulation.
repeat_each <- 1000
#Strength of selection
qqq <- c(0.5, 0.75, 1, 1.25, 1.5)
qq <- rep(qqq, repeat_each)
repeat_times <- length(qq)
#Please set the following path to the path where this R script is located.
path_code <- "C:/Users/jxhan/Desktop/unfinished/evolutionary path2/main.R"
#Please set the following path to Slim.
path_Slim <- "C:/msys64/mingw64/bin/slim"

#Trait number
trait_number <- 1
#Carrying capacity
K <- 500
#Degree of epistasis
m <- 2
#Genome size
genome_size <- 50
#Lambda for mutation effect size
lambda <- 0.1



generation_total <- 5000
#generation_total <- generation_1 + generation_2 + generation_3 + generation_4

Optimal <- array(dim = c(generation_total, trait_number, repeat_times))

scenario_final <-  matrix(c(rep(1, generation_total),rep(0, generation_total*(trait_number-1))), nrow = generation_total)
#scenario_final <- rbind(matrix(c(rep(1, generation_3),rep(0, generation_3*(trait_number-1))),nrow = generation_3), matrix(rep(0, generation_4*trait_number),nrow = generation_4))


idx_1 <- 1
for (idx_1 in 1:repeat_times) {
    Optimal[,,idx_1] <- scenario_final
}


#Build a vector of phenotype names.
phenotype_name <- c(paste0(rep("phenotype", trait_number),(1:trait_number)))
#Setting the main path.
path <- substr(path_code, 1, nchar(path_code)-7)

#Find the number of times it has been reused to avoid new runs overwriting previous data.
file_used_number <- 1
idx_tem <- 1
while (idx_tem > 0){
  fileName <- paste0(path, '/reused_', file_used_number)
  if (file.exists(fileName)){
    file_used_number <- file_used_number + 1
  }
  else{
    dir.create(fileName);
    idx_tem <- 0;
  }
}

#Setting other pathes.
path_figures <- paste0(path, '/reused_', file_used_number, '/figures')
path_result <- paste0(path, '/reused_', file_used_number, '/result')
path_analysis <- paste0(path, '/reused_', file_used_number, '/analysis')
path_mutation <- paste0(path, '/reused_', file_used_number, '/mutation')
if(!file.exists(path_figures)){dir.create(path_figures)}
if(!file.exists(path_result)){dir.create(path_result)}
if(!file.exists(path_analysis)){dir.create(path_analysis)}
if(!file.exists(path_mutation)){dir.create(path_mutation)}

path_simulation <- paste0(path, '/simulation.slim')
path_simulation_tem <- paste0(path, '/reused_', file_used_number,'/simulation_tem.slim')
path_runSlim <- paste0(path, '/runSlim.R')
path_organizeData <- paste0(path, '/organizeData.R')
path_analysisData <- paste0(path, '/analysisData.R')
path_selection <- paste0(path, '/model_selection.R')
path_mu_selection <- paste0(path, '/mu_selection_strength.R')
write.table(cbind("repeat_number", "phenotype"), file = paste0(path_result,"/data_ind.dat"), row.names = FALSE, col.names = FALSE)
# Run simulations.
source(path_runSlim)
#Re-organize data
source(path_organizeData)
#Analysis data
source(path_analysisData)
