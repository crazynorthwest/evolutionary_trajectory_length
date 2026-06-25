#Function to run Slim with specified code
RunSlim <- function(script)
{
  system2(path_Slim, args = shQuote(script), stdout=T, stderr=T)
}


lines_0 <- readLines(path_simulation)
for (idx_1 in 1:repeat_times)
{
  Optimal_list <- Optimal[,,idx_1]
  if (file.exists(path_simulation_tem)){file.remove(path_simulation_tem)}
  lines_1 <- paste0('function (string)addPath(void){path_data = "',paste0(path, '/reused_', file_used_number),'";')
  lines_2 <- 'return path_data;}'
  lines_3 <- paste0('function (void)addParameter(void){if(!exists("n")){defineConstant("n",',trait_number,');}')
  lines_4 <- paste0('if(!exists("filenamenumber")){defineConstant("filenamenumber",',idx_1,');}')
  lines_5 <- paste0('if(!exists("K")){defineConstant("K",',K,');}')
  lines_6 <- paste0('if(!exists("q")){defineConstant("q",',qq[idx_1],');}')
  lines_8 <- paste0('if(!exists("m")){defineConstant("m",',m,');}')
  lines_9 <- paste0('if(!exists("genome_size")){defineConstant("genome_size",',genome_size,');}')
  lines_10 <- paste0('if(!exists("lambda")){defineConstant("lambda",',lambda,');}')
  lines_11 <- paste0('if(!exists("Optimal_list")){defineConstant("Optimal_list",matrix(c(',paste(c(Optimal_list),collapse = ", "),'),nrow=',generation_total,'));}}')
  lines_final <- c(lines_1, lines_2, lines_3, lines_4, lines_5, lines_6, lines_8, lines_9, lines_10, lines_11, lines_0)
  writeLines(lines_final, con = path_simulation_tem, sep = "\n", useBytes = F)
  RunSlim(path_simulation_tem)
}

