
clusteroneR = function(network, pp=8, density_threshold=0, 
                       java_path="../java/cluster_one-1.0.jar",
                       output_path="../data/") {
  
  # 1. write network to file
  # ensure network is symmetric, i.e. contains A-B and B-A
  I = !duplicated(paste(network[,1], network[,2], sep="-"))
  network = network[I, ]
  network2 = network[,c(2,1)]
  names(network2) = names(network)
  network = rbind(network, network2)
  fn_tmpM = paste(output_path,round(runif(1)*1e8),".network.txt", sep="")
  print(fn_tmpM)
  write_tsv(network, fn_tmpM)
  
  # 2. Call clusterone java
  fn_tmpout = paste(output_path,round(runif(1)*1e8),".clusterone.txt", sep="")
  system_call = paste('java -jar ', java_path, ' ', fn_tmpM, ' -s 2 -d ', density_threshold,
                 ' --penalty ', pp, ' > ', fn_tmpout)
  system(system_call)
  
  # 3. Read clusterone java output
  co.results = as.data.frame(read_csv(fn_tmpout))
  co.results = as.list(sapply(co.results, FUN=function(x) gsub("\t", ";", x)))
  
  # 4. delete tmp files
  if (file.exists(fn_tmpM)) file.remove(fn_tmpM)
  if (file.exists(fn_tmpout)) file.remove(fn_tmpout)
  
  return(co.results)
}

