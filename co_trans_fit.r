fit<-function(data_new){
	out_all<-c()
	for(j in 1:length(data_new)){
	  Data_fit<-as.data.frame(data_new[[j]])
	  colnames(Data_fit)<-c("time","inty","position")
	  tmp<-NA
	  
	  # start values 
	  decay=seq(.08, 0.11, by=.02)
	  delay=seq(0, 10, by=.1)
	  k=seq(0.1, 1, 0.2)

	  st1 <- expand.grid(decay=decay, delay=delay, k=k)
	  tryCatch(
		tmp <- nls2(
		  inty ~ I(time < delay)* k/decay +
			(time >= delay)*I(k/decay*(exp(-decay*(time - delay)))),
		  data = Data_fit,
		  algorithm = "port",
		  control = list(warnOnly = T),
		  start = st1,
		  lower = list(decay=0.01, delay=0.001)#
		)
		,error = function(e) {1+1}
	  )
	  if(is.na(tmp)==F){
		out_all<-rbind(out_all,coef(tmp))
	  } else {
		out_all<-rbind(out_all, rep(NA, 3))
	  }
	 }
	return(out_all)
}
