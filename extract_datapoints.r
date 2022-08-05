# creates a data.frame that mimics the results of a typical RNAseq rifampicin timeseries based on a stochastic simulation

#probe_all: result list from the simulation script. Can be also a list of lists for replicate simulations.
#rif_time: time of simulated rifampicin addition
#timesteps: total steps of the simulation
#time_points: sample timepoints in seconds


extract_datapoints<-function(probe_all, rif_time,timesteps, time_points=c(0,1,2,3,4,6,8,10,15,20)*60, reps=1){
	if(reps==1){
		data_new<-c()
		for(j in 1:length(probe_all)){
			temp<-probe_all[[j]]
			temp<-temp[rif_time:timesteps]
			tmp<-temp[time_points+1]
			position<-rep(as.numeric(names(probe_all)[j]),length(tmp))
			out<-cbind((time_points)/60,tmp, position)
			colnames(out)<-c("time","inty","position")
			data_new<-rbind(data_new,out)
		}
		data_new<-as.data.frame(data_new)
	} else{
		data_new<-list()
		for(i in 1:length(probe_all[[1]])){
			Data5<-c()
			for(j in 1:length(probe_all)){
				 temp<-probe_all[[j]][[i]]
				 temp<-temp[rif_time:timesteps]
				 tmp<-temp[time_points+1]
				 position<-rep(as.numeric(names(probe_all[[1]])[1]),length(tmp))
				 out<-cbind((time_points)/60,tmp, position)
				 colnames(out)<-c("time","inty","position")
				 Data5<-rbind(Data5,out)
				  
			}
			data_new[[i]]<-as.data.frame(Data5)
		}
	}
	data_new
}


