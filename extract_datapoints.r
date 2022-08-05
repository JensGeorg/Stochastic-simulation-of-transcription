# creates a data.frame that mimics the results of a typical RNAseq rifampicin timeseries based on a stochastic simulation

#probe_all: result list from the simulation script. Can be also a list of lists for replicate simulations.
#rif_time: time of simulated rifampicin addition
#timesteps: total steps of the simulation
#time_points: sample timepoints in seconds

extract_datapoints<-function(probe_all, rif_time,timesteps, time_points=c(1,c(1,2,3,4,6,8,10,15,20)*60)){
	data_new<-list()
	for(i in 1:length(probe_all[[1]])){
		Data5<-c()
		for(j in 1:length(probe_all)){
			 temp<-probe_all[[j]][[i]]
			 temp<-temp[rif_time:timesteps]
			 tim<-time_points
			 tmp<-temp[tim]
			 position<-rep(as.numeric(names(probe_all[[1]])[1]),length(tmp))
			 group<-rep(1,length(tmp))
			 out<-cbind((tim-1)/60,tmp,group, position)
			 colnames(out)<-c("time","inty","group","position")
			 Data5<-rbind(Data5,out)
			  
		}
		data_new[[i]]<-as.data.frame(Data5)
	}
	data_new
}
