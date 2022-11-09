
# Script to simulate transcription and RNA degradation
# The script simulates the results from a transcriptomics experiment (Microarray or RNA-seq) for single transcripts
# The script mimics the actual cellular processes in a non-model based way

#Arguments:

#timesteps: Duration of simulation in seconds
#rif_time: Timepoint were initiation of synthesis stops
#start_pos: transcriptional start site(s) 
#rna_length: Length of the simulated transcriptional unit [nucleotides]
#probe_pos: Positions in the transcript were the RNA number is stored for each timestep
#pol_freq: initiation rate(s) [Initiations/s] 
#pol_freq2: Secondary initiation rate(s), for a change of initiation rate(s) after a given number of steps
#ini_rate_change_time: Timepoint were the Initiation rate switches from pol_freq to pol_freq2
#deg: Initial degradation constant [1/s]
#deg2: Secondary degradation constant, for a change of the degradation constant after a given number of steps 
#deg_change_time: Timepoint were the degradation constant switches from deg to deg2
#mode_of_decay: Switch between "Co-transcriptional degradation" (co) and "Post-transcriptional degradation" (post).
#pol_speed: Elongation rate [nt/s]
#pol_speed2: Secondary elongation rate for a "speed change" at a defined position.
#pol_speed2_position: Position for the start of the secondary elongation rate 


#term_length: Position for partial transcription termination of the transcript [nt]
#term_prob: Probability for transcription termination at term_length 
#end_of_termination_time: Timepoint were the the termination at term_length stops. Assuming a termination by another RNA polymerase, i.e. transcriptional interference, end_of_termination_time ~ rif_time.
#random_term_prob: Constant probability at each timestep for termination

#Transcriptional interference
#ti_anti_tss: start of antisense transcription in nt distance from sense start
#ti_anti_pol_freq: initiation rate(s) of asRNA(s)
#ti_anti_deg: degradation constant for the asRNA(s)
#ti_prob_sense: probability of termination for the sense transcript upon collision with an asRNA RNAP. The asRNA termination probability is 1-ti_prob_sense
#ti_range: Range in which RNAPs collide which each other [nt]
#ti_anti_usage: switch for explixit numeric modelling of ti. If ti_anti_usage=TRUE, term_prob and term_length should not be used.

simulate<-function(			timesteps=8000,
					rif_time=4000,
					start_pos=c(1,1500),
					rna_length=2000,
					pol_freq=c(0.1,0.01),  # [initiations/s]
					pol_freq2=pol_freq,
					deg=0.01,
					deg2=deg,
					pol_speed=10,
					ini_rate_change_time=timesteps,
					deg_change_time=timesteps,
					term_prob=0,
					term_length=rna_length,
					end_of_termination_time=rif_time,
					random_term_prob=0,
					probe_pos=c(50,350,500,1000,1500,1900),
					mode_of_decay="co",
					pol_speed2 = NA,
					pol_speed2_position = rna_length,
					
					# pausing sites
					pausing_position = 500,
					pausing_probability = 0,
					pausing_off_probability = 0.02, #[1/s]
					
					
					#ti antisense transcription parameters
					ti_anti_usage=FALSE,
					ti_anti_tss=0, 
					ti_anti_pol_freq=pol_freq[1]/5,
					ti_anti_deg=0.01,
					ti_prob_sense=0.1,
					#ti_prob_anti=1-ti_prob_sense 
					ti_rna_length=ti_anti_tss,
					ti_pol_speed=pol_speed,
					ti_range=30
					){
	
if(ti_anti_usage){  #either modelling of ti events by a factor or explixit numerical modelling of ti
	term_prob=0
	term_length=rna_length
}			

rna<-list()
probe<-vector("list",length(probe_pos))
names(probe)<-probe_pos

if(ti_anti_usage){
	probe_as<-vector("list",length(probe_pos))
	names(probe_as)<-probe_pos
	ti<-list()
	full_ti<-list()
}

if(pausing_probability>0){
	pause_list<-list()
}

count<-function(li, pos){
	cou<-length(which(unlist(li)==pos))
	cou
}

ladd<-function(x, nt=1){ #nt = number of nt added each step
	x<-c(x,(max(x)+1):(max(x)+nt))
	x
}

ladd_ti<-function(x, nt=1){ #nt = number of nt added each step
	x<-c(x,(min(x)-1):((min(x)-1)-nt))
	x
}


decay_fun<-function(x, nt=500){
	if(length(x)>0){
		nt<-min(length(x),nt)
		x<-x[nt:length(x)]
	}
	x
}

out<-c()
decay_rnas<-list()
fulllength_rna<-list()

for(i in 1:timesteps){

	# degradation
	if(i>deg_change_time){
			deg<-deg2
	}
	if(mode_of_decay=="co"){
		if(length(fulllength_rna)>0){
			random<-runif(n=length(fulllength_rna),min=0,max=1)
			random<-which(random<=deg)
			if(length(random)>0){
				 fulllength_rna<-fulllength_rna[-random]
			 }
		}
		if(length(rna)>0){
			random<-runif(n=length(rna),min=0,max=1)
			random<-which(random<=deg)
			if(length(random)>0){
				for(j in 1:length(random)){
					rna[[random[j]]]<-decay_fun(rna[[random[j]]],nt=rna_length)
				}
			}
		}
		if(pausing_probability>0){
			if(length(pause_list)>0){
				random<-runif(n=length(pause_list),min=0,max=1)
				random<-which(random<=deg)
				if(length(random)>0){
					for(j in 1:length(random)){
						pause_list[[random[j]]]<-decay_fun(pause_list[[random[j]]],nt=rna_length)
					}
				}
			}		
		}
		if(ti_anti_usage){
			if(length(ti)>0){
				random<-runif(n=length(ti),min=0,max=1)
				random<-which(random<=ti_anti_deg)
				if(length(random)>0){
					for(j in 1:length(random)){
						ti[[random[j]]]<-decay_fun(ti[[random[j]]],nt=rna_length)
					}
				}
			}
			if(length(full_ti)>0){
				 random<-runif(n=length(full_ti),min=0,max=1)
				 random<-which(random<=ti_anti_deg)
				 if(length(random)>0){
					 full_ti<-full_ti[-random]
				 }
			}			
		}		
	}
	
	
	if(mode_of_decay=="post"){
		if(length(fulllength_rna)>0){
			 random<-runif(n=length(fulllength_rna),min=0,max=1)
			 random<-which(random<=deg)
			 if(length(random)>0){
				 fulllength_rna<-fulllength_rna[-random]
			 }
		}
		if(ti_anti_usage){
			if(length(full_ti)>0){
				 random<-runif(n=length(full_ti),min=0,max=1)
				 random<-which(random<=ti_anti_deg)
				 if(length(random)>0){
					 full_ti<-full_ti[-random]
				 }
			}
		}
	}
	
	
	
	
	
	
	# termination 
	rna_l<-unlist(lapply(rna,max))
	fulll<-which(rna_l>=rna_length)
	if(length(fulll)>0){
		fulllength_rna<-c(fulllength_rna,rna[fulll]) # fullength rnas were transfered from the growing rna list
		rna<-rna[-fulll]
	}
	
	if(ti_anti_usage){
		ti_rna_l<-unlist(lapply(ti,min))
		ti_fulll<-which(ti_rna_l<=ti_rna_length)
		if(length(ti_fulll)>0){
			full_ti<-c(full_ti,ti[ti_fulll]) # fullength rnas were transfered from the growing rna list
			ti<-ti[-ti_fulll]  				 # to the fullength list
		}				
	}
	
	# random termination events
	if(random_term_prob>0){
		random<-runif(n=length(rna),min=0,max=1)
		terminated<-which(random<=random_term_prob)
		if(length(terminated)>0){
			fulllength_rna<-c(fulllength_rna,rna[terminated]) # pre-maturely terminated RNAs were transfered to the fullength list
			rna<-rna[-terminated]
		}
	}
	
	if(ti_anti_usage){
		random<-runif(n=length(ti),min=0,max=1)
		terminated<-which(random<=random_term_prob)
		if(length(terminated)>0){
			full_ti<-c(full_ti,ti[terminated]) # pre-maturely terminated RNAs were transfered to the fullength list
			ti<-ti[-terminated]
			
		}
	
	}
	
	#termination with a given factor at a defined position
	if(i<end_of_termination_time){  
		full2<-which(rna_l>=term_length-pol_speed/2 & rna_l<=term_length+pol_speed/2) # transcript in range of termination site
		if(length(full2)>0){			
			random<-runif(n=length(full2),min=0,max=1)
			random<-which(random<=term_prob)
			if(length(random)>0){
				fulllength_rna<-c(fulllength_rna,rna[full2[random]])
				rna<-rna[-full2[random]]
			}			
		}		
	}

	# termination by ti collision mechanisms, numerical simulation
	if(ti_anti_usage){
		rna_l<-unlist(lapply(rna,max)) # polymerase position of sense RNAs
		ti_rna_l<-unlist(lapply(ti,min)) # polymerase position of antisense RNAs
		ov_sense<-c()
		ov_anti<-c()
		# find sense/antisense pairs in collision distance
		for(jj in 1:length(ti_rna_l)){
			ov12<-which(ti_range-(ti_rna_l[jj]-rna_l)>=0 & ti_range-(ti_rna_l[jj]-rna_l)<mean(pol_speed,ti_pol_speed)*2)
			ov12<-na.omit(setdiff(ov12,ov_sense))
			if(length(ov12)>0){
				ov_sense<-c(ov_sense,ov12[1])
				ov_anti<-c(ov_anti,jj)
			}
		}
		na<-which(is.na(ov_sense))
		if(length(na)>0){
			break
		}
		ov_sense<-na.omit(ov_sense)
		if(length(ov_sense)>0){
			random<-runif(n=length(ov_sense),min=0,max=1)
			terminated<-which(random<=ti_prob_sense)			
			if(length(terminated)>0){
				fulllength_rna<-c(fulllength_rna,rna[ov_sense[terminated]])
				rna<-rna[-ov_sense[terminated]]
				ov_sense<-ov_sense[-terminated]	
			}
			
			if(length(terminated)>0){
				ov_anti<-ov_anti[-terminated]
			} 
	
	
			if(length(ov_anti)>0){
				full_ti<-c(full_ti,ti[ov_anti]) 
				ti<-ti[-ov_anti]  		
			}			
		}		
	}
	
	
	# Inititation of transcription
	for(jj in 1:length(start_pos)){
		pol_freq3<-pol_freq[jj]
		if(i>ini_rate_change_time){
			pol_freq3<-pol_freq2[jj]
		}
		repeats<-1
		if(pol_freq3>1){
			repeats<-as.numeric(paste(c(1,rep(0, nchar(as.integer(pol_freq3)))), collapse=""))
			pol_freq3<-pol_freq3/as.numeric(paste(c(1,rep(0, nchar(as.integer(pol_freq3)))), collapse=""))
		}
		if(i<rif_time){
			random<-runif(n=repeats,min=0,max=1)
			le<-length(which(random<=pol_freq3))
			if(le>0){
				new_rna<-rep(list(start_pos[jj]),le)
				rna<-c(rna,new_rna)
			}
		}
	}

	if(ti_anti_usage){
		# initiation of antisense transcript
		for(jj in 1:length(ti_anti_tss)){
			pol_freq3<-ti_anti_pol_freq[jj]
			repeats<-1
			if(pol_freq3>1){
				repeats<-as.numeric(paste(c(1,rep(0, nchar(as.integer(pol_freq3)))), collapse=""))
				pol_freq3<-pol_freq3/as.numeric(paste(c(1,rep(0, nchar(as.integer(pol_freq3)))), collapse=""))
			}
			if(i<rif_time){
				random<-runif(n=repeats,min=0,max=1)
				le<-length(which(random<=pol_freq3))
				if(le>0){
					new_rna<-rep(list(ti_anti_tss[jj]),le)
					ti<-c(ti,new_rna)
				}
			}
		}
	}

	# transcription elongation
	if(is.na(pol_speed2)){
		rna<-lapply(rna,ladd, nt=pol_speed) # extend existing RNA molecules by nt
	
	} else {
		rna_l<-unlist(lapply(rna,max))
		tmp<-which(rna_l>=pol_speed2_position)
		if(length(tmp)>0){
			if(length(tmp)<length(rna)){
				tmp2<-setdiff(1:length(rna),tmp)
				rna1<-lapply(rna[tmp],ladd, nt=pol_speed2)
				rna2<-lapply(rna[tmp2],ladd, nt=pol_speed)
				rna[tmp]<-rna1
				rna[tmp2]<-rna2
				
			} else {
				rna<-lapply(rna,ladd, nt=pol_speed2)
			}

		} else {
			rna<-lapply(rna,ladd, nt=pol_speed) # extend existing RNA molecules by nt
		}
	}
	

	if(ti_anti_usage){ # antisense transcript can only have an unchanged elongation rate
		ti<-lapply(ti,ladd_ti, nt=ti_pol_speed) # extend existing RNA molecules by nt	
	}
	
	
	
	# start of pause
	if(pausing_probability>0){
		rna_l<-unlist(lapply(rna,max))
		paused<-which(rna_l>pausing_position-pol_speed/2 & rna_l<pausing_position+pol_speed/2) # transcript in range of pausing site
		if(length(paused)>0){
			random<-runif(n=length(paused),min=0,max=1)
			random<-which(random<=pausing_probability)
			#print(random)
			if(length(random)>0){
				pause_list<-c(pause_list,rna[paused[random]])
				rna<-rna[-paused[random]]
			}	
		}
	}				
		
	# end of pause 
	if(pausing_probability>0){
		if(length(pause_list)>0){
			random<-runif(n=length(pause_list),min=0,max=1)
			ends<-which(random<=pausing_off_probability)
			if(length(ends)>0){
				rna<-c(rna,pause_list[ends])
				pause_list<-pause_list[-ends]
			}
		}
	}

	
	# counting transcript numbers overlapping the probe/bin positions
	for(j in 1:length(probe)){
		tmp<-count(rna,as.numeric(names(probe)[j]))
		tmp_full<-count(fulllength_rna,as.numeric(names(probe)[j]))
		#tmp_decay<-count(decay_rnas,as.numeric(names(probe)[j]))
		
		tmp<-tmp+tmp_full#+tmp_decay
		if(pausing_probability>0){
			if(length(pause_list)>0){
				tmp_pause<-count(pause_list,as.numeric(names(probe)[j]))
				tmp<-tmp+tmp_pause
			}
		}
		probe[[j]]<-c(probe[[j]],tmp)
		if(ti_anti_usage){
			tmp<-count(ti,as.numeric(names(probe_as)[j]))
			tmp_full<-count(full_ti,as.numeric(names(probe_as)[j]))
			#tmp_decay<-count(decay_rnas,as.numeric(names(probe)[j]))
			tmp<-tmp+tmp_full#+tmp_decay
			probe_as[[j]]<-c(probe_as[[j]],tmp)
			
		}		
	}
}
if(ti_anti_usage){
	return(list(probe,probe_as))
}
return(probe)
}
