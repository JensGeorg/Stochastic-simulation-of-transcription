# Stochastic simulation of transcription

This R function simulates the bacterial RNA transcription and decay. It is predominantly designed to generate time series data after a stop of the transcription initiation, which resembles the effect of the antibiotic Rifampicin. Experimental Rifampicin time series data are widely used to study the RNA stability. Some use cases and the respective visualizations are shown in the following. Due to the stochastic nature, the data and visualizations will change from run to run.

- [Transcription and degradation of a 1000nt RNA (co-transcriptional decay)](#Transcription-and-degradation-of-a-1000nt-RNA-(co-transcriptional-decay))
- [Transcription and degradation of a 1000nt RNA (post-transcriptional decay)](#Transcription-and-degradation-of-a-1000nt-RNA-(post-transcriptional-decay))
- [Transcriptional interference by the collision of sense and antisense transcription](#Transcriptional-interference-by-the-collision-of-sense-and-antisense-transcription)
- [Simulation of a transcriptional pausing site](#Simulation-of-a-transcriptional-pausing-site)
  - [Fitting of decay curves](#Fitting-of-decay-curves)  
  - [Segmentation of the delay coefficients](#Segmentation-of-the-delay-coefficients)

### Transcription and degradation of a 1000nt RNA (co-transcriptional decay)

```
source("simulate.r")
pol_freq=0.6 #[initiations/s]
deg=0.01 #[1/s]
start_pos=1
pos=c(1,seq(200,1000,200)) #sample positions [nt]
pol_speed=10 #[nt/s]
rna_length=1000 #length of transcript [nt]

steady_state=as.integer(log(0.005)/-deg - rna_length/pol_speed) #estimate when steady state is reached for the full-length transcript
## rna_length/pol_speed = delay for the 3'end of the transcript

rif_time=steady_state*1.5  # time when transcription initiation stops
total_time=rif_time + 6 * log(2)/deg # the simualtion stops 6 half-life times after the stop of transcription initiation

dat<-simulate(timesteps=total_time,
              rif_time=rif_time,
              pol_freq=pol_freq,
              deg=deg,
              start_pos=start_pos,
              probe_pos=pos,
              pol_speed=pol_speed,
              rna_length=1000,
              mode_of_decay="co")
```
The counted transcript numbers at the respective positions are stored in a list. Each list entry corresponds to a position and contains a vector with the transcript numbers at each timestep. These data can be used to draw simulated synthesis/decay curves.

```
# visualization of the synthesis and decay phase 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(0, total_time), xlab="time [s]", ylab="molecules", main="co-transcriptional decay")
for(i in 1:length(dat)){
  points(1:total_time,dat[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)


# visualization of the decay curves with analytical solution 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(rif_time, total_time), xlab="time [s]", ylab="molecules", main="co-transcriptional decay")
for(i in 1:length(dat)){
  points(rif_time:total_time,dat[[i]][rif_time:total_time], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)

for(i in 1:length(pos)){
  curve(pol_freq/deg *x/x, from=0 , to= pos[i]/pol_speed+rif_time, add=TRUE, col=i, lty=2)
  curve(pol_freq/deg * exp(-deg*(x-(pos[i]/pol_speed+rif_time))), from=pos[i]/pol_speed+rif_time , to= total_time, add=TRUE, col=i, lty=2)
}
```
The analytical solution regardind the co-transcriptional decay model is shown as broken line.

<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/co_transcriptional_full.png" width="350"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/co_transcriptional_short.png" width="350"/>
</p>

### Transcription and degradation of a 1000nt RNA (post-transcriptional decay)
The script can be also used to simulate the hypothetical post-transcriptional decay

```
source("simulate.r")
pol_freq=0.6 #[initiations/s]
deg=0.01 #[1/s]
start_pos=1
pos=c(1,seq(200,1000,200)) #sample positions [nt]
pol_speed=10 #[nt/s]
rna_length=1000 #length of transcript [nt]

steady_state=as.integer(log(0.005)/-deg - rna_length/pol_speed) #estimate when steady state is reached for the full-length transcript
## rna_length/pol_speed = delay for the 3'end of the transcript

rif_time=steady_state*1.5  # time when transcription initiation stops
total_time=rif_time + 6 * log(2)/deg # the simualtion stops 6 half-life times after the stop of transcription initiation

dat<-simulate(timesteps=total_time,
              rif_time=rif_time,
              pol_freq=pol_freq,
              deg=deg,
              start_pos=start_pos,
              probe_pos=pos,
              pol_speed=pol_speed,
              rna_length=1000,
              mode_of_decay="post")
              
# visualization of the synthesis and decay phase 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(0, total_time), xlab="time [s]", ylab="molecules", main="post-transcriptional decay")
for(i in 1:length(dat)){
  points(1:total_time,dat[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)


# visualization of the decay curves with analytical solution 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(rif_time, total_time), xlab="time [s]", ylab="molecules", main="post-transcriptional decay")
for(i in 1:length(dat)){
  points(rif_time:total_time,dat[[i]][rif_time:total_time], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)

for(i in 1:length(pos)){
  curve(pol_freq/deg*x/x+pol_freq*(max(pos)-pos[i])/(pol_speed),from=0,to=pos[i]/(pol_speed)+rif_time, add=T, type="l", col=i,lty=2)
  curve(pol_freq/deg+pol_freq*(max(pos)/(pol_speed)-(x-rif_time)),from=pos[i]/(pol_speed)+rif_time,to=max(pos)/(pol_speed)+rif_time, add=T,  col=i,lty=2)	
  curve(pol_freq/deg * exp(-deg*(x-(max(pos)/pol_speed+rif_time))), from=max(pos)/pol_speed+rif_time , to= total_time, add=T, col=i, lty=2)
}
```
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/post_transcriptional_full.png" width="350"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/post_transcriptional_short.png" width="350"/>
</p>

### Transcriptional interference by the collision of sense and antisense transcription
For this scenario the transcription of a 3000nt sense RNA and an asRNA starting at position 2000 in reverse orientation is simulated.

```
source("simulate.r")
pol_freq=0.7 #[initiations/s]
ti_anti_pol_freq=0.7 #[initiations/s]
deg=0.01 #[1/s]
ti_anti_tss=2000
ti_anti_deg=0.01 #[1/s]
ti_prob_sense=0.29
start_pos=1
ti_rna_length=1000
pos=c(1,seq(200,3000,500),seq(1900,1990,10)) #sample positions [nt]
pol_speed=10 #[nt/s]
rna_length=3000 #length of transcript [nt]

rif_time=steady_state * 2   # time when transcription initiation stops
total_time=rif_time + 10 * log(2)/deg 

dat<-simulate(timesteps=total_time,
              rif_time=rif_time,
              pol_freq=pol_freq,
              ti_anti_pol_freq=ti_anti_pol_freq,
              ti_anti_deg=ti_anti_deg,
              ti_anti_tss=ti_anti_tss,
              ti_rna_length=ti_rna_length,
              deg=deg,
              ti_prob_sense=ti_prob_sense,
              start_pos=start_pos,
              probe_pos=pos,
              pol_speed=pol_speed,
              rna_length=rna_length,
              ti_anti_usage=TRUE,
              mode_of_decay="co")

# visualization of the full synthesis and decay curves 
plot(1,1,type="n", ylim=c(0,max(unlist(dat_sense))), xlim=c(0, total_time), xlab="time [s]", ylab="molecules", main="collision TI (sense transcript)")
for(i in 1:length(dat_sense)){
  points(1:total_time,dat_sense[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)

plot(1,1,type="n", ylim=c(0,max(unlist(dat_anti))), xlim=c(0, total_time), xlab="time [s]", ylab="molecules", main="collision TI (antisense transcript)")
for(i in 1:length(dat_anti)){
  points(1:total_time,dat_anti[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)
```
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/collision_ti_sense_full.png" width="350"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/collision_ti_anti_full.png" width="350"/>
</p>

```
# visualization of the decay curves (after rifampicin addition)
dat_sense<-dat[[1]][c(1,4,6,7,8)]
dat_anti<-dat[[2]][c(17,16,15,14,13,12,11,10,9)]

ti<-ti_anti_pol_freq*ti_prob_sense
pos2<-as.numeric(names(dat_sense))

plot(1,1,type="n", ylim=c(0,max(unlist(dat_sense))), xlim=c(rif_time, total_time), xlab="time [s]", ylab="molecules", main="collision TI (sense transcript)")
for(i in 1:length(dat_sense)){
  points(rif_time:total_time,dat_sense[[i]][rif_time:total_time], type="l", col=i)
}
legend("topright", bty="n", legend=names(dat_sense), text.col=1:length(pos2))
abline(v=rif_time)

for(i in 1:length(dat_sense)){
  if(pos2[i]>=ti_anti_tss){
    curve((pol_freq/deg-ti/deg*x/x),from=0,to=(pos2[i]-ti_anti_tss)/pol_speed+rif_time, add=T,col=i,lty=2) 
    curve(pol_freq/deg - ti/deg * exp(-deg*(x-(pos2[i]-ti_anti_tss)/pol_speed-rif_time)),from=(pos2[i]-ti_anti_tss)/pol_speed+rif_time,to=pos2[i]/pol_speed+rif_time, add=T, col=i,lty=2) 
    curve( (pol_freq/deg  - ti/deg * exp(-deg*(ti_anti_tss/pol_speed)))*exp(-deg*(x-pos2[i]/pol_speed-rif_time)),from=pos2[i]/pol_speed+rif_time,to=total_time,lty=2, add=T, col=i)
    
  }

  if(pos2[i]<ti_anti_tss){
    curve(pol_freq/deg *x/x, from=0 , to= pos2[i]/pol_speed+rif_time, add=TRUE, col=i, lty=2)
    curve(pol_freq/deg * exp(-deg*(x-(pos2[i]/pol_speed+rif_time))), from=pos2[i]/pol_speed+rif_time , to= total_time, add=TRUE, col=i, lty=2)
  }
}
```
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/collision_sense_short.png" width="350"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/collision_ti_anti_short.png" width="350"/>
</p>

The positional RNA abundance drops after the start of the asRNA due to termination by transcriptional interference

```
dat_sense<-dat[[1]]
dat_anti<-dat[[2]]
plot(1,1,type="n", xlim=c(0,max(pos)), ylim=c(0, max(unlist(dat))), xlab="position [nt]",ylab="molecules", main="positional RNA abundance in steady state" )

for(i in 1:length(dat_sense)){
  points(as.numeric(names(dat_sense)[i]),dat_sense[[i]][rif_time], pch=19)
  points(as.numeric(names(dat_anti)[i]),dat_anti[[i]][rif_time], pch=19, col=2)
}
legend("topright", bty="n", legend=c("sense","anti"), pch=c(19,19), col=c(1,2), text.col=c(1,2))

```
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/ti_positional_steady_state_abundance.png" width="350"/>
 </p>
 
### Simulation of a transcriptional pausing site
The transcription of a 2000nt transcript is simulated. At position 1000 is a pausing site with an off rate of 0.02 1/s (mean lifetime 50s).

```
source("simulate.r")
pol_freq=0.6 #[initiations/s]
deg=0.01 #[1/s]
start_pos=1
pos=seq(5,2000,50) #sample positions [nt]
pol_speed=10 #[nt/s]
rna_length=1000 #length of transcript [nt]

steady_state=as.integer(log(0.005)/-deg - rna_length/pol_speed) #estimate when steady state is reached for the full-length transcript
## rna_length/pol_speed = delay for the 3'end of the transcript

rif_time=steady_state*1.5  # time when transcription initiation stops
total_time=rif_time + 10 * log(2)/deg # the simualtion stops 6 half-life times after the stop of transcription initiation

dat<-simulate(timesteps=total_time,
              rif_time=rif_time,
              pol_freq=pol_freq,
              deg=deg,
              start_pos=start_pos,
              probe_pos=pos,
              pol_speed=pol_speed,
              rna_length=2000,
              mode_of_decay="co",
              
              pausing_position = 1000,
              pausing_probability = 1,
					    pausing_off_probability = 0.02 
              )
              
dat2<-dat[seq(1,40,5)]
pos2<-as.numeric(names(dat2))
# visualization of the synthesis and decay phase 
plot(1,1,type="n", ylim=c(0,max(unlist(dat2))), xlim=c(0, total_time), xlab="time [s]", ylab="molecules", main="co-transcriptional decay (with pausing site)")
for(i in 1:length(dat2)){
  points(1:total_time,dat2[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=names(dat2), text.col=1:length(pos2))
abline(v=rif_time)


# visualization of the decay curves with analytical solution 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(rif_time, total_time), xlab="time [s]", ylab="molecules", main="co-transcriptional decay (with pausing site)")
for(i in 1:length(dat2)){
  points(rif_time:total_time,dat2[[i]][rif_time:total_time], type="l", col=i)
}

legend("topright", bty="n", legend=names(dat2), text.col=1:length(pos2))
abline(v=rif_time)    
              
```
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/pausing_long.png" width="350"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/pausing_short.png" width="350"/>
</p>


#### Fitting of decay curves  
#### Segmentation of the delay coefficients
