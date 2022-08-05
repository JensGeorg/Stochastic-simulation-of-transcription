# Stochastic simulation of transcription

This R function simulates the bacterial RNA transcription and decay. It is predominantly designed to generate time series data after a stop of the transcription initiation, which resembles the effect of the antibiotic Rifampicin. Experimental Rifampicin time series data are widely used to study the RNA stability. 
Various aspects of transcription can be investigated.

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
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(0, total_time), xlab="molecules", ylab="time [s]", main="co-transcriptional decay")
for(i in 1:length(dat)){
  points(1:total_time,dat[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)


# visualization of the decay curves with analytical solution 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(rif_time, total_time), xlab="molecules", ylab="time [s]", main="co-transcriptional decay")
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
<p float="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/co_transcriptional_1000nt.png" width="400"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/co_transcriptional_1000nt_with_analyt.png" width="400"/>
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
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(0, total_time), xlab="molecules", ylab="time [s]", main="post-transcriptional decay")
for(i in 1:length(dat)){
  points(1:total_time,dat[[i]], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)


# visualization of the decay curves with analytical solution 
plot(1,1,type="n", ylim=c(0,max(unlist(dat))), xlim=c(rif_time, total_time), xlab="molecules", ylab="time [s]", main="post-transcriptional decay")
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
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/post_trans.png" width="500"/>
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/post_trans_with_analyt.png" width="500"/>
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
pos=c(1,seq(200,3000,500),seq(1900,1990,10)) #sample positions [nt]
pol_speed=10 #[nt/s]
rna_length=3000 #length of transcript [nt]

rif_time=steady_state*10  # time when transcription initiation stops
total_time=rif_time + 8 * log(2)/deg 

dat<-simulate(timesteps=total_time,
              rif_time=rif_time,
              pol_freq=pol_freq,
              ti_anti_pol_freq=ti_anti_pol_freq,
              ti_anti_deg=ti_anti_deg,
              deg=deg,
              ti_prob_sense=ti_prob_sense,
              start_pos=start_pos,
              probe_pos=pos,
              pol_speed=pol_speed,
              rna_length=rna_length,
              ti_anti_usage=TRUE,
              mode_of_decay="co")

dat_sense<-dat[[1]]
dat_anti<-dat[[2]]

# visualization of the decay curves 
plot(1,1,type="n", ylim=c(0,max(unlist(dat_sense))), xlim=c(rif_time, total_time), xlab="molecules", ylab="time [s]", main="collision TI (sense transcript)")
for(i in 1:length(dat_sense)){
  points(rif_time:total_time,dat_sense[[i]][rif_time:total_time], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)

# visualization of the decay curves 
plot(1,1,type="n", ylim=c(0,max(unlist(dat_anti))), xlim=c(rif_time, total_time), xlab="molecules", ylab="time [s]", main="collision TI (sense transcript)")
for(i in 1:length(dat_anti)){
  points(rif_time:total_time,dat_anti[[i]][rif_time:total_time], type="l", col=i)
}
legend("topright", bty="n", legend=pos, text.col=1:length(pos))
abline(v=rif_time)
```
