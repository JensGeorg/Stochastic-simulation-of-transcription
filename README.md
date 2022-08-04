# Stochastic-simulation-of-transcription

This R function simulates the bacterial RNA transcription and decay. It is predominantly designed to generate time series data after a stop of the transcription initiation, which resembles the effect of the antibiotic Rifampicin. Experimental Rifampicin time series data are widely used to study the RNA stability. 
Various aspects of transcription can be investigated.

## Transcription and degradation of a 1000nt RNA (co-transcriptional decay)

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
<p align="center">
  <img src="https://github.com/JensGeorg/Stochastic-simulation-of-transcription/blob/main/simulate_figs/co_transcriptional_1000nt.png"/>
</p>
